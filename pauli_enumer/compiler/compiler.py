from functions import *
import synthesis_SC
from synthesis_sd import *
import time
ctime = time.time
from qiskit import QuantumCircuit

def simple_qubit_mapping(np):
    pi = {}
    for i in range(np):
        pi[i] = i
    return pi

class Board:
    def __init__(self, graph, cs):
        self.graph = graph
        self.color = []  # 一种颜色占据的位置
        self.grid = []  # 一个位置的颜色
        self.edge = []  # 边缘位置
        # self.weight = [max(1/(1.1**l), 0.000000001) for l in layers]
        self.ct = -1
        md = 10000
        for p1 in range(len(graph.C)):
            d = 0
            for p2 in range(len(graph.C)):
                d += self.graph.C[p1][p2]
            if d < md:
                md = d
                self.ct = p1
        self.edge.append(self.ct)
        for i in range(len(graph.G)):
            self.grid.append([])
        for c in cs:
            self.color.append([])

    def move(self, c, pos):
        self.grid[pos] = c
        self.edge.remove(pos)
        for i in c:
            self.color[i].append(pos)
        for a in self.graph[pos].adj:
            if len(self.grid[a]) == 0 and a not in self.edge:
                self.edge.append(a)

    def score(self, c, p):
        td = 0
        for ci in c:
            mdci = 10000
            if len(self.color[ci]) == 0:
                mdci = 0
            for pci in self.color[ci]:
                mdci = min(mdci, self.graph.C[pci][p])
            td += mdci  # * self.weight[ci]
        td = td * 100 + self.graph.C[self.ct][p]
        return td

    def unconnected_degree(self, pss, pi):
        ud = 0
        for ps in pss:
            connected = ps[:1]
            remain = ps[1:]
            for i in range(len(ps) - 1):
                md = 10000
                nn = -1
                for q1 in connected:
                    for q2 in remain:
                        if self.graph.C[pi[q1]][pi[q2]] < md:
                            md = self.graph.C[pi[q1]][pi[q2]]
                            nn = q2
                ud += md - 1
                connected.append(nn)
                remain.remove(nn)
        return ud

class QScheduler:
    def __init__(self, bs, nq):
        pc = []
        for i in range(nq):
            pc.append([])
        i = 0
        for b in bs:
            for q in b:
                pc[q].append(i)
            i += 1
        self.pieces = []
        for p in pc:
            self.pieces.append([False] + p)
        self.placed = []

    def move(self, i):
        self.pieces[i][0] = True
        self.placed.append(i)

    def cdd_pieces(self):
        indp = []
        nq = len(self.pieces)
        for i in range(nq):
            if self.pieces[i][0] or len(self.pieces[i]) == 1:
                continue
            flag = True
            for j in range(nq):
                if j == i or self.pieces[j][0] or len(self.pieces[i]) >= len(self.pieces[j]):
                    continue
                if all(it in self.pieces[j][1:] for it in self.pieces[i][1:]):
                    flag = False
            if flag:
                indp.append(i)

        res = []
        mpriority1 = -1
        mpriority2 = -1
        for pi in indp:
            priority1 = 0
            priority2 = len(self.pieces[pi])
            for pp in self.placed:
                for c in self.pieces[pp][1:]:
                    if c in self.pieces[pi]:
                        priority1 += 1
            if priority1 > mpriority1 or (priority1 == mpriority1 and priority2 > mpriority2):
                res.clear()
                res.append(pi)
                mpriority1 = priority1
                mpriority2 = priority2
            if priority1 == mpriority1 and priority2 == mpriority2:
                res.append(pi)
        return res

class PH_compiler:
    def __init__(self, pauli_blocks, ac):
        self.phycir_path = ''
        self.graph = ac

        self.blocks = []
        for bk in pauli_blocks:
            self.blocks.append(compute_block_cover(bk))
        print(self.blocks[:3])
        self.pauli_blocks = pauli_blocks
        self.lnq = len(pauli_blocks[0][0])
        self.ne = 2

    def set_phycir_path(self, path):
        self.phycir_path = path

    def block_scheduling(self):
        self.board = Board(self.graph, self.blocks)
        self.scheduler = QScheduler(self.blocks, self.lnq)
        self.pauli_layers = [[b] for b in self.pauli_blocks]
    
    def compile(self, opt = 0):
        self.block_scheduling()
        t0 = ctime()
        self.my_pi = simple_qubit_mapping(self.lnq)
        if opt == 1:
            self.my_pi = self.initial_mapping()
        qc = QuantumCircuit(len(self.graph.C))
        for i in range(self.ne):
            qc.x(i)
        # qc.reset(4)
        qc, inner, outer = synthesis_SC.block_opt_SC(self.pauli_layers, graph=self.graph, pauli_map=self.my_pi, qc=qc)
        t1 = ctime()
        if len(self.phycir_path) > 0:
            f = open(self.phycir_path + '_opt' + str(opt) + '_origin.txt', mode='w+')
            f.write(qc.qasm())
            f.close()
        # ncx, nsg, qc = print_qc(qc, opt_level=3)
        if len(self.phycir_path) > 0:
            f = open(self.phycir_path + '_opt' + str(opt) + '.txt', mode='w+')
            f.write(qc.qasm())
            f.close()
        print('ph result: inner: ', inner, ', outer: ', outer, ', depth: ', qc.depth())
        return [inner + outer, outer, qc.depth(), t1 - t0], qc
    
class My_compiler:
    def __init__(self, pauli_blocks, ac):
        self.phycir_path = ''
        self.graph = ac

        self.blocks = []
        for bk in pauli_blocks:
            self.blocks.append(compute_block_cover(bk))
        # print(self.blocks[:3])
        self.pauli_blocks = pauli_blocks
        self.lnq = len(pauli_blocks[0][0])
        self.ne = 2

    def set_phycir_path(self, path):
        self.phycir_path = path
    

    def initial_mapping(self):
        self.my_pi = {}
        while True:
            cdd_p = self.scheduler.cdd_pieces()
            if len(cdd_p) == 0:
                break
            c = -1
            pos = -1
            mscore = 1000000000
            for posi in self.board.edge:
                for ci in cdd_p:
                    score = self.board.score(self.scheduler.pieces[ci][1:], posi)
                    if score < mscore:
                        mscore = score
                        c = ci
                        pos = posi
            self.scheduler.move(c)
            self.board.move(self.scheduler.pieces[c][1:], pos)
            self.my_pi[c] = pos
            print(c, ' ', end="")
        print('\n')
        return self.my_pi
    
    def block_scheduling(self):
        self.board = Board(self.graph, self.blocks)
        self.scheduler = QScheduler(self.blocks, self.lnq)
        self.pauli_layers = [[b] for b in self.pauli_blocks]
    
    def compile(self, opt=2):
        self.block_scheduling()
        t0 = ctime()
        self.initial_mapping()
        # print(self.my_pi)
        if opt == 1:
            self.my_pi = simple_qubit_mapping(self.lnq)
        qc = QuantumCircuit(len(self.graph.C))
        # qc.reset(4)
        for i in range(self.ne):
            qc.x(i)
        qc, inner, outer = synthesis_SC.block_opt_SC(self.pauli_layers, graph=self.graph, pauli_map=self.my_pi, qc=qc, synthesis_opt=True)
        t1 = ctime()
        if len(self.phycir_path) > 0:
            f = open(self.phycir_path + '_opt' + str(opt) + '_origin.txt', mode='w+')
            f.write(qc.qasm())
            f.close()
        # ncx, nsg, qc = print_qc(qc, opt_level=3)
        if len(self.phycir_path) > 0:
            f = open(self.phycir_path + '_opt' + str(opt) + '.txt', mode='w+')
            f.write(qc.qasm())
            f.close()
        print('my result: inner: ', inner, ', outer: ', outer, ', depth: ', qc.depth())
        return [inner + outer, outer, qc.depth(), t1 - t0], qc