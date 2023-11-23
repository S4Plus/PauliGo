from functions import *
from arch import *
import synthesis_SC
from tools import print_qc
from parallel_bl import depth_oriented_scheduling

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

class Compiler:
    def __init__(self, pauli_blocks, ac):
        G, C = load_graph(ac, dist_comp=True)
        self.graph = pGraph(G, C)

        blocks = []
        for bk in pauli_blocks:
            blocks.append(compute_block_cover(bk))
        print(blocks[:5])
        self.board = Board(self.graph, blocks)
        self.scheduler = QScheduler(blocks, len(pauli_blocks[0][0]))

        self.lnq = len(pauli_blocks[0][0])
        self.pauli_layers = depth_oriented_scheduling(pauli_blocks, length=self.lnq // 2, maxiter=30)

    def set_phycir_path(self, path):
        self.phycir_path = path
    
    def ph_compile(self, opt=0):
        self.my_pi = dummy_mapping(self.lnq)
        if opt==1:
            self.initial_mapping()
        qc, inner, outer = synthesis_SC.block_opt_SC(self.pauli_layers, graph=self.graph, pauli_map=self.my_pi)
        f = open(self.phycir_path + '_origin.txt', mode='w+')
        f.write(qc.qasm())
        f.close()
        # ncx, nsg, qc = print_qc(qc, opt_level=3)
        # f = open(self.phycir_path + '.txt', mode='w+')
        # f.write(qc.qasm())
        # f.close()
        print('ph result: inner: ', inner, ': outer: ', outer, '\n')
        return [0, 0, inner + outer, outer, qc.depth()]

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
    
    def go_compile(self, opt=2):
        self.initial_mapping()
        if opt == 1:
            self.my_pi = dummy_mapping(self.lnq)
        qc, inner, outer = synthesis_SC.block_opt_SC(self.pauli_layers, graph=self.graph, pauli_map=self.my_pi, synthesis_opt=True)
        f = open(self.phycir_path + '_opt' + str(opt) + '_origin.txt', mode='w+')
        f.write(qc.qasm())
        f.close()
        # ncx, nsg, qc = print_qc(qc, opt_level=3)
        # f = open(self.phycir_path + '_opt' + str(opt) + '.txt', mode='w+')
        # f.write(qc.qasm())
        # f.close()
        print('my result: inner: ', inner, ', outer: ', outer, ', depth: ', qc.depth())
        return [0, 0, inner + outer, outer, qc.depth()]