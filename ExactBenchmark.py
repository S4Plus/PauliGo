import random
from arch import *
from synthesis_SC import tree
from functions import generate1, print_tree1
from compiler import Compiler
from benchmark.mypauli import pauliString
from qiskit import QuantumCircuit, QuantumRegister
class ExactBenchmark:
    def __init__(self, G):
        self.G = G
    def get_bm(self, nb, ms, Ms):
        bm = []
        size = range(ms, Ms)
        for i in range(nb):
            nq = random.choice(size)
            seed = random.choice(range(len(self.G)))
            b = []
            b.append(seed)
            for j in range(nq - 1):
                leafs = []
                for q in b:
                    for ad in self.G[q]:
                        if ad in b:
                            continue
                        if ad in leafs:
                            continue
                        leafs.append(ad)
                #print(leafs)
                b.append(random.choice(leafs))
            bm.append(b)
        qset = set()
        for b in bm:
            qset = qset | set(b)
        pqs = list(qset)  # all physical qubits
        nq = len(pqs)
        random.shuffle(pqs)
        PI = {}
        for i in range(nq):
            PI[pqs[i]] = i  # physical to logical
        pi = {}
        for a, b in PI.items():
            pi[b] = a  # logical to physical
        lbs = []
        for b in bm:
            lb = []
            for q in b:
                lb.append(PI[q])
            lbs.append(lb)
        return lbs, pi

G, C = load_graph('grid0303', dist_comp=True)
graph = []
for i in range(len(G)):
    graph.append([])
for i in range(len(G)):
    for j  in range(len(G)):
        if G[i][j] > 0:
            graph[i].append(j)
eb = ExactBenchmark(graph)
import time
random.seed(time.time())
blocks, pi = eb.get_bm(5, 5, 9)
# blocks = [[6, 7, 1, 4], [0, 8, 6, 7, 5, 1], [3, 4, 5, 7], [3, 4, 2, 8], [1, 7, 4]]
blocks = [[0, 2, 4, 5, 6, 1, 3, 8], [6, 8, 1, 4, 3, 7, 0], [4, 0, 5, 1, 7, 6], [6, 8, 4, 7, 0, 5], [3, 5, 1, 4, 0, 6, 2, 7]]
print(blocks)
# print(pi)
pauli_blocks = []
for b in blocks:
    ps = ''
    for i in range(9):
        if i in b:
            ps += random.choice(['X','Y','Z'])
        else:
            ps += 'I'
    print(ps)
    pauli_blocks.append([pauliString(ps, 1)])
compiler = Compiler(pauli_blocks, 'grid0303')
# qc = compiler.go_compile()
my_pi = compiler.initial_mapping()
print(my_pi)
graph = pGraph(G, C)
def l2p(b, pi):
    return [pi[i] for i in b]

def go_synthesis(graph, tree, qc:QuantumCircuit, PI):
    if len(tree.children) > 0:
        for i in range(len(tree.children)):
            go_synthesis(graph, tree.children[i], qc, PI)
        tree.children.sort(key=lambda x: x.t)
        tree.t = max([tree.children[i].t + i + 1 for i in range(len(tree.children))])
        for child in tree.children:
            qc.cx(PI[child.number], PI[tree.number])
    else:
        tree.t = 0

def go_synthesis1(graph, qc, nodes, pi):
    res = 0
    PI = {}
    for k,v in pi.items():
        PI[v] = k
    tree = generate1(graph, nodes)
    go_synthesis(graph, tree, qc, PI)
    return res


i = 0
for b in blocks:
    Q = QuantumRegister(9, name='')
    qc = QuantumCircuit(Q)
    print(b)
    go_synthesis1(graph, qc, l2p(b, my_pi), my_pi)
    # qc.draw('latex', filename='./draw/cir_{}.png'.format(i), idle_wires=False)
    qc.draw('latex_source', filename='./draw/circuits/cir_{}.tex'.format(i), idle_wires=False)
    i += 1
    print(qc.qasm())
    input()

# print(qc.qasm())
# print(compiler.my_pi)
# compiler.ph_compile()