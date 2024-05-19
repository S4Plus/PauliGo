import random
from arch import *
from synthesis_SC import tree
from functions import generate, print_tree1
from compiler import Compiler
from benchmark.mypauli import pauliString
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
blocks, pi = eb.get_bm(7, 3, 7)
blocks = [[7, 2, 3, 4, 0], [3, 2, 7, 4], [1, 4, 7], [5, 1, 2, 0, 6], [4, 0, 1, 7], [7, 2, 3], [3, 2, 1, 4, 6, 7]]
# blocks = [[7, 3, 5, 2, 0], [2, 7, 3, 0, 8, 4], [7, 4, 3, 6, 5, 1], [6, 2, 3, 7, 5, 4], [5, 3, 7, 0], [5, 8, 0, 1], [6, 2, 1, 0, 8, 7]]
print(blocks)
# print(pi)
pauli_blocks = []
for b in blocks:
    ps = ''
    for i in range(9):
        if i in b:
            ps += 'Z'
        else:
            ps += 'I'
    print(ps)
    pauli_blocks.append([pauliString(ps, 1)])
compiler = Compiler(pauli_blocks, 'grid0303')
my_pi = compiler.initial_mapping()
print(my_pi)
graph = pGraph(G, C)
def l2p(b, pi):
    return [pi[i] for i in b]
for b in blocks:
    print(l2p(b, my_pi))
    generate(graph, l2p(b, my_pi))
# compiler.go_compile()
# print(compiler.my_pi)
# compiler.ph_compile()