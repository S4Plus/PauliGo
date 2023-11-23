import random

class ExactBenchmark:
    def __init__(self, G):
        self.G = G
    def get_bm(self, nb, ms):
        bm = []
        size = range(ms, len(self.G) + 1)
        for i in range(nb):
            nq = random.choice(size)
            seed = random.choice(range(len(self.G)))
            b = []
            b.append(seed)
            for j in range(nq - 1):
                leafs = []
                for q in b:
                    for adj in self.G[q]:
                        if adj in b:
                            continue
                        if adj in leafs:
                            continue
                        leafs.append(adj)
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