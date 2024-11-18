
class gate():
    def __init__(self, name, q):
        self.q = q
        self.name = name

def tikz(gate_list, nq, f):
    m = [[] for i in range(nq)]
    for g in gate_list:
        if g.name == 'cx':
            q1, q2 = g.q[0], g.q[1]
            d1 = max(len(m[q1]),len(m[q2]))
            for q in g.q:
                while len(m[q]) < d1:
                    m[q].append(0)
            m[q1].append('\\ctrl{{{}}}'.format(q2-q1))
            m[q2].append('\\targ{}')
        elif g.name == 'swap':
            q1, q2 = g.q[0], g.q[1]
            d1 = max(len(m[q1]),len(m[q2]))
            for q in g.q:
                while len(m[q]) < d1:
                    m[q].append(0)
            m[q1].append('\\swap{{{}}}'.format(q2-q1))
            m[q2].append('\\targX{}')
        else:
            m[g.q[0]].append(g.name)
    d2 = max([len(mi) for mi in m])
    for i in range(nq):
        while len(m[i]) < d2:
            m[i].append(0)
    # f = open(fname, 'w+')
    f.write('\n\n')
    k = 0
    for i in m:
        f.write('\\lstick{{${}$}}'.format(k))
        for j in i:
            if j == 0:
                f.write(' &')
            else:
                f.write(' & ' + j)
        k += 1
        f.write(' & \\\\\n')
    # f.close()

class circuit():
    def __init__(self):
        self.gate_list=[]
    def cx(self, q1, q2):
        self.gate_list.append(gate('cx',[q1,q2]))
    def swap(self, q1, q2):
        self.gate_list.append(gate('swap',[q1,q2]))
    def rz(self, q):
        self.gate_list.append(gate('\\gate{{R_z(\\theta)}}',[q]))
    def many_cx(self, gl):
        for g in gl:
            self.cx(g[0],g[1])
    def exchange(self, q1, q2):
        for g in self.gate_list:
            g.q = [q2 if i == q1 else i for i in g.q]
cir = circuit()
f = open('./example.txt', 'w+')
cir.many_cx([[5, 6] , [1, 3] , [2, 0] , [3, 0] , [6, 4] , [0, 4]])
# cir.exchange(0, 2)
# cir.exchange(1, 0)
# cir.exchange(2, 1)
tikz(cir.gate_list, 9, f)
cir = circuit()
cir.many_cx([[5, 7] , [2, 0] , [1, 4] , [7, 4] , [0, 4]])
tikz(cir.gate_list, 9, f)
cir = circuit()
cir.many_cx([[6, 8] , [3, 1] , [8, 1] , [0, 4] , [7, 4] , [1, 4]])
tikz(cir.gate_list, 9, f)
cir = circuit()
cir.many_cx([[0, 4] , [1, 4] , [6, 4]])
tikz(cir.gate_list, 9, f)
cir = circuit()
cir.many_cx([[2, 7] , [4, 7] , [5, 7]])
tikz(cir.gate_list, 9, f)
f.close()