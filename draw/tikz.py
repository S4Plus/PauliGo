
class gate():
    def __init__(self, name, q):
        self.q = q
        self.name = name

def tikz(gate_list, nq, fname):
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
    f = open(fname, 'w+')
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
    f.close()

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
cir = circuit()
cir.many_cx([[5,4], [4,3]])
tikz(cir.gate_list, 6, './cir2.txt')
cir = circuit()
cir.many_cx([[5,4], [4,3], [1, 0], [2, 0], [3, 0]])
tikz(cir.gate_list, 6, './cir3.txt')
cir = circuit()
cir.many_cx([[5,4], [4,3], [3, 0], [2, 0], [1, 0]])
tikz(cir.gate_list, 6, './cir4.txt')