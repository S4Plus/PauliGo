from benchmark.offline import *
from compiler import Compiler
from functions import compute_block_cover
from qiskit import transpile
from time import time


uccsd_index = [['u{}'.format(i), 'uccsd_{}.pickle'.format(i)] for i in[8, 12, 16, 20, 24, 28]] # 8, 12, 16, 20, 24, 28
random_index = [['r{}'.format(i), 'random_{}.pickle'.format(i)] for i in range(8,22,2)]
molecules = [[str(i), '{}.pickle'.format(i)] for i in ['H2', 'LiH', 'H2O', 'NH3', 'CH4']]
bms = molecules + random_index + uccsd_index # random_index + molecules + uccsd_index
cps = ['ucc'] # 'tk', 'ph', 'go'
acs = ['manhattan', 'sycamore']

def to_txt(bm, path):
    op_list = load_benchmark(bm)
    print(type(op_list[0][0]), op_list[0][0])
    f = open(path, 'w+')
    f.write(str(op_list))
    f.close()

for bm in bms:
    path = './benchmark/text/{}'.format(bm[1].replace('.pickle', '.txt'))
    to_txt(bm[1], path)
exit()

def program_prep(origin_blocks):
    bn = pn = 0
    blocks = []
    for bk in origin_blocks:
        if (len(compute_block_cover(bk))) > 0:
            blocks.append(bk)
            pn += len(bk)
    bn = len(origin_blocks)
    return blocks # , bn, pn

def compile(bm, ac, cp, opt=0):
    op_list = load_benchmark(bm)
    blocks = program_prep(op_list)
    compiler = Compiler(blocks, ac)
    qc = compiler.start(cp, opt)
    count_res = {}
    count_res['swap'] = qc.count_ops()['swap']
    t0 = time()
    qc = transpile(qc, basis_gates=['cx', 'u3'], optimization_level=3)
    print("Qiskit L3:", time()-t0)
    count_res['cx'] = qc.count_ops()['cx']
    count_res['depth'] = qc.depth()
    print(count_res, '\n')
    return count_res

# f = open('./data/tk_res.pickle', '+wb')
res = {}
for bm in bms:
    res[bm[0]] = {}
    for ac in acs:
        res[bm[0]][ac] = {}
        for cp in cps:
            print(bm[0], ' ', ac, ' ', cp)
            res[bm[0]][ac][cp] = compile(bm[1], ac, cp)
f1 = open('./data/tk_res1.txt', 'w+')
f1.write(str(res))
f1.close()
