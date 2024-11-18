from benchmark.offline import *
from compiler import Compiler
from functions import compute_block_cover
from qiskit import transpile
from time import time
from copy import deepcopy
from benchmark.mypauli import *


def program_prep(origin_blocks):
    bn = pn = 0
    blocks = []
    for bk in origin_blocks:
        blocks.append([pauliString(ps[0], ps[1]) for ps in bk])
    bn = len(origin_blocks)
    return blocks # , bn, pn

def compile():
    with open('./benchmark/H2.pickle', 'rb') as f:
        op_list = pickle.load(f)
    # print(op_list)
    # exit()
    blocks = program_prep(op_list)
    print(blocks)
    compiler = Compiler(blocks, 'all4')
    qc = compiler.start('ph')
    # count_res = {}
    # count_res['swap'] = qc.count_ops()['swap']
    # t0 = time()
    qc = transpile(qc, basis_gates=['cx', 'u3'], optimization_level=3)
    # print("Qiskit L3:", time()-t0)
    count_cx = qc.count_ops()['cx']
    # count_res['depth'] = qc.depth()
    # print(count_res, '\n')
    qc.draw('mpl', filename = './H2.png')
    print(count_cx)

compile()

