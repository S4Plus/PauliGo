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
    return blocks

def compile():
    with open('./benchmark/H2.pickle', 'rb') as f:
        op_list = pickle.load(f)
    blocks = program_prep(op_list)
    compiler = Compiler(blocks, 'grid0505')
    qc = compiler.start('ucc')
    print(qc.count_ops())
    print(qc.qasm())
    # nSWAP = qc.count_ops()['swap']
    # depth = qc.depth()
    # print(nSWAP, ' ', depth)

compile()

