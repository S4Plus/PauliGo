import os
import random

from functions import *
from benchmark.offline import *
from synthesis_sd import *
from parallel_bl import *
import synthesis_SC
from qiskit import transpile
from tools import print_qc

<<<<<<< HEAD
output_circuit_path = '../data/physical circuit/manhattan/'
# moles = ['LiH', 'BeH2', 'CH4', 'MgH', 'LiCl', 'CO2']
=======
output_circuit_path = './data/phycir/'
moles = ['LiH', 'BeH2', 'CH4', 'MgH', 'LiCl', 'CO2']  # 
>>>>>>> f0dd2b047923bfaf2931d32c0ad192f0dd14a3b2
# orbital = [8, 12, 16, 20, 24, 28]
# moles = ['H2S']
ac = 'manhattan'
G, C = load_graph(ac, dist_comp=True)
graph = pGraph(G, C)
pnq = len(G)


def start(parr, compiler='go', opt=2, show_mapping = False):
    res = []
    print('block number: ', len(parr))
    lnq = len(parr[0][0])
    print('qubit: ', lnq)
    length = lnq // 2
    a2 = depth_oriented_scheduling(parr, length=length, maxiter=30)
    res.append(lnq)
    res.append(len(parr))
    blocks = []
    block_wei = []
    layer = 0
    for bl in a2:
        for bk in bl:
            blocks.append(compute_block_cover(bk))
            block_wei.append(layer)
        layer += 1
    print(blocks[:5])
    board = Board(graph, blocks, block_wei)
    scheduler = QScheduler(blocks, lnq, block_wei)
    my_pi = {}
    while True:
        cdd_p = scheduler.cdd_pieces(board)
        if len(cdd_p) == 0:
            break
        c = -1
        pos = -1
        mscore = 1000000000
        for pi in board.edge:
            for ci in cdd_p:
                score = board.score(scheduler.pieces[ci][1:], pi)
                if score < mscore:
                    mscore = score
                    c = ci
                    pos = pi
        scheduler.move(c)
        board.move(scheduler.pieces[c][1:], pos)
        my_pi[c] = pos
        print(c, ' ', end="")
    print('\n')
    paulis = []
    for bl in a2:
        for bk in bl:
            for p in bk:
                paulis.append(ps2nodes(p.ps))
    print('pauli number: ', len(paulis))
    res.append(len(paulis))

    unused = [i for i in range(pnq)]
    for k, v in my_pi.items():
        unused.remove(v)
    for i in range(lnq):
        use = -1
        if i not in my_pi.keys():
            use = unused[0]
            unused.remove(use)
        if use != -1:
            my_pi[i] = use

    if show_mapping:
        my_PI = {}
        for i in range(pnq):
            my_PI[i] = -1
        for k, v in my_pi.items():
            my_PI[v] = k
        print_mapping(my_PI, ac)

    if compiler == 'go':
        unconnect = board.unconnected_degree(paulis, my_pi)
        print('my unconnected degree: ', unconnect)
        res.append(unconnect)
        if opt == 1:
            my_pi = dummy_mapping(lnq)
        qc, inner, outer = synthesis_SC.block_opt_SC(a2, arch=ac, pauli_map=my_pi, synthesis_opt=True)
        f = open(output_circuit_path + 'go/' + name + '_' + ac + '_opt' + str(opt) + '_origin.txt', mode='w+')
        f.write(qc.qasm())
        f.close()
        ncx, nsg, qc = print_qc(qc, opt_level=3)
        f = open(output_circuit_path + 'go/' + name + '_' + ac + '_opt' + str(opt) + '.txt', mode='w+')
        f.write(qc.qasm())
        f.close()
        print('my result: inner: ', inner, ', outer: ', outer, ', depth: ', qc.depth())
        res.append(ncx)
        res.append(inner + outer)
        res.append(outer)
        res.append(qc.depth())
    else:
        dummy_pi = dummy_mapping(lnq)
        unconnect = board.unconnected_degree(paulis, dummy_pi)
        print('dummy unconnected degree: ', unconnect)
        res.append(unconnect)

        qc, inner, outer = synthesis_SC.block_opt_SC(a2, arch=ac, pauli_map=dummy_pi)
        f = open(output_circuit_path + 'ph/' + name + '_' + ac + '_origin.txt', mode='w+')
        f.write(qc.qasm())
        f.close()
        ncx, nsg, qc = print_qc(qc, opt_level=3)
        f = open(output_circuit_path + 'ph/' + name + '_' + ac + '.txt', mode='w+')
        f.write(qc.qasm())
        f.close()
        print('ph result: inner: ', inner, ': outer: ', outer, '\n')
        res.append(ncx)
        res.append(inner + outer)
        res.append(outer)
        res.append(qc.depth())
    return res


bm_path = './benchmark/data'
<<<<<<< HEAD
=======
comp_name = 'ph'
opt_pass = 2
>>>>>>> f0dd2b047923bfaf2931d32c0ad192f0dd14a3b2
d = {}
columns = ['name', 'qubits', 'blocks', 'paulis',
           'uc', 'gate', 'swaps', 'outer', 'depth']
for col in columns:
    d[col] =[]
bms = os.listdir(bm_path)
# for bm in bms:
    # if 'UCCSD' not in bm:
    #     continue
    # print(bm)
    # name = bm.replace('.pickle', '')
    # parr1 = load_benchmark(bm)
for bm in moles:
    parr1 = load_benchmark(bm + '_UCCSD.pickle')
    name = 'uccsd_' + str(len(parr1[0][0]))
    parr = []
    for bk in parr1:
        if (len(compute_block_cover(bk))) > 0:
            parr.append(bk)
    item = [name]
    item += start(parr, compiler=comp_name, opt=opt_pass, show_mapping=False)
    for i in range(len(item)):
        d[columns[i]].append(item[i])
    print('\n')
# D = pd.DataFrame(d)
f = open('./data/v2/uccsd_' + ac + '_' + comp_name + '_opt' + str(opt_pass) + '.txt', mode='w+')
for i in range(len(d['name'])):
    for j in columns:
        f.write(str(d[j][i]) + ' ')
    f.write('\n')
f.close()
