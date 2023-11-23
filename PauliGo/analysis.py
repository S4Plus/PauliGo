from benchmark.offline import *
from parallel_bl import *
from qiskit import QuantumCircuit, transpile
import synthesis_SC

from tools import *
from arch import *
from synthesis_sd import *
import pandas as pd
import time, sys, os
from functions import *
from benchmark.hami import gene_random_oplist

moles = ['LiH', 'BeH2', 'CH4', 'MgH', 'LiCl', 'CO2']
orbital = [8, 12, 16, 20, 24, 28]
test_method = 'random_code'
ac = 'grid'

def get_oplist(t, n):
    if t == 'uccsd':
        return load_oplist(n, benchmark='uccsd')
    if t == 'random':
        return gene_random_oplist(n)
if test_method == 'intstr':
    i = 0
    print('UCCSD:', orbital[i])
    parr = load_oplist(moles[i], benchmark='uccsd')
    lnq = len(parr[0][0])
    length = lnq // 2
    a2 = depth_oriented_scheduling(parr, length=length, maxiter=30)
    intstr = [0] * lnq * lnq
    for i1 in a2:
        for i2 in i1:
            lcover = compute_block_cover(i2)
            if len(lcover) == 1:
                continue
            for j in range(len(lcover) - 1):
                for k in range(j + 1, len(lcover)):
                    intstr[lcover[j] * lnq + lcover[k]] += 1
                    intstr[lcover[k] * lnq + lcover[j]] += 1
    for i1 in range(lnq):
        for i2 in range(lnq):
            print(intstr[i1 * lnq + i2], ' ', end = "")
        print('\n')

if test_method == 'swap':
    d = {}
    d['name'] = []
    d['nps'] = []
    d['inner'] = []
    d['outer'] = []
    d['depth'] = []
    for i in range(6):
        print('UCCSD:', orbital[i])
        d['name'].append('UCCSD:'+str(orbital[i]))
        parr = load_oplist(moles[i], benchmark='uccsd')
        lnq = len(parr[0][0])
        length = lnq // 2
        a2 = depth_oriented_scheduling(parr, length=length, maxiter=30)
        nps = 0
        for i1 in a2:
            for i2 in i1:
                nps += len(i2)
        d['nps'].append(nps)
        qc, inner, outer = synthesis_SC.block_opt_SC(a2, arch=ac)
        d['depth'].append(circtui_depth(qc))
        d['inner'].append(inner)
        d['outer'].append(outer)
    D = pd.DataFrame(d)
    D.to_excel('D:/DAC/data/metric/uccsd.xlsx')
    print(D)

if test_method == 'random':
    d = {}
    d['name'] = []
    d['nps'] = []
    d['inner'] = []
    d['outer'] = []
    d['depth'] = []
    for i in range(6):
        print('random:', orbital[i])
        d['name'].append('random:' + str(orbital[i]))
        parr = gene_random_oplist(orbital[i])
        lnq = len(parr[0][0])
        length = lnq // 2
        a2 = depth_oriented_scheduling(parr, length=length, maxiter=30)
        nps = 0
        for i1 in a2:
            for i2 in i1:
                nps += len(i2)
        print('nps: ', nps)
        d['nps'].append(nps)
        qc, inner, outer = synthesis_SC.block_opt_SC(a2, arch=ac)
        d['depth'].append(circtui_depth(qc))
        d['inner'].append(inner)
        d['outer'].append(outer)
    D = pd.DataFrame(d)
    D.to_excel('D:/DAC/data/metric/random.xlsx')
    print(D)

if test_method == 'code':
    for i in range(6):
        print('UCCSD:', orbital[i])
        parr = load_oplist(moles[i], benchmark='uccsd')
        lnq = len(parr[0][0])
        length = lnq // 2
        a2 = depth_oriented_scheduling(parr, length=length, maxiter=30)
        f = open('D:/DAC/data/uccsd_'+str(orbital[i])+'_str.txt', mode='w')
        for i1 in a2:
            for i2 in i1:
                lcover = compute_block_cover(i2)
                for j in range(orbital[i]):
                    if j in lcover:
                        f.write('*')
                    else:
                        f.write(" ")
                f.write('\n')
                for i3 in i2:
                    lcover = ps2nodes(i3.ps)
                    for j in range(orbital[i]):
                        if j in lcover:
                            f.write('+')
                        else:
                            f.write("-")
                    f.write('\n')
        f.close()

if test_method == 'random_code':
    for i in range(6):
        print('random:', orbital[i])
        parr = gene_random_oplist(orbital[i])
        lnq = len(parr[0][0])
        length = lnq // 2
        a2 = depth_oriented_scheduling(parr, length=length, maxiter=30)
        f = open('D:/DAC/data/random_'+str(orbital[i])+'_str.txt', mode='w')
        for i1 in a2:
            for i2 in i1:
                for i3 in i2:
                    lcover = ps2nodes(i3.ps)
                    for j in range(orbital[i]):
                        if j in lcover:
                            f.write('+')
                        else:
                            f.write("-")
                    f.write('\n')
        f.close()