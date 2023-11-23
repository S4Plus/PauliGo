from benchmark.offline import *
from compiler import Compiler
from functions import compute_block_cover
from molecules import molecule_oplist
from benchmark.hami import gene_random_oplist

phycir_path = './data/phycir/'
uccsd_moles = ['LiH', 'BeH2', 'CH4', 'MgH', 'LiCl', 'CO2']  # 

def program_prep(origin_blocks):
    bn = pn = 0
    blocks = []
    for bk in origin_blocks:
        if (len(compute_block_cover(bk))) > 0:
            blocks.append(bk)
            pn += len(bk)
    bn = len(origin_blocks)
    return blocks, bn, pn

def uccsd(ac, compiler_name, opt=2):
    f = open('./data/v3/uccsd_' + ac + '_' + compiler_name + '_opt' + str(opt) + '.txt', mode='w+')
    for bm in uccsd_moles:
        origin_blocks = load_benchmark(bm + '_UCCSD.pickle')
        name = 'uccsd_' + str(len(origin_blocks[0][0]))
        blocks, bn, pn = program_prep(origin_blocks)
        compiler = Compiler(blocks, ac)
        item = [name, bn, pn]
        compiler.set_phycir_path(phycir_path + compiler_name + '/' + name + '_' + ac)  #  + compiler_name + '/'
        if compiler_name == 'ph':
            item += compiler.ph_compile(opt=opt)
        if compiler_name == 'go':
            item += compiler.go_compile(opt=opt)
        print(item)
        for i in item:
            f.write(str(i) + ' ')
        f.write('\n')
    f.close()

def molecule(ac, compiler_name, opt=2):
    f = open('./data/v3/molecule_' + ac + '_' + compiler_name + '_opt' + str(opt) + '.txt', mode='w+')
    from molecules import molecules
    for mole in molecules:
        origin_blocks = molecule_oplist(mole[1])
        name = mole[0] + '_' + str(len(origin_blocks[0][0]))
        print('\n', name)
        blocks, bn, pn = program_prep(origin_blocks)
        compiler = Compiler(blocks, ac)
        item = [name, bn, pn]
        compiler.set_phycir_path(phycir_path + compiler_name + '/' + name + '_' + ac)  #  + compiler_name + '/'
        if compiler_name == 'ph':
            item += compiler.ph_compile(opt=opt)
        if compiler_name == 'go':
            item += compiler.go_compile(opt=opt)
        print(item, '\n')
        for i in item:
            f.write(str(i) + ' ')
        f.write('\n')
    f.close()

def random_hami(ac, compiler_name, opt=2):
    f = open('./data/v3/random_' + ac + '_' + compiler_name + '_opt' + str(opt) + '.txt', mode='w+')
    seeds = [83, 193, 239, 367, 439, 571, 661, 743, 881, 967, 1049, 1153, 1201, 1367, 1489, 1543, 1621, 1787, 1889, 1949]
    for nq in range(8, 22, 2):
        origin_blocks = gene_random_oplist(nq, 5*nq*nq, seed=seeds[nq//2-4])
        name = 'random_' + str(len(origin_blocks[0][0]))
        blocks, bn, pn = program_prep(origin_blocks)
        compiler = Compiler(blocks, ac)
        item = [name, bn, pn]
        compiler.set_phycir_path(phycir_path + compiler_name + '/' + name + '_' + ac)  #  + compiler_name + '/'
        if compiler_name == 'ph':
            item += compiler.ph_compile(opt=opt)
        if compiler_name == 'go':
            item += compiler.go_compile(opt=opt)
        print(item)
        for i in item:
            f.write(str(i) + ' ')
        f.write('\n')
    f.close()

acs = ['sycamore', 'manhattan']
compiler_names = ['ph', 'go']
benchmarks = ['random']
opt_levels = [1,2]

for ac in acs:
    for compiler_name in compiler_names:
        for benchmark in benchmarks:
            for opt_level in opt_levels:
                if benchmark == 'molecule':
                    molecule(ac, compiler_name, opt=opt_level)
                if benchmark == 'uccsd':
                    uccsd(ac, compiler_name, opt=opt_level)
                if benchmark == 'random':
                    random_hami(ac, compiler_name, opt=opt_level)
