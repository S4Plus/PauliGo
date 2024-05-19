cps = ['tk', 'ucc', 'ph', 'go'] # , 'ph', 'go'
acs = ['manhattan', 'sycamore']
uccsd_index = ['u{}'.format(i) for i in[8, 12, 16, 20, 24, 28]] # 8, 12, 16, 20, 24, 28
random_index = ['r{}'.format(i) for i in range(8,22,2)]
molecules = ['H2', 'LiH', 'H2O', 'NH3', 'CH4']
bms = uccsd_index + random_index + molecules
metrics = ['swap', 'cx', 'depth']

result = open('../data/all_res.txt', 'w+')
def print_dict_tree(dictionary, indent=0):
    for key, value in dictionary.items():
        if isinstance(value, dict):
            result.write('  ' * indent + str(key) + ':' + '\n')
            print_dict_tree(value, indent + 1)
        else:
            result.write('  ' * indent + str(key) + ': ' + str(value) + '\n')

import ast

def parser(f):
    res = {}
    for ac in acs:
        res[ac] = {}
        for cp in cps:
            res[ac][cp] = {}
            for bm in bms:
                res[ac][cp][bm] = {}

    indexs = []
    for line in f:
        if len(line) == 0:
            continue
        # print(line)
        # input()
        if line[0] not in ['Q', 'T', '{']:
            indexs = line.split()
            # print(indexs)
        if line[0] == '{':
            d = ast.literal_eval(line)
            res[indexs[1]][indexs[2]][indexs[0]] = d
    return res

f = open('../data/tk.txt', 'r')
ori_data = parser(f)
f = open('../data/tk_res1.pickle', 'r')
for line in f:
    m1 = ast.literal_eval(line)
    for k1,v1 in m1.items():
        for k2,v2 in v1.items():
            for k3,v3 in v2.items():
                for k4, v4 in v3.items():
                    ori_data[k2][k3][k1][k4] = v4
    # ori_data.update(ast.literal_eval(line))
# print_dict_tree(ori_data)
# input()
f = open('../data/tk_res1.txt', 'r')
for line in f:
    m1 = ast.literal_eval(line)
    for k1,v1 in m1.items():
        for k2,v2 in v1.items():
            for k3,v3 in v2.items():
                for k4, v4 in v3.items():
                    ori_data[k2][k3][k1][k4] = v4

result.write(str(ori_data))
result.close()
input()

f = open('./table.txt', 'w+')
for bm in bms:
    f.write(bm)
    for ac in acs:
        for metric in ['swap', 'depth']:
            x = []
            for cp in cps:
                x.append(ori_data[ac][cp][bm][metric])
            for xi in x:
                if xi == min(x):
                    f.write(' & \\textbf{' + str(xi) + '}')
                else :
                    f.write(' & ' + str(xi))
    f.write('\\\\\n')
f.close()
# input()

import matplotlib.pyplot as plt
import numpy as np

def draw_grouped_bar_chart(xlabel, data, labels, save_path):
    labels = ['PauliSimp', 'UCCSynthesis', 'ph', 'ours']
    n_bars = len(data)
    colors = ['deepskyblue', 'darkviolet', 'black']
    # for i in [xlabel, data, labels]:
    #     print(i)
    x = np.arange(len(xlabel))
    width = 0.2
    fig, ax = plt.subplots(figsize=(6, 4))
    # y = n_bars*width/2 + 0*width
    # print(y, '\n', x - y)
    styles = [{'edgecolor' : 'black', 'hatch' : '//', 'fill' : False}, 
              {'edgecolor' : 'black', 'hatch' : '\\\\', 'fill' : False},
              {'edgecolor' : 'black', 'fill' : False}, 
              {'color' : 'black'}]
    for i in range(n_bars):
        ax.bar(x - n_bars*width/2 + i*width, data[i], width, label=labels[i], **styles[i])
    ax.set_xlabel('# qubits', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(xlabel, fontsize=12)
    ax.legend(loc='upper left', fontsize=12)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

bms = [uccsd_index, random_index, molecules]
for ac in acs:
    for metric in ['swap', 'depth']:
        for bm in bms:
            data = []
            for cp in cps:
                d = []
                for bi in bm:
                    # print('{} {} {} {}'.format(ac, cp, bi, metric))
                    d.append(ori_data[ac][cp][bi][metric])
                data.append(d)
            bm_name = bm[0][0]
            if bm[0][0] == 'H':
                bm_name = 'm'
            draw_grouped_bar_chart(bm, data, cps, './{}_{}_{}.png'.format(ac, metric, bm_name))
