import ast
import matplotlib.pyplot as plt
import numpy as np

cps = ['tk', 'ucc', 'ph', 'go'] # , 'ph', 'go'
acs = ['manhattan', 'sycamore']
uccsd_index = ['u{}'.format(i) for i in[8, 12, 16, 20, 24, 28]] # 8, 12, 16, 20, 24, 28
random_index = ['r{}'.format(i) for i in range(8,22,2)]
molecules = ['H2', 'LiH', 'H2O', 'NH3', 'CH4']
bms = uccsd_index + random_index + molecules
metrics = ['cx', 'depth']

f = open('../all_res.txt', 'r')
for line in f:
    data = ast.literal_eval(line)

f = open('./table.txt', 'w+')
for ac in acs:
    # f.write(ac)
    for bm in bms:
        # f.write(' & ' + bm)
        # for metric in metrics:
        #     x = []
        #     for cp in cps:
        #         x.append(data[ac][cp][bm][metric])
        #     for xi in x:
        #         if xi == min(x):
        #             f.write(' & \\textbf{' + str(xi) + '}')
        #         else :
        #             f.write(' & ' + str(xi))
        # f.write('\\\\\n')
        f.write(bm)
        for cp in cps:
            for metric in metrics:
                f.write(' & ' + str(data[ac][cp][bm][metric]))
        f.write('\\\\\n')
f.close()
# exit()

def draw_grouped_bar_chart(xlabel, Y, save_path, leg_pos = 'upper left'):
    labels = ['PauliSimp', 'UCCSynthesis', 'Paulihedral', 'POMCAS']
    n_bars = len(Y)
    colors = ['deepskyblue', 'darkviolet', 'black']
    # for i in [xlabel, Y, labels]:
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
        ax.bar(x - n_bars*width/2 + i*width, Y[i], width, label=labels[i], **styles[i])
        for j in range(len(x)):
            ax.text(x[j] - n_bars*width/2 + i*width, Y[i][j] + 2 * 10**4, '{:.2}'.format(Y[i][j]/10**6), fontsize=10, ha='center')
    # ax.set_xlabel('# qubits', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(xlabel, fontsize=12)
    ax.tick_params(axis='y', labelsize=12)
    ax.legend(loc=leg_pos, fontsize=12)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

# overview
Y = []
for cp in cps:
    y = []
    for metric in metrics:
        yi = 0
        for ac in acs:
            for bm in bms:
                yi += data[ac][cp][bm][metric]
        y.append(yi)
    Y.append(y.copy())

for i, metric in enumerate(metrics):
    print(metric + ':')
    y0 = min(Y[0][i], Y[1][i])
    print('r1: {:.0%}'.format(1 - Y[2][i]/y0))
    print('r2: {:.0%}'.format(1 - Y[3][i]/Y[2][i]))

path = './overview/overview.png'
draw_grouped_bar_chart(['CNOT count', 'depth'], Y, path, 'upper right')

# effect of arch
for metric in metrics:
    print('\n' + metric + ':')
    Y = []
    for cp in cps:
        y = []
        for ac in acs:
            yi = 0
            for bm in bms:
                yi += data[ac][cp][bm][metric]
            y.append(yi)
        Y.append(y.copy())

    for i, ac in enumerate(acs):
        print(ac + ':')
        y0 = min(Y[0][i], Y[1][i])
        print('r1: {:.0%}'.format(1 - Y[2][i]/y0))
        print('r2: {:.0%}'.format(1 - Y[3][i]/Y[2][i]))


    path = './overview/arch_{}.png'.format(metric)
    draw_grouped_bar_chart(acs, Y, path, 'upper right')

# effect of benchmark
for metric in metrics:
    print('\n' + metric + ':')
    Y = []
    for cp in cps:
        y = []
        for s in [uccsd_index, molecules, random_index]:
            yi = 0
            for ac in acs:
                for bm in s:
                    yi += data[ac][cp][bm][metric]
            y.append(yi)
        Y.append(y.copy())

    for i, bm in enumerate(['ucc', 'molecule', 'random']):
        print(bm + ':')
        y0 = min(Y[0][i], Y[1][i])
        print('r1: {:.0%}'.format(1 - Y[2][i]/y0))
        print('r2: {:.0%}'.format(1 - Y[3][i]/Y[2][i]))

    path = './overview/bm_{}.png'.format(metric)
    draw_grouped_bar_chart(['uccsd', 'molecule', 'random'], Y, path, 'upper right')