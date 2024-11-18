import matplotlib.pyplot as plt
x = range(50)
Y = []
f = open('../data/H2_data.txt', 'r')
for line in f:
    Y.append([float(j) for j in line.split()])
plt.clf()
plt.rcParams['figure.figsize'] = (12, 4)
labels = ['PauliForest' + r'($loss_1$)', 'Paulihedral' + r'($loss_2$)', 'Reference' + r'($loss_0$)']
linestyles = ['-', '--', '-.']
markers = ['s', '^', 'o']
for i, y in enumerate(Y):
    # plt.scatter()
    plt.plot(x, y, label = labels[i], linestyle=linestyles[i], markersize=5, color = 'black', marker=markers[i]) # , marker=markers[i]
plt.xlabel(r'$\theta$', fontsize=12)
plt.xticks(ticks = [0, 24, 49], labels = [r'$0$', r'$\pi$', r'$2\pi$'], fontsize=12)
plt.ylabel(r'$loss_i=\langle \Psi_i (\theta)|H|\Psi_i (\theta)\rangle$ [Hartree]', fontsize=12)
# plt.title('Energy functions')
plt.legend()
plt.savefig('../data/figs/vn.png')