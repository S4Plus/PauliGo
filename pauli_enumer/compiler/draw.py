import matplotlib.pyplot as plt

f = open('./data.txt', 'r')
y1 = [float(y) for y in f.readline().split()]
y2 = [float(y) for y in f.readline().split()]

f = open('./data_ref.txt', 'r')
y3 = [float(y) for y in f.readline().split()]
y4 = [float(y) for y in f.readline().split()]
ref = [-1.10115 for i in range(125)]

x = [x for x in range(1, 126)]

plt.rcParams['figure.figsize'] = (8, 4)
plt.plot(x, y1)
plt.plot(x, y2)
plt.plot(x, y3)
plt.plot(x, y4)
plt.plot(x, ref)
plt.xlabel('Iteration count')
plt.ylabel('Energy')
plt.title('Convergence of estimate of H2 energy')
plt.legend(['Ours+noise', 'Paulihedral+noise', 'Ours', 'Paulihedral', 'Reference'])
# plt.savefig('./vn.png')
plt.savefig('./_vn.png')