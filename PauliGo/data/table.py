
if __name__ == "__main__":
    share_path = "./v3/"
    benchmarks = ['uccsd', 'molecule', 'random']
    hardwares = ['sycamore', 'manhattan']
    compilers = ['go', 'ph']
    columns = [1, 5, 7]
    tab = open('./tab.txt', mode='w')
    tab1 = open('./tab1.txt', mode='w')
    for h in hardwares:
        xi = 0
        d1 = [[0, 0, 0], [0, 0, 0]]
        d2 = [[0, 0, 0], [0, 0, 0]]
        for b in benchmarks:
            with open(share_path + b + '_' + h + '_' + compilers[0] + '_opt2.txt') as f1,\
                    open(share_path + b + '_' + h + '_' + compilers[1] + '_opt2.txt') as f2:
                for l1, l2 in zip(f1, f2):
                    it1 = l1.split()
                    it2 = l2.split()
                    #tab.write('& ')
                    name = it1[0].replace('_', '-')
                    nq = name[name.find('-')+1:]
                    if b == 'molecule':
                        name = name[:name.find('-')]
                    tab.write(name + ' & ' + nq)
                    for c in columns:
                        tab.write(' & ' + it1[c])
                    for c in columns[1:]:
                        tab.write(' & ' + it2[c])
                    tab.write('\\\\\n\\hline\n')
                    d1[0][xi] += int(it1[5])
                    d1[1][xi] += int(it1[7])
                    d2[0][xi] += int(it2[5])
                    d2[1][xi] += int(it2[7])
            xi += 1
        metrics = ['SWAP', 'depth']
        tab1.write('\multirow{2}{*}{Sycamore}')
        for i in range(2):
            tab1.write(' & ' + metrics[i])
            for j in range(3):
                tab1.write(' & ' + str(int((d2[i][j] - d1[i][j]) / d2[i][j] * 100)) + '\\%')
            if i == 1:
                tab1.write('\\\\\n\\hline\n')
            else:
                tab1.write('\\\\\n\\cline{2-5}\n')
        tab1.write('\n')
    tab1.close()