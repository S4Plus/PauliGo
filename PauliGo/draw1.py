import matplotlib.pyplot as plt
import numpy as np

def read_data_from_file(file_path, column):
    data = []
    if 'random' in file_path:
        column -= 1
    with open(file_path, 'r') as file:
        for line in file:
            item = line.strip().split()
            name = item[0]
            name = name[name.find('_') + 1:]
            data.append((name, int(item[column])))
    return data

def draw_grouped_bar_chart(data_go, data_ph, save_path):
    names = [item[0] for item in data_go]
    data_set1 = [item[1] for item in data_go]
    data_set2 = [item[1] for item in data_ph]

    x = np.arange(len(names))
    width = 0.25
    mean1 = int(np.mean(data_set1))
    mean2 = int(np.mean(data_set2))
    print(mean1, ' ', mean2, ' ', (mean2 - mean1) / mean2)
    fig, ax = plt.subplots(figsize=(6, 4))
    bars1 = ax.bar(x - width/2, data_set1, width, label='Ini+Syn', color='deepskyblue')
    bars2 = ax.bar(x + width/2, data_set2, width, label='Paulihedral', color='darkviolet')
    ax.axhline(mean1, color='deepskyblue', linestyle='dashed', linewidth=2)
    ax.axhline(mean2, color='darkviolet', linestyle='dashed', linewidth=2)
    ax.annotate('{}'.format(mean1),
                xy=(0, mean1),
                xytext=(0, -20),  # 3 points vertical offset
                textcoords="offset points",
                ha='center', va='bottom', fontsize=12, color='deepskyblue')
    ax.annotate('{}'.format(mean2),
                xy=(0, mean2),
                xytext=(0, 6),  # 3 points vertical offset
                textcoords="offset points",
                ha='center', va='bottom', fontsize=12, color='darkviolet')
    ax.set_xlabel('# qubits', fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=12)
    ax.legend(loc='upper left', fontsize=12)

    # Attach values to the bars
    # for bar1, bar2 in zip(bars1, bars2):
    #     height1 = bar1.get_height()
    #     height2 = bar2.get_height()
    #     ax.annotate('{}'.format(height1),
    #                 xy=(bar1.get_x() + bar1.get_width() / 2, height1),
    #                 xytext=(0, 3),  # 3 points vertical offset
    #                 textcoords="offset points",
    #                 ha='center', va='bottom')
    #     ax.annotate('{}'.format(height2),
    #                 xy=(bar2.get_x() + bar2.get_width() / 2, height2),
    #                 xytext=(0, 3),  # 3 points vertical offset
    #                 textcoords="offset points",
    #                 ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()

def draw_rate(data_go, data_ph, save_path):

    names = [item[0] for item in data_go]
    data_set1 = []
    for i in range(len(data_go)):
        if data_ph[i][1] > 0:
            data_set1.append(data_go[i][1]/data_ph[i][1])
        else:
            data_set1.append(data_go[i][1])

    x = np.arange(len(names))
    width = 0.5

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(x, data_set1, width, label='pg/ph', color='deepskyblue')
    ax.axhline(np.mean(data_set1), color='darkviolet', linestyle='dashed', linewidth=2)

    ax.set_xlabel('# qubits', fontsize=12)
    ax.set_ylim(0.1, 1.0)
    # ax.set_ylabel('pg/ph')
    # ax.set_title('Grouped Bar Chart')
    ax.set_xticks(x)
    ax.set_xticklabels(names)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.legend(loc='upper left', fontsize=12)

    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
if __name__ == "__main__":
    share_path = "./data/v3/"
    benchmarks = ['molecule', 'uccsd', 'random']
    hardwares = ['manhattan']  # 'manhattan',
    compilers = ['go', 'ph']
    columns = [['swap', -3], ['depth', -1]]
    save_path0 = "./data/figs/"  # Replace "path_to_folder" with the actual folder path where you want to save the picture
    for b in benchmarks:
        for h in hardwares:
            for column in columns:
                print(h, ' ', b, ' ', column[0])
                data = []
                for c in compilers:
                    data.append(read_data_from_file(share_path + b + '_' + h + '_' + c + '_opt2.txt', column[1]))
                save_path = save_path0 + b + '_' + h + '_' + column[0] + '.png'
                draw_grouped_bar_chart(data[0], data[1], save_path)
