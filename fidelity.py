from compiler import Compiler
import pickle
from benchmark.mypauli import pauliString
from qiskit.providers.aer import QasmSimulator
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeMelbourne, FakeVigo
import numpy as np
from qiskit.opflow import PauliSumOp, StateFn
from qiskit.quantum_info import SparsePauliOp
import datetime
import warnings

def print_time():
    print('[' + (datetime.datetime.utcnow() + datetime.timedelta(hours=8)).strftime("%Y-%m-%d %H:%M:%S") + '] ', end="")

print_time()
print("START EXECUTION")
warnings.simplefilter('ignore', category=DeprecationWarning)
warnings.simplefilter('ignore', category=FutureWarning)

device_backend = FakeVigo()
noise_model = NoiseModel.from_backend(device_backend)
device = QasmSimulator.from_backend(device_backend)
coupling_map = device.configuration().coupling_map
n_qubits = device.configuration().n_qubits
print_time()
print('n_PQubit: ', n_qubits)
noise_model = NoiseModel.from_backend(device)
basis_gates = noise_model.basis_gates
ac = 'grid0104'
with open('./benchmark/H2.pickle', 'rb') as f:
    block_list = pickle.load(f)
blocks = []
for b in block_list:
    blocks.append([pauliString(i[0], coeff=(i[1].real + i[1].imag)) for i in b])
    # for i in b:
    #     print(i[0])
    # print('\n')
print_time()
print('n_block: ', len(blocks))
print_time()
print(block_list)
cp = Compiler(blocks, ac)
qc1 = cp.start('go', opt = 1)
qc2 = cp.start('ph')
qcs = [qc1, qc2]

from qiskit import transpile
from qiskit.qasm3 import dumps
from qiskit.quantum_info import Statevector
def qasm_to_qcis(qasm_str: str):
    qasm_lines = qasm_str.strip().split("\n")
    qcis_lines = []
    for line in qasm_lines:
        if line.startswith(("OPENQASM", "include", "input", "qubit")):
            continue
        line = line.upper()
        line = line.replace("[", "")
        line = line.replace("]", "")
        line = line.replace(",", "")
        line = line.replace(";", "")
        line = line.replace(")", "")
        line = line.replace("(", " ")
        cmd = line.strip().split(" ")
        length = len(cmd)
        if length == 2:
            qcis_lines.append("{} {}".format(cmd[0], cmd[1]))
        elif length == 3:
            if cmd[1] == "PI/2":
                if cmd[0] == "RX":
                    qcis_lines.append("X2P {}".format(cmd[2]))
                    continue
                elif cmd[0] == "RY":
                    qcis_lines.append("Y2P {}".format(cmd[2]))
                    continue
            elif cmd[1] == "-PI/2":
                if cmd[0] == "RX":
                    qcis_lines.append("X2M {}".format(cmd[2]))
                    continue
                elif cmd[0] == "RY":
                    qcis_lines.append("Y2M {}".format(cmd[2]))
                    continue
            if cmd[0] == "CZ":
                qcis_lines.append("{} {} {}".format(cmd[0], cmd[1], cmd[2]))
                continue
            alpha = np.pi
            if '-' in cmd[1]:
                alpha = -alpha
            if '/' in cmd[1]:
                alpha = alpha / int(cmd[1][cmd[1].find('/') + 1:])
            qcis_lines.append("{} {} {}".format(cmd[0], cmd[2], alpha))
        else:
            raise ValueError("QASM格式无法识别")
    return "\n".join(qcis_lines)
i = 0
for qc in qcs:
    qc = transpile(qc, basis_gates=['cz', 'rx', 'ry', 'rz', 'h'], optimization_level=3)
    _qc = qc.assign_parameters([np.pi])
    _qc.qasm(filename='data/qc{}_bv.qasm'.format(i))
    i += 1
    qasm_str = dumps(_qc)
    qcis_str = qasm_to_qcis(qasm_str)
    for j, k in zip(list(range(4)), [6, 0, 7, 1]):
        qcis_str = qcis_str.replace('Q' + str(j), 'Q' + str(k))
    for j in [6, 0, 7, 1]:
        qcis_str += '\nM Q{}'.format(j)
    output_path = "./data/qc" + str(i) + "_bv.txt"
    with open(output_path, "w") as f:
        f.write(qcis_str)
    state = Statevector.from_instruction(_qc)
    # print(state)
# exit()
qc1.draw('mpl', filename='./data/figs/qc1_bv.png')
qc2.draw('mpl', filename='./data/figs/qc2_bv.png')



import matplotlib.pyplot as plt
from qiskit import Aer
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit.algorithms import VQE, NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import SPSA

seed = 170
iterations = 30
backend = Aer.get_backend('aer_simulator')
x = [[] for i in range(3)]
y = [[] for i in range(3)]
algorithm_globals.random_seed = seed


def store_intermediate_result1(eval_count, parameters, mean, std, flag):
    # print_time()
    # print(eval_count, ' ', mean)
    counts1.append(eval_count)
    values1.append(mean)

def store_intermediate_result2(eval_count, parameters, mean, std, flag):
    # print_time()
    # print(eval_count, ' ', mean)
    counts2.append(eval_count)
    values2.append(mean)

def store_intermediate_result(eval_count, mean, qc_index):
    x[qc_index].append(eval_count)
    y[qc_index].append(mean)

def vqe(ansatz, op_list, qc_index, noise = True):
    # https://qiskit.org/documentation/stubs/qiskit.utils.QuantumInstance.html
    if noise:
        qi = QuantumInstance(backend=backend, seed_simulator=seed, seed_transpiler=seed,
                         coupling_map=coupling_map, noise_model=noise_model) # , initial_layout=list(range(n_qubits))
    else:
        qi = QuantumInstance(backend=backend, seed_simulator=seed, seed_transpiler=seed,
                         coupling_map=coupling_map)
    # npme = NumPyMinimumEigensolver()
    # result = npme.compute_minimum_eigenvalue(operator=op_list)
    # ref_value = result.eigenvalue.real
    # print_time()
    # print(f'Reference value: {ref_value:.5f}')
    if qc_index == 1:
        spsa = SPSA(maxiter=iterations)
    else:
        spsa = SPSA(maxiter=iterations)
    vqe = VQE(ansatz, optimizer=spsa, quantum_instance=qi, initial_point=np.array([0.5 * np.pi])) #, callback=store_intermediate_result1
    energy_evaluation, expectation = vqe.get_energy_evaluation(operator=op_list, return_expectation=True)
    print(energy_evaluation([np.pi]))
    # state = Statevector.from_instruction(ansatz.assign_parameters([np.pi]))
    # print(state)
    # for i in range(50):
    #     theta = 2 * np.pi * i / 50
    #     y = energy_evaluation(np.array([[theta]]))
    #     print(i, ' -> ', y)
    #     store_intermediate_result(i, y, qc_index)
    return
    



ham1 = pickle.load(open('./benchmark/H2_hamiltonian_bv.pickle', 'rb'))
energy1 = 0
m = {'Z':-1,'I':1}
for p in ham1:
    if 'X' in p[0] or 'Y' in p[0]:
        continue
    else:
        energy1 += m[p[0][2]] * p[1]
ham2 = pickle.load(open('./benchmark/H2_hamiltonian.pickle', 'rb'))[15:20]
energy2 = 0
for p in ham2:
    if 'X' in p[0] or 'Y' in p[0]:
        continue
    else:
        energy2 += m[p[0][0]] * m[p[0][1]] * p[1]
print(energy1, ' ', energy2)

ham = pickle.load(open('./benchmark/H2_hamiltonian.pickle', 'rb'))[15:20]
print(ham)

pauli_map = {}
for p, v in ham:
    if p in pauli_map.keys():
        pauli_map[p] += v
    else:
        pauli_map[p] = v
observable = sum([PauliSumOp(SparsePauliOp(p, v), coeff=1) for p, v in pauli_map.items()])

# def convert_pi(pi):
#     res = {}
#     for k, v in pi.items():
#         if v < 6:
#             res[k] = v + 1
#         else:
#             res[k] = v + 2
#     return res
# vqe(qc1, observable, 0)
# vqe(qc2, observable, 1)
vqe(qc1, observable, 2, False)
exit()

f = open('./data/H2_data_bv.txt', 'w')
# for v in values1:
#     f.write(str(v) + ' ')
# f.write('\n')
# for v in values2:
#     f.write(str(v) + ' ')
for b in y:
    for bi in b:
        f.write(str(bi) + ' ')
    f.write('\n')
exit()
plt.clf()
plt.rcParams['figure.figsize'] = (12, 4)
for a, b in zip(x, y):
    plt.plot(a, b)
plt.xlabel(r'$\theta$', fontsize=12)
plt.xticks(ticks = [0, 24, 49], labels = [r'$0$', r'$\pi/2$', r'$\pi$'], fontsize=12)
plt.ylabel(r'$\langle \Psi (\theta)|H|\Psi (\theta)\rangle$', fontsize=12)
plt.title('Convergence with noise')
plt.savefig('./data/figs/vn_bv.png')