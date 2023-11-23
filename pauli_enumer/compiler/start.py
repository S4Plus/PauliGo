from compiler import PH_compiler, My_compiler
from arch import load_graph
import pickle
from mypauli import pauliString
from functions import print_qc, fake_device_coupling_graph
from qiskit.providers.aer import QasmSimulator
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.fake_provider import FakeMelbourne, FakeVigo
import numpy as np
from qiskit.opflow import PauliSumOp, StateFn
from qiskit.quantum_info import SparsePauliOp

device_backend = FakeVigo()
noise_model = None
device = QasmSimulator.from_backend(device_backend)
coupling_map = device.configuration().coupling_map
n_qubits = device.configuration().n_qubits
print('n_PQubit: ', n_qubits)
noise_model = NoiseModel.from_backend(device)
basis_gates = noise_model.basis_gates
coupling_map = [[0,1],[1,2],[2,3]]
ac = fake_device_coupling_graph(coupling_map, 4)
with open('./benchmark/H2.pickle', 'rb') as f:
    block_list = pickle.load(f)
blocks = []
for b in block_list:
    blocks.append([pauliString(i[0], coeff=(i[1].real + i[1].imag)) for i in b])
    # for i in b:
    #     print(i[0])
    # print('\n')
print('n_block: ', len(blocks))
cp1 = My_compiler(blocks, ac)
cp2 = PH_compiler(blocks, ac)
res1, qc1 = cp1.compile()
res2, qc2 = cp2.compile()

print('our opt_level3: ', end="")
qc1 = print_qc(qc1, opt_level = 3)
print('ph opt_level3: ', end="")
qc2 = print_qc(qc2, opt_level = 3)
print('n_parameter: ', qc1.num_parameters)


import matplotlib.pyplot as plt
from qiskit import Aer
from qiskit.utils import QuantumInstance, algorithm_globals
from qiskit.algorithms import VQE, NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import SPSA

seed = 170
iterations = 125
backend = Aer.get_backend('aer_simulator')
counts1 = []
values1 = []
counts2 = []
values2 = []
algorithm_globals.random_seed = seed
# https://qiskit.org/documentation/stubs/qiskit.utils.QuantumInstance.html
qi = QuantumInstance(backend=backend, seed_simulator=seed, seed_transpiler=seed,
                     coupling_map=coupling_map, noise_model=noise_model) # , initial_layout=list(range(n_qubits))

def store_intermediate_result1(eval_count, parameters, mean, std, flag):
    print(eval_count, ' ', mean)
    counts1.append(eval_count)
    values1.append(mean)

def store_intermediate_result2(eval_count, parameters, mean, std, flag):
    print(eval_count, ' ', mean)
    counts2.append(eval_count)
    values2.append(mean)

def vqe(ansatz, op_list, qc_index):
    npme = NumPyMinimumEigensolver()
    result = npme.compute_minimum_eigenvalue(operator=op_list)
    ref_value = result.eigenvalue.real
    print(f'Reference value: {ref_value:.5f}')
    if qc_index == 1:
        spsa = SPSA(maxiter=iterations, callback=store_intermediate_result1)
    else:
        spsa = SPSA(maxiter=iterations, callback=store_intermediate_result2)
    vqe = VQE(ansatz, optimizer=spsa, quantum_instance=qi, initial_point=np.array([0.5])) #, callback=store_intermediate_result1
    # energy_evaluation, expectation = vqe.get_energy_evaluation(operator=op_list, return_expectation=True)
    # for i in range(50):
    #     theta = np.pi * i / 50
    #     print(i, ' -> ', energy_evaluation(np.array([[theta]])))
    # input()
    result1 = vqe.compute_minimum_eigenvalue(operator=op_list)
    print(f'VQE on Aer qasm simulator (no noise): {result1.eigenvalue.real:.5f}')
    print(f'Delta from reference energy value is {(result1.eigenvalue.real - ref_value):.5f}')
    


f = open('./benchmark/H2_hamiltonian.pickle', 'rb')
ham = pickle.load(f)
pauli_map = {}
for p, v in ham:
    if p in pauli_map.keys():
        pauli_map[p] += v
    else:
        pauli_map[p] = v
observable = sum([PauliSumOp(SparsePauliOp(p, v), coeff=1) for p, v in pauli_map.items()])

# pi = cp2.my_pi
# m_qubit_op = []
# for ps, v in pauli_map.items():
#     P = ['I'] * n_qubits
#     k = 0
#     for c in ps:
#         if c != 'I':
#             P[pi[k]] = c
#         k += 1
#     m_qubit_op.append(PauliSumOp(SparsePauliOp(''.join(P), coeffs=v), coeff=1))
# observable = sum(m_qubit_op)

state_obs = StateFn(observable, is_measurement=True)

def loss1(x):
    bound = qc1.assign_parameters(x)
    return np.real((state_obs @ StateFn(bound)).eval())

def loss2(x):
    bound = qc2.assign_parameters(x)
    return np.real((state_obs @ StateFn(bound)).eval())

# for i in range(50):
#     theta = 2 * np.pi * i / 50
#     param_dict = {qc2.parameters[0] : theta}
#     print(i, ' -> ', loss2(param_dict))
# input()
vqe(qc2, observable, 2)
vqe(qc1, observable, 1)

f = open('./data.txt', 'w')
for v in values1:
    f.write(str(v) + ' ')
f.write('\n')
for v in values2:
    f.write(str(v) + ' ')
plt.rcParams['figure.figsize'] = (12, 4)
plt.plot(counts1, values1)
plt.plot(counts2, values2)
plt.xlabel('Eval count')
plt.ylabel('Energy')
plt.title('Convergence with noise')
plt.savefig('./vn.png')