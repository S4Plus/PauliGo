import numpy as np
from qiskit.algorithms.optimizers import SPSA
from qiskit.circuit.library import PauliTwoDesign
from qiskit.opflow import Z, StateFn, PauliSumOp
from qiskit.quantum_info import SparsePauliOp
import pickle
from compiler import PH_compiler, My_compiler
from arch import load_graph
from mypauli import pauliString
from functions import print_qc, fake_device_coupling_graph
from qiskit.algorithms import NumPyMinimumEigensolver
from qiskit.providers.fake_provider import FakeVigo

coupling_map = [[0,1],[1,2],[2,3]]
# ac = load_graph('all04')
ac = fake_device_coupling_graph(coupling_map, 4)
with open('./benchmark/H2.pickle', 'rb') as f:
    block_list = pickle.load(f)
blocks = []
for b in block_list:
    blocks.append([pauliString(i[0], coeff=(i[1].real + i[1].imag)) for i in b])

print('n_block: ', len(blocks))
cp1 = My_compiler(blocks, ac)
cp2 = PH_compiler(blocks, ac)
res1, qc1 = cp1.compile(opt=1)
res2, qc2 = cp2.compile()

print('our opt_level3: ', end="")
qc1 = print_qc(qc1, opt_level = 3)
print('ph opt_level3: ', end="")
qc2 = print_qc(qc2, opt_level = 3)
print('n_parameter: ', qc1.num_parameters)

f = open('./benchmark/H2_hamiltonian.pickle', 'rb')
ham = pickle.load(f)
pauli_map = {}
for p, v in ham:
    if p in pauli_map.keys():
        pauli_map[p] += v
    else:
        pauli_map[p] = v
# observable = SparsePauliOp.from_list([(p, v) for p, v in pauli_map.items()])
observable = sum([PauliSumOp(SparsePauliOp(p, v), coeff=1) for p, v in pauli_map.items()])
state_obs = StateFn(observable, is_measurement=True)
print('number of paulis: ', len(pauli_map))
npme = NumPyMinimumEigensolver()
result = npme.compute_minimum_eigenvalue(operator=observable)
ref_value = result.eigenvalue.real
print(f'Reference value: {ref_value:.5f}')

ansatz1 = qc1
ansatz2 = qc2
# ansatz = PauliTwoDesign(4, reps=1, seed=2)
# observable = Z ^ Z
initial_point = np.random.random(ansatz1.num_parameters)

def store_intermediate_result(eval_count, parameters, mean, std, flag):
    print(eval_count, ' ', mean)

def loss1(x):
    bound = ansatz1.assign_parameters(x)
    energy = np.real((state_obs @ StateFn(bound)).eval())
    # print(x, energy)
    return energy

def loss2(x):
    bound = ansatz2.assign_parameters(x)
    energy = np.real((state_obs @ StateFn(bound)).eval())
    return energy

# for i in range(50):
#     theta = np.pi * i / 50
#     param_dict = {ansatz1.parameters[0] : theta}
#     print(i, ' -> ', loss1(param_dict))

spsa = SPSA(maxiter=125, callback=store_intermediate_result)
result1 = spsa.optimize(ansatz1.num_parameters, loss1, initial_point=initial_point)
# result2 = spsa.optimize(ansatz2.num_parameters, loss2, initial_point=initial_point)
print(result1)
# print(result2)