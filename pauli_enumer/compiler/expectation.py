from qiskit.opflow import StateFn, PauliExpectation, CircuitSampler
from qiskit import Aer, QuantumCircuit
from qiskit.opflow import CircuitStateFn, PauliSumOp
from qiskit.quantum_info import SparsePauliOp
from time import time
import pickle

f = open('./benchmark/LiH_hamiltonian.pickle', 'rb')
ham = pickle.load(f)
pauli_map = {}
for p, v in ham:
    if p in pauli_map.keys():
        pauli_map[p] += v
    else:
        pauli_map[p] = v
op = SparsePauliOp.from_list([(p, v) for p, v in pauli_map.items()])

state = QuantumCircuit(12)
for i in range(2):
    state.x(i)
state.x(10)
state.x(11)
opflow_op = PauliSumOp(op)
opflow_state = CircuitStateFn(state)

# hamiltonian info
print('number of paulis: ', len(pauli_map))
t1 = time()
ham_matix = opflow_op.to_matrix()
print('time cost for convert to matrix: ', time() - t1)

t1 = time()
# Define the state to sample
import numpy as np
import spsa
def energy(H : np.ndarray, phi : np.matrix):
    return np.array(phi @ H @ phi.H)[0][0].real

iter = 0
def loss1(theta : np.ndarray) -> float:
    return np.sum(theta) ** 2
    print(theta[0], end="")
    x = [complex(0 ,0)] * 4096
    x[15] = theta[0]
    x[3 + 2**11 + 2**10] = np.sqrt(1 - theta[0])
    phi = np.matrix([x])
    e =  energy(ham_matix, phi)
    print(' -> ', e)
    return e

min_energy = spsa.maximize(loss1, [1])

