import os
os.environ["OMP_NUM_THREADS"] = "1"
from VQE import vqe
from ansatz_source import mole_ansatz
from reference_state import single_reference_state
# from function import *
from hamiltonian_source import init_scf
# from scipy.linalg import *
# from sparse_matrices_helper import *
from geometry import mole
n_procs=8

geometry=mole(d=1).LiH()
molecule, hamiltonian_fermOp, n_qubits, occ_indices=init_scf(geometry, basis="sto-3g")
sys_information=[molecule.nelectron, hamiltonian_fermOp, n_qubits, occ_indices]
ref_states=single_reference_state(sys_information).HF_state()

op_pool=mole_ansatz(sys_information).UCCGSD()
terminal=vqe(op_pool, ref_states, sys_information, adapt_tol = 0.04, adapt_maxiter=100, n_procs=n_procs)

energy_i, psi_adapt_array, new_operator_pool_idx, new_parameters, error = terminal.adapt_vqe(add_para_num=1)

from compiler.mypauli import pauliString
blocks = []
for i in new_operator_pool_idx:
    op = op_pool[1][i]
    b = []
    for k, v in op.terms.items():
        P = ['I'] * n_qubits
        for ik in k:
            P[ik[0]] = ik[1]
        ps = ''.join(P)
        b.append([ps, v])
    blocks.append(b.copy())
print(blocks[:2])
import pickle
with open('./compiler/benchmark/LiH.pickle', 'wb') as f:
    pickle.dump(blocks, f)

from openfermion import jordan_wigner, bravyi_kitaev
qubit_hamiltonian = [jordan_wigner(op) for op in hamiltonian_fermOp]
qubit_hamiltonian_terms = []
for op in qubit_hamiltonian:
    for k, v in op.terms.items():
        P = ['I'] * n_qubits
        for ik in k:
            P[ik[0]] = ik[1]
        ps = ''.join(P)
        qubit_hamiltonian_terms.append([ps, v])
with open('./compiler/benchmark/LiH_hamiltonian.pickle', 'wb') as f:
    pickle.dump(qubit_hamiltonian_terms, f)