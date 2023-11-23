from compiler import PH_compiler, My_compiler
from arch import load_graph
import pickle
from mypauli import pauliString

ac = load_graph('grid0403')
with open('./benchmark/LiH.pickle', 'rb') as f:
    block_list = pickle.load(f)
blocks = []
for b in block_list:
    blocks.append([pauliString(i[0], coeff=(i[1].real + i[1].imag)) for i in b])
print('n_block: ', len(blocks))
cp1 = My_compiler(blocks, ac)
cp2 = PH_compiler(blocks, ac)
res1, qc1 = cp1.compile()
res2, qc2 = cp2.compile()
from qiskit import transpile
qc1 = transpile(qc1, basis_gates=['cx', 'u', 'rz'], optimization_level=3)
qc2 = transpile(qc2, basis_gates=['cx', 'u', 'rz'], optimization_level=3)

from qiskit.circuit import Parameter
def store_cir(qc, path):
    f = open(path, mode='w')
    for op in qc.data:
        f.write(op.operation.name + ' ')
        for q in op.qubits:
            f.write('q[{0}] '.format(q.index))
        if len(op.operation.params) == 0:
            f.write('\n')
            continue
        if isinstance(op.operation.params[0], Parameter):
            f.write(op.operation.params[0].name)
        else:
            f.write(str(op.operation.params[0]))
        f.write('\n')
    f.close()

# 存储线路
store_cir(qc1, './qc1.qasm')
store_cir(qc2, './qc2.qasm')

with open('./benchmark/LiH_hamiltonian.pickle', 'rb') as f:
    qubit_op = pickle.load(f)

from pyqpanda import *
from qp_ex import prepareInitialState

# get molcule hamiltonian
def get_mol_pauli(qubit_op):
    mol_pauli = []
    for it in qubit_op:
        mol_pauli.append(PauliOperator(it[0], it[1]))
    return mol_pauli

# qiskit QuantumCircuit -> qpanda VariationalQuantumCircuit
def get_ansatz(qlist, qiskit_cir):
    vqc = VariationalQuantumCircuit()
    var_para = []
    for i in range(qiskit_cir.num_parameters):
        var_para.append(var(0.5, True))
    para_i = {}
    index = 0
    for para in qiskit_cir.parameters:
        para_i[para] = index
        index += 1
    for op in qiskit_cir.data:
        if op.operation.name == 'u':
            qi = op.qubits[0].index
            params = op.operation.params
            vqc.insert(U3(qlist[qi]), params[0], params[1], params[2])
        if op.operation.name == 'rz':
            qi = op.qubits[0].index
            para = op.operation.params[0]
            if para in para_i.keys():
                para = var_para[para_i[para]]
                vqc.insert(VariationalQuantumGate_RZ(qlist[qi]), para)
            else:
                vqc.insert(RZ(qlist[qi]), para)
        if op.operation.name == 'cx':
            qi = [q.index for q in op.qubits]
            vqc.insert(CNOT(qlist[qi[0]], qlist[qi[1]]))
    return vqc

def GradientDescent(mol_pauli, qiskit_cir, n_qubit, n_en, iters):
    machine=init_quantum_machine(QMachineType.CPU)
    qlist = machine.qAlloc_many(n_qubit)

    vqc = VariationalQuantumCircuit()
    vqc.insert(prepareInitialState(qlist, n_en))
    vqc.insert(get_ansatz(qlist, qiskit_cir))
    convert_qprog_to_originir(vqc, machine)

    loss = qop(vqc, mol_pauli, machine, qlist)
    gd_optimizer = MomentumOptimizer.minimize(loss, 0.1, 0.9)
    leaves = gd_optimizer.get_variables()

    min_energy=float('inf')
    for i in range(iters):
        gd_optimizer.run(leaves, 0)
        loss_value = gd_optimizer.get_loss()
    
        print(loss_value)
        if loss_value < min_energy:
            min_energy = loss_value
    
    return min_energy
if __name__=="__main__": 
    mol_pauli = get_mol_pauli(qubit_op)
    GradientDescent(mol_pauli, qc1, 12, 4, 300)
    GradientDescent(mol_pauli, qc2, 12, 4, 300)