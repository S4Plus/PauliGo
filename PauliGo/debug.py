from qiskit import QuantumCircuit, transpile
# qc = QuantumCircuit.from_qasm_file('./data/phycir/ph/uccsd_8_manhattan_opt2_origin.txt')
qc = QuantumCircuit.from_qasm_file('./data/phycir/go/debug.txt')
qc3 = transpile(qc, basis_gates=['cx', 'u3'], optimization_level=3)