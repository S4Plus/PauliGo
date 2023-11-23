import pdflatex
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.tools.visualization import circuit_drawer
q = QuantumRegister(3, name='q')
circuit = QuantumCircuit(q)
circuit.x(q[1])
circuit.cx(q[0], q[1])
circuit_drawer(circuit, output='latex')