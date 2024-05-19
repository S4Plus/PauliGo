# from qiskit import *
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import BravyiKitaevMapper
from benchmark.mypauli import pauliString

H2 = "H 0 0 0; H 0 0 0.735"
HF = "F 0.000000    0.000000    0.092693; H 0.000000    0.000000   -0.834240"
LiH = "Li 0.523344   -0.189189    0.000000; H -1.080371   -0.189189    0.000000"
H2O = "O -0.000000   -0.000000    0.120915; H -0.000000    0.756655   -0.483660; H -0.000000   -0.756655   -0.483660"
CH4 = "C 0.000000    0.000000    0.000000; H 0.635186    0.635186    0.635186; H -0.635186   -0.635186    0.635186; H -0.635186    0.635186   -0.635186; H 0.635186   -0.635186   -0.635186"
NH3 = "N -0.000000    0.000000    0.126158; H 0.000000    0.934569   -0.294370; H -0.809360   -0.467284   -0.294370; H 0.809360   -0.467284   -0.294370"
NH2 = "N -0.000000    0.000000    0.147085; H 0.000000    0.802314   -0.514796; H -0.000000   -0.802314   -0.514796"
CH2 = "C 0.000000    0.000000    0.180663; H -0.000000   -0.864762   -0.541989; H 0.000000    0.864762   -0.541989"
molecules = [['H2', H2], ['HF', HF], ['LiH', LiH], ['H2O', H2O], ['NH3', NH3], ['CH4', CH4]]  # , ['NH2', NH2], ['CH2', CH2]

def molecule_oplist(mole):
    driver = PySCFDriver(
        atom=mole,
        basis="sto3g",
        # charge=0,
        # spin=0,
        unit=DistanceUnit.ANGSTROM,
    )

    problem = driver.run()
    fermionic_op = problem.hamiltonian.second_q_op()
    mapper = BravyiKitaevMapper()
    qubit_jw_op = mapper.map(fermionic_op)
    oplist = []
    for op in qubit_jw_op:
        ps = str(op.primitive.paulis[0])
        if all(p == 'I' for p in ps):
            continue
        oplist.append([pauliString(ps, coeff=op.coeffs)])
    return oplist
# paulis = qubit_jw_op.primitive.paulis
# coeffs = qubit_jw_op.primitive.coeffs
# print(qubit_jw_op)
import pickle
for mole in molecules:
    op_list = molecule_oplist(mole[1])
    f = open('./benchmark/data/{}.pickle'.format(mole[0]), '+wb')
    pickle.dump(op_list, f)