import functools
import numpy
import pyscf
# import pyscf.lo
# import pyscf.symm
# import pyscf.cc
# import pyscf.fci
from openfermion import FermionOperator, normal_ordered
from function import phy_operator

def get_mo_integrals_from_molecule_and_hf_orb(mol: pyscf.gto.Mole, mo_coeff: numpy.ndarray, debug: bool = False):
    """
    Notes:
        For two_body_mo, the current order is (ps|qr). A transpose like
        numpy.moveaxis(two_body_mo, [0, 2, 3, 1], [0, 1, 2, 3]) is necessary
        to get the (p+, q+, r, s) order or [p, q, r, s] indexing.
    """

    hcore = mol.intor("int1e_nuc") + mol.intor("int1e_kin")
    one_body_mo = functools.reduce(numpy.dot, (mo_coeff.T, hcore, mo_coeff))
    eri = mol.intor("int2e")
    two_body_mo = None
    try:
        import opt_einsum
        two_body_mo = opt_einsum.contract(
            "ijkl, ip, js, kq, lr->psqr",
            eri, mo_coeff.conj(), mo_coeff,
            mo_coeff.conj(), mo_coeff)
    except ModuleNotFoundError:
        try:
            two_body_mo = numpy.einsum(
                "ijkl, ip, js, kq, lr->psqr",
                eri, mo_coeff.conj(), mo_coeff,
                mo_coeff.conj(), mo_coeff)
        except ValueError:
            print(
                "The system is too large. Please install opt_einsum and re-run init_scf().")
            raise ValueError

    if debug:
        two_body_mo_check = pyscf.ao2mo.restore(1, pyscf.ao2mo.get_mo_eri(
            mol, mo_coeff, compact=False),
            mol.nao_nr()
        )
        error = numpy.linalg.norm(two_body_mo - two_body_mo_check)
        assert(numpy.isclose(error, 0.0))
    return one_body_mo, two_body_mo


def inter_pqrs(one_body_mo, two_body_mo, eps):
    two_body_mo_pqrs = numpy.moveaxis(two_body_mo, [0, 2, 3, 1], [0, 1, 2, 3])
    pq_tuple=zip(*((abs(one_body_mo) > eps).nonzero()))
    pqrs_tuple=zip(*((abs(two_body_mo_pqrs) > eps).nonzero()))
    return one_body_mo, two_body_mo_pqrs, pq_tuple, pqrs_tuple


def get_hamiltonian_ferm_op_from_mo_ints(one_body_mo: numpy.ndarray, two_body_mo: numpy.ndarray, eps: float = 0.0):
    """
    Construct the one- and two-body terms of the Hamiltonian for a given
    one-electron MO integral in (p+, q) order and a give two-electron MO
    integral in (p+, s, q+, r) order (PySCF's ordering).

    Args:
        one_body_mo (numpy.ndarray): one-electron integral in (p+, q) order.
        two_body_mo (numpy.ndarray): two-electron integral in
            (p+, s, q+, r) order.
        eps (float): cut-off threshold.

    Notes:
        The integrals are for spin-orbitals.
    """
    two_body_mo_pqrs = numpy.moveaxis(two_body_mo, [0, 2, 3, 1], [0, 1, 2, 3])
    hamiltonian_ferm_op_1 =FermionOperator()
    hamiltonian_ferm_op_2 =FermionOperator()
    for (p, q) in zip(*((abs(one_body_mo) > eps).nonzero())):
        p = int(p)
        q = int(q)
        pa = p * 2
        pb = p * 2 + 1
        qa = q * 2
        qb = q * 2 + 1
        hamiltonian_ferm_op_1 +=FermionOperator(((pa, 1), (qa, 0)), one_body_mo[p][q])
        hamiltonian_ferm_op_1 +=FermionOperator(((pb, 1), (qb, 0)), one_body_mo[p][q])
    for (p, q, r, s) in zip(*((abs(two_body_mo_pqrs) > eps).nonzero())):
        p = int(p)
        q = int(q)
        r = int(r)
        s = int(s)
        pa = p * 2
        pb = p * 2 + 1
        qa = q * 2
        qb = q * 2 + 1
        ra = r * 2
        rb = r * 2 + 1
        sa = s * 2
        sb = s * 2 + 1
        hamiltonian_ferm_op_2 +=FermionOperator(
            ((pa, 1), (qa, 1), (ra, 0), (sa, 0)),
            two_body_mo_pqrs[p][q][r][s] * 0.5
        )
        hamiltonian_ferm_op_2 +=FermionOperator(
            ((pb, 1), (qb, 1), (rb, 0), (sb, 0)),
            two_body_mo_pqrs[p][q][r][s] * 0.5
        )
        hamiltonian_ferm_op_2 +=FermionOperator(
            ((pa, 1), (qb, 1), (rb, 0), (sa, 0)),
            two_body_mo_pqrs[p][q][r][s] * 0.5
        )
        hamiltonian_ferm_op_2 +=FermionOperator(
            ((pb, 1), (qa, 1), (ra, 0), (sb, 0)),
            two_body_mo_pqrs[p][q][r][s] * 0.5
        )

    hamiltonian_ferm_op_1 = normal_ordered(hamiltonian_ferm_op_1)
    hamiltonian_ferm_op_2 = normal_ordered(hamiltonian_ferm_op_2)
    return hamiltonian_ferm_op_1, hamiltonian_ferm_op_2


def init_scf(geometry, basis="sto-6g", spin=0, run_fci: bool = False, run_rccsd: bool = False):
    eps = 0.0
    molecule = pyscf.gto.M(atom=geometry, basis=basis, spin=spin)
    mf = pyscf.scf.RHF(molecule)
    print("Running RHF...")
    mf.kernel()
    mo_coeff = mf.mo_coeff

    if run_rccsd:
        print("Running RCCSD")
        mf_cc = pyscf.cc.RCCSD(mf)
        mf_cc.kernel()
        energy_RCCSD = mf_cc.e_tot
        print("CCSD energy: %20.16f Ha" % (energy_RCCSD))

    if run_fci:
        mf_fci = pyscf.fci.FCI(mf)
        energy_fci = mf_fci.kernel()[0]
        print("FCI energy: %20.16f Ha" % (energy_fci))

    energy_nuc = molecule.energy_nuc()
    n_orb = mo_coeff.shape[1]
    n_electrons=molecule.nelectron 
    one_body_mo, two_body_mo = get_mo_integrals_from_molecule_and_hf_orb(molecule, mo_coeff)
    core_correction = 0.0

    particle_num=n_electrons
    hamiltonian_ferm_op_1, hamiltonian_ferm_op_2 = get_hamiltonian_ferm_op_from_mo_ints(one_body_mo, two_body_mo, eps)
    n_qubits=n_orb*2

    s2_fermion=phy_operator(n_qubits)[0][2]
    occ_indices = [i for i in range(particle_num)]
    hamiltonian_ferm_op = hamiltonian_ferm_op_1 + hamiltonian_ferm_op_2
    hamiltonian_ferm_op += energy_nuc + core_correction+ 0.5*s2_fermion

    head="{0:<26}{1:^10}"
    table="{0:<1}{1:<30}{2:^1}{3:<30}{4:>1}"
    print(head.format(" ",  "Molecule Information"))
    print("+-------------------------------------------------------------+")
    print(table.format("|", "Qubits",           "|",    str(n_qubits),             "|"))
    print(table.format("|", "Electrons",        "|",    str(particle_num),       "|"))
    print("+-------------------------------------------------------------+")
    print('\n')

    return molecule, hamiltonian_ferm_op, n_qubits, occ_indices


