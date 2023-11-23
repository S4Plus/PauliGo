from openfermion import FermionOperator, QubitOperator
from spin_adapt_source import *

def SFEx_generator(i,j):
    ia = i * 2 + 0
    ib = i * 2 + 1
    ja = j * 2 + 0
    jb = j * 2 + 1
    term1 = FermionOperator(((ia, 1), (ja, 0)), 1.0)
    term2 = FermionOperator(((ib, 1), (jb, 0)), 1.0)
    tpq_list = [term1 , term2]
    return tpq_list


def DFEx_generator(p,q,r,s):
    """
    p < q is vir ; r < s is occ, r < p < s < q
    """
    pa = p * 2 + 0
    pb = p * 2 + 1
    qa = q * 2 + 0
    qb = q * 2 + 1
    ra = r * 2 + 0
    rb = r * 2 + 1
    sa = s * 2 + 0
    sb = s * 2 + 1
    
    # term1=FermionOperator(((sa,0), (pa,1), (ra,0), (qa,1)), 1.0)
    # term2=FermionOperator(((sb,0), (pb,1), (rb,0), (qb,1)), 1.0)
    # term3=FermionOperator(((sa,0), (pb,1), (rb,0), (qa,1)), 1.0)
    # term4=FermionOperator(((sb,0), (pa,1), (ra,0), (qb,1)), 1.0)

    # term1=FermionOperator(((sa,0), (pa,1), (ra,0), (qa,1)), 1.0)
    # term2=FermionOperator(((sb,0), (pb,1), (rb,0), (qb,1)), 1.0)
    # term3=FermionOperator(((sa,0), (pa,1), (rb,0), (qb,1)), 1.0)
    # term4=FermionOperator(((sb,0), (pb,1), (ra,0), (qa,1)), 1.0)

    # term1=FermionOperator(((sa,0), (pb,1), (rb,0), (qa,1)), 1.0)
    # term2=FermionOperator(((sb,0), (pa,1), (ra,0), (qb,1)), 1.0)
    # term3=FermionOperator(((sa,0), (pa,1), (rb,0), (qb,1)), 1.0)
    # term4=FermionOperator(((sb,0), (pb,1), (ra,0), (qa,1)), 1.0)


    term1=FermionOperator(((sa,0), (pa,1), (ra,0), (qa,1)), 1.0)
    term2=FermionOperator(((sb,0), (pb,1), (rb,0), (qb,1)), 1.0)
    term3=FermionOperator(((sa,0), (pa,1), (rb,0), (qb,1)), 1.0)
    term4=FermionOperator(((sb,0), (pb,1), (ra,0), (qa,1)), 1.0)
    term5=FermionOperator(((sa,0), (pb,1), (rb,0), (qa,1)), 1.0)
    term6=FermionOperator(((sb,0), (pa,1), (ra,0), (qb,1)), 1.0)
    # Uex_pqrs=[term1,term2,term3,term4]
    Uex_pqrs=[term1,term2,term3,term4,term5,term6]
    return Uex_pqrs


def SQEx_generator(i,j):
    term1=QubitOperator([(j,"Y"), (i, "X")], -0.5j)
    term2=QubitOperator([(j,"X"), (i, "Y")], 0.5j)
    single=term1+term2
    return single


def DQEx_generator(i,j,k,l):
    term1=QubitOperator([(i,"X"), (j, "Y"), (k, "X"), (l, "X")], 0.125j)
    term2=QubitOperator([(i,"Y"), (j, "X"), (k, "X"), (l, "X")], 0.125j)
    term3=QubitOperator([(i,"Y"), (j, "Y"), (k, "Y"), (l, "X")], 0.125j)
    term4=QubitOperator([(i,"Y"), (j, "Y"), (k, "X"), (l, "Y")], 0.125j)
    term5=QubitOperator([(i,"X"), (j, "X"), (k, "Y"), (l, "X")], -0.125j)
    term6=QubitOperator([(i,"X"), (j, "X"), (k, "X"), (l, "Y")], -0.125j)
    term7=QubitOperator([(i,"Y"), (j, "X"), (k, "Y"), (l, "Y")], -0.125j)
    term8=QubitOperator([(i,"X"), (j, "Y"), (k, "Y"), (l, "Y")], -0.125j)
    double=term1+term2+term3+term4+term5+term6+term7+term8
    return double


def SPEx_generator(i,j):
    term1=QubitOperator([(j,"Y"), (i, "X")], 1.j)
    term2=QubitOperator([(j,"X"), (i, "Y")], 1.j)
    single=[term1,term2]
    return single


def DPEx_generator(i,j,k,l):
    term1=QubitOperator([(i,"X"), (j, "Y"), (k, "X"), (l, "X")], 1.j)
    term2=QubitOperator([(i,"Y"), (j, "X"), (k, "X"), (l, "X")], 1.j)
    term3=QubitOperator([(i,"Y"), (j, "Y"), (k, "Y"), (l, "X")], 1.j)
    term4=QubitOperator([(i,"Y"), (j, "Y"), (k, "X"), (l, "Y")], 1.j)
    term5=QubitOperator([(i,"X"), (j, "X"), (k, "Y"), (l, "X")], 1.j)
    term6=QubitOperator([(i,"X"), (j, "X"), (k, "X"), (l, "Y")], 1.j)
    term7=QubitOperator([(i,"Y"), (j, "X"), (k, "Y"), (l, "Y")], 1.j)
    term8=QubitOperator([(i,"X"), (j, "Y"), (k, "Y"), (l, "Y")], 1.j)
    double=[term1,term2,term3,term4,term5,term6,term7,term8]
    return double


def SA_SFEx_generator(i, j):
    """
    Spin-adapted single excitation operators.

    Args:
        i (int): index of the spatial orbital which the
            creation operator will act on.
        j (int): index of the spatial orbital which the
            annihilation operator will act on.

    Returns:
        tpq_list (list): Spin-adapted single excitation operators.

    Reference:
        Scuseria, G. E. et al., J. Chem. Phys. 89, 7382 (1988)
    """
    ia = i * 2 + 0
    ib = i * 2 + 1
    ja = j * 2 + 0
    jb = j * 2 + 1
    term1 = FermionOperator(((ia, 1), (ja, 0)), 1.0)
    term2 = FermionOperator(((ib, 1), (jb, 0)), 1.0)
    tpq_list = [term1 + term2]
    return tpq_list


def SA_DFEx_generator(creation_list, annihilation_list):
    """
    Spin-adapted double excitation operators.

    Args:
        creation_list (list): list of spatial orbital indices which the
            creation operator will act on.
        annihilation_list (list): list of spatial orbital indices which the
            annihilation operator will act on.

    Returns:
        tpqrs_list (list): Spin-adapted double excitation operators.

    Reference:
        Igor O. Sokolov et al., J. Chem. Phys. 152, 124107 (2020)
        Ireneusz W. Bulik et al., J. Chem. Theory Comput. 11, 3171-3179 (2015)
        Scuseria, G. E. et al., J. Chem. Phys. 89, 7382 (1988)
    """
    p = creation_list[0]
    r = annihilation_list[0]
    q = creation_list[1]
    s = annihilation_list[1]
    tpqrs1 = Pij_dagger(p, q) * Pij(r, s)
    tpqrs2 = Qij_vec_inner(p, q, r, s)
    tpqrs_list = [tpqrs1, tpqrs2]
    return tpqrs_list


