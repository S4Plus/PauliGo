import openfermion
import numpy

def Pij(i: int, j: int):
    ia = i * 2 + 0
    ib = i * 2 + 1
    ja = j * 2 + 0
    jb = j * 2 + 1
    term1 = openfermion.FermionOperator(
        ((ja, 0), (ib, 0)),
        1.0
    )
    term2 = openfermion.FermionOperator(
        ((ia, 0), (jb, 0)),
        1.0
    )
    return numpy.sqrt(0.5) * (term1 + term2)


def Pij_dagger(i: int, j: int):
    return openfermion.hermitian_conjugated(Pij(i, j))


def Qij_plus(i: int, j: int):
    ia = i * 2 + 0
    ib = i * 2 + 1
    ja = j * 2 + 0
    jb = j * 2 + 1
    term = openfermion.FermionOperator(
        ((ja, 0), (ia, 0)),
        1.0
    )
    return term


def Qij_minus(i: int, j: int):
    ia = i * 2 + 0
    ib = i * 2 + 1
    ja = j * 2 + 0
    jb = j * 2 + 1
    term = openfermion.FermionOperator(
        ((jb, 0), (ib, 0)),
        1.0
    )
    return term


def Qij_0(i: int, j: int):
    ia = i * 2 + 0
    ib = i * 2 + 1
    ja = j * 2 + 0
    jb = j * 2 + 1
    term1 = openfermion.FermionOperator(
        ((ja, 0), (ib, 0)),
        1.0
    )
    term2 = openfermion.FermionOperator(
        ((ia, 0), (jb, 0)),
        1.0
    )
    return numpy.sqrt(0.5) * (term1 - term2)


def Qij_vec(i: int, j: int):
    return [Qij_plus(i, j), Qij_minus(i, j), Qij_0(i, j)]


def Qij_vec_dagger(i: int, j: int):
    return [openfermion.hermitian_conjugated(i) for i in Qij_vec(i, j)]


def Qij_vec_inner(a: int, b: int, i: int, j: int):
    vec_dagger = Qij_vec_dagger(a, b)
    vec = Qij_vec(i, j)
    return sum([vec[i] * vec_dagger[i] for i in range(len(vec))])






