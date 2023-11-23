import openfermion
import numpy
import scipy
import math

def create_initial_hartree_fock_state(
        occ_indices_spin: list,
        n_qubits: int):
    r"""
    Create the Hartree-Fock reference state using JW transformation:
        |HF>=\prod_{i} {a^{\dagger}_{i}} |0>

    Args:
        occ_indices_spin (list): Indices of occupied spin orbitals.
        n_qubits: Number of qubits.

    Returns:
        hf_ref_state (scipy.sparse.csc.csc_matrix): The sparse matrix (vector)
            representation of HF wavefunction.
    """
    creation_op_list = [(i, 1) for i in occ_indices_spin]
    hf_creation_fermOp = openfermion.FermionOperator( tuple(creation_op_list), 1)
    hf_creation_spMat = openfermion.get_sparse_operator(hf_creation_fermOp, n_qubits)
    zero_state = openfermion.jw_configuration_state( [], n_qubits).reshape([-1, 1])
    zero_state_spMat = scipy.sparse.csc_matrix( (zero_state[zero_state.nonzero()], (zero_state.nonzero()[0], zero_state.nonzero()[1])), shape=zero_state.shape)
    hf_ref_state = hf_creation_spMat.dot(zero_state_spMat)
    return hf_ref_state


def _qubitOpList_to_spMatList_mp_worker(args: tuple):
    start_idx = args[0]
    end_idx = args[1]
    qubitOp_list = args[2]
    qubitOp_list_worker = qubitOp_list[start_idx:end_idx]
    n_qubits = args[3]
    n_params_worker = len(qubitOp_list_worker)
    spMat_list_worker = []
    for i in range(n_params_worker):
        spMat_list_worker.append(openfermion.get_sparse_operator(
            qubitOp_list_worker[i], n_qubits=n_qubits))
    return spMat_list_worker


def qubitOpList_to_spMatList_mp(
        qubitOp_list: list,
        n_qubits: int,
        n_procs: int = 1):
    """
    Convert openfermion's QubitOperators to sparce matrices in CSC format.

    Args:
        qubitOp_list (list): A list of openfermion's QubitOperator
        n_qubits (int): Number of qubits.
        n_procs (int): Number of processes to use.

    Returns:
        spMat_list (list): A list of corresponding sparse matrices.
    """
    spMat_list = []
    n_terms = len(qubitOp_list)
    n_workers = min(n_terms, n_procs)
    # if (n_workers != n_procs):
    #     print("Warning: change n_procs to %d" % (n_workers))

    chunk_size = n_terms // n_workers
    chunk_list = [chunk_size for i in range(n_workers)]
    for i in range(n_terms - chunk_size * n_workers):
        chunk_list[i] += 1

    import multiprocessing
    current_manager = multiprocessing.Manager()
    fermOp_list_shared = current_manager.list(qubitOp_list)
    args_workers = []
    start_idx = 0
    end_idx = 0
    for i in range(n_workers):
        start_idx = end_idx
        end_idx += chunk_list[i]
        args_workers.append((start_idx, end_idx, fermOp_list_shared, n_qubits))
    Pool = multiprocessing.Pool(n_workers)
    map_result = Pool.map(_qubitOpList_to_spMatList_mp_worker, args_workers)
    Pool.close()
    Pool.join()

    for i in range(n_workers):
        spMat_list += map_result[i]

    return spMat_list


def expm_multiply_quspin(A, vec):
    """
    Calculate exp(A) dot vec efficiently using QuSpin.
    If QuSpin is not available, the conventional
        scipy.sparse.linalg.expm_multiply() will be used instead.

    Args:
        A (scipy.sparse.csc.csc_matrix): A sparse matrix in CSC format.
        vec (scipy.sparse.csc.csc_matrix): A vector.

    Returns:
        expmv_result (numpy.ndarray): exp(A) dot vec

    Examples:
        >>> import numpy
        >>> import scipy.sparse
        >>> from sparse_matrices_helper import expm_multiply_quspin
        >>> m=1024
        >>> n=1024
        >>> A = scipy.sparse.random(m, n, density=0.001, format="csc")
        >>> vec = scipy.sparse.random(n, 1, density=0.001, format="csc")
        >>> result = expm_multiply_quspin(A, vec)
        >>> print(result)
        (128, 0)	(0.6830400798324581+0j)
        (708, 0)	(0.22998394355784096+0j)
        >>> from scipy.sparse.linalg import expm
        >>> result_2 = expm(A).dot(vec)
        >>> print(result_2)
        (708, 0)	0.22998394355784096
        (128, 0)	0.6830400798324581

    Notes:
        Make sure the OpenMP parallelism of this function
        DO NOT conflict with other MPI-based or pipe-based parallelism.
    """
    global expm_routine_msg
    try:
        import quspin.tools.evolution
    except ModuleNotFoundError:
        return scipy.sparse.linalg.expm_multiply(A, vec)
    try:
        import cuda_expm_multiply
        import cupyx.scipy.sparse
        import cupy
        A_gpu = cupyx.scipy.sparse.csc_matrix(A)
        vec_gpu = None
        if type(vec) is numpy.ndarray:
            vec_gpu = cupy.asarray(vec)
        elif type(vec) is scipy.sparse.csc_matrix:
            vec_gpu = cupyx.scipy.sparse.csc_matrix(vec)
        else:
            raise TypeError("Not supported type {}.".format(type(vec)))
        res = cuda_expm_multiply.expm_multiply(A_gpu, vec_gpu)
        if type(res) is cupy.ndarray:
            res = res.asnumpy()
        elif type(res) is cupyx.scipy.sparse.csc_matrix:
            res = res.get()
        return res
    except (ImportError, ModuleNotFoundError):
        pass
    dtype_real = numpy.float64
    dtype_cmplx = numpy.result_type(dtype_real, numpy.complex64)
    psi = vec.toarray().flatten().copy().astype(numpy.complex128)
    work_array = numpy.zeros((2 * psi.shape[0],), dtype=psi.dtype)
    eA = quspin.tools.evolution.expm_multiply_parallel(
        A, a=1.0, dtype=dtype_cmplx)
    eA.dot(psi, work_array=work_array, overwrite_v=True)
    expmv_result = scipy.sparse.csc_matrix(psi.reshape((psi.shape[0], 1)))
    return expmv_result


def tf(params, pauli):
    I=scipy.sparse.identity(pauli.shape[0])
    node1=math.sin(params)
    node2=math.cos(params)
    term=node2*I+node1*pauli
    return term


def A_psi(amplitues, clusterOp_spMat_list, hf_ref_state):
    n_params=len(amplitues)
    A=hf_ref_state
    for i in range(n_params):
        tf_term=tf(amplitues[i], clusterOp_spMat_list[i])
        A=scipy.sparse.csc_matrix.dot(tf_term, A)
    return A


def single_trotter_psi_spMat(
        amplitues: numpy.ndarray,
        clusterOp_spMat_list: list,
        hf_ref_state: scipy.sparse.csc.csc_matrix):
    """
    Calculate the single-Trotter UCC wavefunction:
        e^(i tn An) ... e^(i t3 A3) e^(i t2 A2) e^(i t1 A1)|psi>

    Args:
        amplitudes (numpy.ndarray): Variational parameters. Namely,
            coefficients of the UCC cluster operators.
        clusterOp_spMat_list (list): A list of sparse matrices in CSC format.
            The matrices are converted from UCC cluster operators.
        hf_ref_state (scipy.sparse.csc.csc_matrix): Sparse matrix (vector) of
            the Hartree-Fock reference state.

    Returns:
        Psi_spMat (scipy.sparse.csc.csc_matrix): Sparse matri (vector) of the
            single-Trotter UCC wavefunction.

    Notes:
        In this function, the imaginary "i" is supposed to be
        already in each matrix in clusterOp_spMat_list.
    """
    Psi_spMat = hf_ref_state
    n_params = len(amplitues)
    for i in range(n_params):
        Psi_spMat = expm_multiply_quspin(
            amplitues[i] * clusterOp_spMat_list[i], Psi_spMat)
    return Psi_spMat


def expm_energy_objective(
        amplitues: numpy.ndarray,
        clusterOp_spMat_list: list,
        hf_ref_state: scipy.sparse.csc.csc_matrix,
        hamiltonian_spMat: scipy.sparse.csc.csc_matrix,
        with_grad: bool = True,
        clusterOp_spMat_list_inv: list = None):
    """
    Evaluates energy of <psi|H|psi> under single-step Trotterization:
        |psi>=exp(cn tn) ... exp(c1 t1) exp(c1 t1)|HF>

    Args:
        amplitudes (numpy.ndarray): Variational parameters. Namely,
            coefficients of the UCC cluster operators.
        clusterOp_spMat_list (list): A list of sparse matrices in CSC format.
            The matrices are converted from UCC cluster operators.
        hf_ref_state (scipy.sparse.csc.csc_matrix): Sparse matrix (vector) of
            the Hartree-Fock reference state.
        hamiltonian_spMat (scipy.sparse.csc.csc_matrix): Sparse matrix of the
            system hamiltonian in CSC format.
        with_grad (bool): Whether to return the analytic gradients of energy
            w.r.t to variational parameters.
        clusterOp_spMat_list_inv (list): A list of sparse matrices in CSC format.
            The matrices are converted from UCC cluster operators * -1. If is not
            None, <psi| = <HF| expm(t0 inv0) expm(t1 inv1)... will be used instead
            of |psi>'s hermitian conjugated.

    Returns:
        energy (numpy.ndarray): Energy, E = <psi|H|psi>
        grad (numpy.ndarray): Gradients of energy w.r.t to
            variational parameters.

    Notes:
        The order of the Trotterized sequence from right to left
        is the same as clusterOp_spMat_list.
    """
    Psi_spMat = single_trotter_psi_spMat(
        amplitues, clusterOp_spMat_list, hf_ref_state)
    Psi_spMat_inv = None
    if clusterOp_spMat_list_inv is None:
        Psi_spMat_inv = Psi_spMat.T.conj()
    else:
        clusterOp_spMat_list_inv_dag = [
            i.T.conj() for i in clusterOp_spMat_list_inv]
        Psi_spMat_inv = single_trotter_psi_spMat(
            amplitues, clusterOp_spMat_list_inv_dag, hf_ref_state)
        Psi_spMat_inv = Psi_spMat_inv.T.conj()
    node1 = scipy.sparse.csc_matrix.dot(hamiltonian_spMat, Psi_spMat)
    node2 = scipy.sparse.csc_matrix.dot(Psi_spMat_inv, node1)
    normalization_factor = Psi_spMat_inv.dot(Psi_spMat)[0, 0]
    energy = node2[0, 0] / normalization_factor
    assert(numpy.isclose(numpy.imag(energy), 0.))
    if (with_grad is False):
        return energy.real

    if not numpy.allclose(
            Psi_spMat_inv.todense(), Psi_spMat.T.conj().todense()):
        raise ValueError("Grad is only implemented for Trottered UCC!")
    left_i = scipy.sparse.csc_matrix.dot(hamiltonian_spMat, Psi_spMat)
    right_i = Psi_spMat
    grad = numpy.zeros(len(amplitues), amplitues.dtype)
    for i in range(len(amplitues), 0, -1):
        grad_i = 2 * (left_i.T.conj().dot(
            clusterOp_spMat_list[i - 1].dot(right_i))).real
        grad[i - 1] = grad_i[0, 0]
        left_i = expm_multiply_quspin(
            -amplitues[i - 1] * clusterOp_spMat_list[i - 1], left_i)
        right_i = expm_multiply_quspin(
            -amplitues[i - 1] * clusterOp_spMat_list[i - 1], right_i)
    assert(numpy.isclose(numpy.linalg.norm(grad.imag), 0.0))

    return energy.real, grad


def tf_energy_objective(amplitues, clusterOp_spMat_list, hf_ref_state, hamiltonian_spMat, with_grad: bool):
    n_params=len(amplitues)
    psi=A_psi(amplitues, clusterOp_spMat_list, hf_ref_state)
    n1=scipy.sparse.csc_matrix.dot(hamiltonian_spMat, psi)
    n2=scipy.sparse.csc_matrix.dot(psi.T.conj(), n1)
    energy=n2[0,0].real
    if with_grad==False:
        return energy
    else:
        fun_array=numpy.zeros(n_params)
        left_i = scipy.sparse.csc_matrix.dot(hamiltonian_spMat, psi)
        right_i = psi
        for i in range(n_params, 0, -1):
            node1 =scipy.sparse.csc_matrix.dot(clusterOp_spMat_list[i-1], right_i)
            node2 =scipy.sparse.csc_matrix.dot(left_i.T.conj(), node1)
            fun_array[i-1] = node2[0,0].real*2
            tf_term=tf(amplitues[i-1], clusterOp_spMat_list[i-1])
            left_i = scipy.sparse.csc_matrix.dot(tf_term.T.conj(), left_i)
            right_i = scipy.sparse.csc_matrix.dot(tf_term.T.conj(), right_i)
        
        return energy, fun_array

