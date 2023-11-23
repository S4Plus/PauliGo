import numpy
import scipy
import scipy.sparse
import scipy.optimize
import copy
import itertools
import pickle
from openfermion import FermionOperator, get_sparse_operator
from sparse_matrices_helper import A_psi, single_trotter_psi_spMat, tf_energy_objective, expm_energy_objective,expm_multiply_quspin


def check_odd_even_num(idx_tuple):
    odd=0
    even=0
    for i in idx_tuple:
        if i % 2==1:
            odd+=1
        if i % 2==0:
            even+=1
    return (odd, even)


def get_all_spin_conservation_orbitals(n_qubits):
    idx_list=[i for i in range(n_qubits)]
    combine_spin_idx= list(itertools.combinations(idx_list, 2))+list(itertools.combinations(idx_list,4))
    spin_conv_idx=[]
    for i in combine_spin_idx:
        if check_odd_even_num(i) in [(0,2), (2,0), (0,4), (4,0), (2,2)]:
            spin_conv_idx.append(i)
    return spin_conv_idx


def two_spa_combine(p,q):
    """
    p < q is vir ; r < s is occ, r < p < s < q
    """
    orb_idx_list=[]
    pa = p * 2 + 0
    pb = p * 2 + 1
    qa = q * 2 + 0
    qb = q * 2 + 1
    orb_idx_list.extend([[pa,qa],[pb,qb]])
    return orb_idx_list


def four_spa_combine(p,q,r,s):
    """
    p < q is vir ; r < s is occ, r < p < s < q
    """
    orb_idx_list=[]
    pa = p * 2 + 0
    pb = p * 2 + 1
    qa = q * 2 + 0
    qb = q * 2 + 1
    ra = r * 2 + 0
    rb = r * 2 + 1
    sa = s * 2 + 0
    sb = s * 2 + 1
    orb_idx_list.extend([[pa,qa,ra,sa],[pb,qb,rb,sb],[pb,qb,ra,sa],[pb,qa,rb,sa],[pa,qb,ra,sb],[pa,qa,rb,sb]])
    return orb_idx_list


def phy_operator(n_qubits):
    """
    Get N and spin S2 and Sz fermion operator
    
    Args:
        norb: number of activate orbital.

    Return:
        N:FermionOperator.
        n_list:FermionOperator.
        s2:FermionOperator.
        sz:FermionOperator.
    """
    norb=int(n_qubits/2)
    N_fermion = FermionOperator()
    N_fermion_list=[]
    N_spMat_list=[]
    for i in range(n_qubits):
        Ni_fermion=FermionOperator(((i,1),(i,0)),1.0)
        Ni_spMat=get_sparse_operator(Ni_fermion, n_qubits)

        N_fermion += Ni_fermion
        N_fermion_list.append(Ni_fermion)
        N_spMat_list.append(Ni_spMat)

    N_spMat=get_sparse_operator(N_fermion, n_qubits)

    Sx = FermionOperator()
    for i in range(norb):
        Sx += FermionOperator(((2*i,1),(2*i+1,0)),0.5) + FermionOperator(((2*i+1,1),(2*i,0)),0.5)

    Sy = FermionOperator()
    for i in range(norb):
        Sy += FermionOperator(((2*i,1),(2*i+1,0)),-0.5*1j) + FermionOperator(((2*i+1,1),(2*i,0)),0.5*1j)

    Sz = FermionOperator()
    for i in range(norb):
        Sz += FermionOperator(((2*i,1),(2*i,0)),0.5) - FermionOperator(((2*i+1,1),(2*i+1,0)),0.5)

    S2 = Sx*Sx + Sy*Sy + Sz*Sz
    S2_spMat=get_sparse_operator(S2, n_qubits)
    Sz_spMat=get_sparse_operator(Sz, n_qubits)

    phy_fermion_list=[N_fermion, N_fermion_list, S2, Sz]
    phy_spMat_list=[N_spMat, N_spMat_list, S2_spMat, Sz_spMat]

    return phy_fermion_list, phy_spMat_list


def select_Nst_max_element(array_or_list, N):
    t= copy.deepcopy(array_or_list)
    max_element=[]
    max_idx=[]
    for i in range(N):
        select_idx= numpy.argmax(t)
        select_element=array_or_list[select_idx]
        max_idx.append(select_idx)
        max_element.append(select_element)
        t[select_idx]=0
    return max_idx, max_element


def check_FCI(hamiltonian_spMat, ref_states):
    e_FCI, v_FCI = scipy.sparse.linalg.eigsh(hamiltonian_spMat, k=1, which="SA", v0=ref_states)
    return e_FCI, v_FCI


def end_warning(strings1, strings2, val1, val2):
    print("\n")
    head="{0:<26}{1:^10}"
    table="{0:<1}{1:<30}{2:^1}{3:<30}{4:>1}"
    print(head.format(" ",  "Information"))
    print("+-------------------------------------------------------------+")
    print(table.format("|", strings1,       "|",    val1,               "|"))
    print(table.format("|", strings2,       "|",    val2,               "|"))
    print("+-------------------------------------------------------------+")
    print('\n')
    print("===============================================================")
    print("=============================END===============================")
    print("===============================================================")
    print('\n')
    return


def measure_expection(operator, psi_arrary):
    if operator ==None:
        expection=scipy.sparse.csc_matrix.dot(psi_arrary.T.conj(), psi_arrary)[0,0].real
    else:
        node= scipy.sparse.csc_matrix.dot(operator, psi_arrary)
        expection=scipy.sparse.csc_matrix.dot(psi_arrary.T.conj(), node)[0,0].real
    return expection


def minimize(fun, x0: numpy.ndarray, args: tuple = (),
             method: str = "Adam", options: dict = {},
             jac: bool = True,
             callback=None):
    """
    A wrapper for scipy and pyTorch's optimizers.

    Args:
        fun (callable): A function which can evaluate fun(x).
        x0 (numpy.ndarray): Initial guess.
        args (tuple): Additional parameters for fun.
        method (str): Optimization methods.
        options (dict): Options for the optimizer.
        jac (bool): Whether the jacobian is calculated by fun.
        callback (callable): A call-back function. Not used
            for PyTorch's optimizers.

    Supported methods:
        "Adam",
        "Adadelta",
        "Adagrad",
        "AdamW",
        "Adamax",
        "ASGD",
        "NAdam",
        "RAdam",
        "RMSProp",
        "Rprop",
        "SGD"

    Settings in options:
        "lr" (float): default 0.1
        "momentum" (float): default 0.0, only used when method="SGD".
        "maxiter" (int): default 999.
        "ftol" (float): default 1e-8.
        "gtol" (float): default 1e-5.
        "disp" (bool): default False.
    """

    scipy_methods = [
        "Nelder-Mead",
        "Powell",
        "CG",
        "BFGS",
        "Newton-CG",
        "L-BFGS-B",
        "TNC",
        "COBYLA",
        "SLSQP",
        "trust-constr",
        "dogleg",
        "trust-ncg",
        "trust-exact",
        "trust-krylov"
    ]
    if method in scipy_methods:
        import scipy.optimize
        result_scipy = scipy.optimize.minimize(
            fun=fun,
            x0=x0,
            args=args,
            method=method,
            options=options,
            jac=jac,
            callback=callback
        )
        return result_scipy

    if jac is not True:
        raise ValueError("jac must be True.")

    import torch

    class result_class(object):
        def __init__(self, x: numpy.ndarray, fun: numpy.float64):
            self.x = x
            self.fun = fun

    x = torch.from_numpy(x0)
    x.requires_grad_(True)
    x.grad = torch.zeros(x.shape, dtype=x.dtype)

    lr = 0.1
    maxiter = 999
    ftol = 1e-8
    gtol = 1e-5
    disp = False
    if "lr" in options:
        lr = options["lr"]
    if "maxiter" in options:
        maxiter = options["maxiter"]
    if "ftol" in options:
        ftol = options["ftol"]
    if "gtol" in options:
        gtol = options["gtol"]
    if "disp" in options:
        disp = options["disp"]

    optimizer = None
    support_methods = [
        "Adam",
        "Adadelta",
        "Adagrad",
        "AdamW",
        "Adamax",
        "ASGD",
        "NAdam",
        "RAdam",
        "RMSProp",
        "Rprop",
        "SGD",
        "LBFGS"
    ]
    if method == "Adam":
        optimizer = torch.optim.Adam([x], lr=lr)
    elif method == "Adadelta":
        optimizer = torch.optim.Adadelta([x], lr=lr)
    elif method == "Adagrad":
        optimizer = torch.optim.Adagrad([x], lr=lr)
    elif method == "AdamW":
        optimizer = torch.optim.AdamW([x], lr=lr)
    elif method == "Adamax":
        optimizer = torch.optim.Adamax([x], lr=lr)
    elif method == "ASGD":
        optimizer = torch.optim.ASGD([x], lr=lr)
    elif method == "NAdam":
        optimizer = torch.optim.NAdam([x], lr=lr)
    elif method == "RAdam":
        optimizer = torch.optim.RAdam([x], lr=lr)
    elif method == "RMSProp":
        optimizer = torch.optim.RMSprop([x], lr=lr)
    elif method == "Rprop":
        optimizer = torch.optim.Rprop([x], lr=lr)
    elif method == "SGD":
        momentum = 0.
        if "momentum" in options:
            momentum = options["momentum"]
        optimizer = torch.optim.SGD([x], lr=lr, momentum=momentum)
    elif method == "LBFGS":
        optimizer = torch.optim.LBFGS(
            [x], lr=lr, line_search_fn="strong_wolfe")
    else:
        raise NotImplementedError("method must be one of these: {}\
".format(support_methods))

    iter_count = 0
    f_diff = ftol * 9999
    f_last = None
    f_val_dict = {}
    g_norm_dict = {}

    def _closure():  # ??
        global f_last
        optimizer.zero_grad()
        x_func = x.detach().numpy().copy()
        f, df = fun(x_func, *args)
        x.grad = torch.from_numpy(df)
        f_val_dict.update({hash(x_func.tobytes()): f})
        g_norm_dict.update(
            {hash(x_func.tobytes()): torch.linalg.norm(x.grad).item()})
        return f

    finish_type = 0
    while iter_count < maxiter:
        if iter_count > 0:
            if grad_last <= gtol:
                finish_type = 1
                break
            if f_diff <= ftol:
                finish_type = 2
                break
        x_func_last = x.detach().numpy().copy()
        optimizer.step(_closure)
        grad_last = g_norm_dict[hash(x_func_last.tobytes())]
        f_cur = f_val_dict[hash(x_func_last.tobytes())]
        if f_last is not None:
            f_diff = abs(f_cur - f_last)
        f_last = f_cur
        if disp:
            print("Iter %5d f=%20.16f |g|=%20.16f" %
                  (iter_count, f_cur, grad_last))
        iter_count += 1

    if disp:
        finish_reason = ""
        if finish_type == 0:
            finish_reason = "maxiter"
        elif finish_type == 1:
            finish_reason = "gtol"
        elif finish_type == 2:
            finish_reason = "ftol"
        print("Finished due to %s" % (finish_reason))

    result = result_class(x_func_last, f_last)

    return result


def grad_cal(clusterOp_list, psi_arrary, operator):
    n_params=len(clusterOp_list)
    pre_est_grad=numpy.zeros(n_params)
    pre_est_grad_abs=numpy.zeros(n_params)
    for i in range(n_params):
        ipsi_adapt_array    =scipy.sparse.csc_matrix.dot( clusterOp_list[i], psi_arrary)                                            
        Hipsi_adapt_array   =scipy.sparse.csc_matrix.dot( operator, ipsi_adapt_array)                            
        pre_est_grad_i      =scipy.sparse.csc_matrix.dot( psi_arrary.T.conj(), Hipsi_adapt_array)[0,0].real*-2
        pre_est_grad[i]=(pre_est_grad_i)
        pre_est_grad_abs[i] =numpy.abs(pre_est_grad_i) 
    return pre_est_grad, pre_est_grad_abs


def _grad_cal_mp_worker(args_worker: list):
    start_idx = args_worker[0]
    end_idx = args_worker[1]
    size_operator_pool_worker = end_idx - start_idx
    pre_est_grad_worker = numpy.zeros(size_operator_pool_worker)
    pre_est_grad_abs_worker= numpy.zeros(size_operator_pool_worker)
    for i in range(start_idx, end_idx):
        op_i = operator_pool_sparse_mat_global[i]
        ipsi_adapt_array = scipy.sparse.csc_matrix.dot(op_i,psi_adapt_global)
        Hipsi_adapt_array = scipy.sparse.csc_matrix.dot(operator_spMat_global,ipsi_adapt_array)
        pre_est_grad_i = scipy.sparse.csc_matrix.dot(psi_adapt_global.T.conj(), Hipsi_adapt_array).real[0,0]*2*-1
        pre_est_grad_worker[i - start_idx] = pre_est_grad_i
        pre_est_grad_abs_worker[i - start_idx] = numpy.abs(pre_est_grad_i)
    return pre_est_grad_worker, pre_est_grad_abs_worker


def grad_cal_mp(operator_pool_sparse_mat: list, psi_adapt: scipy.sparse.csc_matrix, operator_spMat: scipy.sparse.csc_matrix, n_procs: int = 1) -> numpy.ndarray:
    size_operator_pool = len(operator_pool_sparse_mat)
    n_terms = size_operator_pool
    n_workers = min(n_terms, n_procs)
    if (n_workers != n_procs):
        print("Warning: change n_procs to %d" % (n_workers))

    chunk_size = n_terms // n_workers
    chunk_list = [chunk_size for i in range(n_workers)]
    for i in range(n_terms - chunk_size * n_workers):
        chunk_list[i] += 1

    global operator_pool_sparse_mat_global
    global psi_adapt_global
    global operator_spMat_global
    operator_pool_sparse_mat_global = operator_pool_sparse_mat
    psi_adapt_global = psi_adapt
    operator_spMat_global = operator_spMat

    start_idx = 0
    end_idx = 0
    args_workers = []
    for i in range(n_workers):
        start_idx = end_idx
        end_idx += chunk_list[i]
        args_workers.append((start_idx, end_idx))

    import multiprocessing
    Pool = multiprocessing.Pool(n_workers)
    map_result = Pool.map(_grad_cal_mp_worker, args_workers)
    Pool.close()
    Pool.join()

    pre_est_grad = None
    pre_est_grad_abs = None
    for i in range(n_workers):
        if i == 0:
            pre_est_grad=map_result[i][0]
            pre_est_grad_abs = map_result[i][1]
        else:
            pre_est_grad = numpy.concatenate((pre_est_grad, map_result[i][0]))
            pre_est_grad_abs = numpy.concatenate((pre_est_grad_abs, map_result[i][1]))
    return pre_est_grad, pre_est_grad_abs


def grad_select(Operator_spMat_list, psi_spMat, hamiltonian_spMat, add_para_num, n_procs):
    pre_est_grad_abs=gradient_calculation(Operator_spMat_list, psi_spMat, hamiltonian_spMat,n_procs)[1]
    operator_list_idx=select_Nst_max_element(pre_est_grad_abs, add_para_num)[0]
    return operator_list_idx


def single_para_energy_objective(amplitues, clusterOp_spMat, psi, hamiltonian_spMat):
    psi=expm_multiply_quspin(amplitues[0]*clusterOp_spMat, psi)
    n1=scipy.sparse.csc_matrix.dot(hamiltonian_spMat, psi)
    n2=scipy.sparse.csc_matrix.dot(psi.T.conj(), n1)
    energy=n2[0,0].real
    left_i = scipy.sparse.csc_matrix.dot(hamiltonian_spMat, psi)
    right_i = psi
    node1 =scipy.sparse.csc_matrix.dot(clusterOp_spMat, right_i)
    node2 =scipy.sparse.csc_matrix.dot(left_i.T.conj(), node1)
    fun_array = node2[0,0].real*2
    fun_array= numpy.array([fun_array])
    return energy, fun_array


def delta_E_calculator(clusterOp_spMat_list, psi, hamiltonian_spMat, E0):
    delta_E=[]
    for i in clusterOp_spMat_list:
        amplitues=[0.]
        amplitues_array=numpy.array(amplitues)
        res=minimize(single_para_energy_objective, amplitues_array, args=(i, psi, hamiltonian_spMat), jac=True, method="BFGS")
        delta_E.append(E0-res.fun)
    return delta_E


def _delta_E_calculator_mp_worker(args_worker: list):
    start_idx = args_worker[0]
    end_idx = args_worker[1]
    size_operator_pool_worker = end_idx - start_idx
    pre_est_delta_E_worker = numpy.zeros(size_operator_pool_worker)
    for i in range(start_idx, end_idx):
        op_i = operator_pool_sparse_mat_global[i]
        res=scipy.optimize.minimize(single_para_energy_objective, numpy.array([0.]), args=(op_i, psi_adapt_global, operator_spMat_global), jac=True, method="BFGS", options={'gtol':1e-8})
        pre_est_delta_E_worker[i - start_idx] = E0_global-res.fun
    return pre_est_delta_E_worker


def delta_E_calculator_mp(operator_pool_sparse_mat: list, psi_adapt: scipy.sparse.csc_matrix, operator_spMat: scipy.sparse.csc_matrix, E0:float, n_procs: int = 1) -> numpy.ndarray:
    size_operator_pool = len(operator_pool_sparse_mat)
    n_terms = size_operator_pool
    n_workers = min(n_terms, n_procs)
    if (n_workers != n_procs):
        print("Warning: change n_procs to %d" % (n_workers))

    chunk_size = n_terms // n_workers
    chunk_list = [chunk_size for i in range(n_workers)]
    for i in range(n_terms - chunk_size * n_workers):
        chunk_list[i] += 1

    global operator_pool_sparse_mat_global
    global psi_adapt_global
    global operator_spMat_global
    global E0_global
    operator_pool_sparse_mat_global = operator_pool_sparse_mat
    psi_adapt_global = psi_adapt
    operator_spMat_global = operator_spMat
    E0_global=E0

    start_idx = 0
    end_idx = 0
    args_workers = []
    for i in range(n_workers):
        start_idx = end_idx
        end_idx += chunk_list[i]
        args_workers.append((start_idx, end_idx))

    import multiprocessing
    Pool = multiprocessing.Pool(n_workers)
    map_result = Pool.map(_delta_E_calculator_mp_worker, args_workers)
    Pool.close()
    Pool.join()

    pre_est_delta_E = None
    for i in range(n_workers):
        if i == 0:
            pre_est_delta_E=map_result[i]
        else:
            pre_est_delta_E = numpy.concatenate((pre_est_delta_E, map_result[i]))
    return pre_est_delta_E


def deltaE_select(Operator_spMat_list, psi_spMat, hamiltonian_spMat, E0, add_para_num, n_procs):
    e_list=delta_E_calculator_mp(Operator_spMat_list, psi_spMat, hamiltonian_spMat, E0, n_procs)
    operator_list_idx=select_Nst_max_element(e_list, add_para_num)[0]
    return operator_list_idx

    
def gradient_calculation(Operator_spMat_list, psi_spMat, hamiltonian_spMat, n_procs):
    if n_procs ==None:
        pre_est_grad= grad_cal(Operator_spMat_list, psi_spMat, hamiltonian_spMat)
    else:
        pre_est_grad= grad_cal_mp(Operator_spMat_list, psi_spMat, hamiltonian_spMat, n_procs=n_procs)
    return pre_est_grad[0], pre_est_grad[1]


def check_convergence(delta_H, tol):
    convergence=False
    if delta_H < tol:
        convergence=True
    return convergence


def VQE_para_optimize(method, para_array, operator_list, hf_state, hamiltonian, grad=True):
    if method=='expm':
        res= minimize(expm_energy_objective,            
                    para_array,        
                    args=(operator_list, hf_state, hamiltonian, grad),
                    jac=True,
                    method="LBFGS",
                    options={'lr':1.0,'disp':False,'gtol':1e-5})
        optimize_para = list(res.x)
        optimize_array = single_trotter_psi_spMat(optimize_para, operator_list, hf_state)
        optimize_enery=res.fun

    if method=='tri':
        res = minimize(tf_energy_objective,            
                    para_array,        
                    args=(operator_list, hf_state, hamiltonian, grad),
                    jac=True,
                    method="LBFGS",
                    options={'lr':1.0,'disp':False,'gtol':1e-5})
        optimize_para = list(res.x)
        optimize_array = A_psi(optimize_para, operator_list, hf_state)
        optimize_enery=res.fun
    return optimize_enery, optimize_para, optimize_array


def VQE_hamiltonian_optimize(method, para_array, operator_list, hf_state, hamiltonian, grad=True):
    if method=='expm':
        res= minimize(expm_energy_objective,            
                    para_array,        
                    args=(operator_list, hf_state, hamiltonian, grad),
                    jac=True,
                    method="LBFGS",
                    options={'lr':1.0,'disp':False})
        optimize_para = list(res.x)
        optimize_array = single_trotter_psi_spMat(optimize_para, operator_list, hf_state)
        optimize_enery=res.fun

    if method=='tri':
        res = minimize(tf_energy_objective,            
                    para_array,        
                    args=(operator_list, hf_state, hamiltonian, grad),
                    jac=True,
                    method="LBFGS",
                    options={'lr':1.0,'disp':False})
        optimize_para = list(res.x)
        optimize_array = A_psi(optimize_para, operator_list, hf_state)
        optimize_enery=res.fun
    return optimize_enery, optimize_para, optimize_array


def save_variable(v,filename, method='wb'):
  f=open(filename, method)
  pickle.dump(v,f)
  f.close()
  return filename


def load_variavle(filename):
  f=open(filename,'rb')
  r=pickle.load(f)
  f.close()
  return r
    
    
def allclose_v2(A, B, atol = 1e-8):
    # If you want to check matrix shapes as well
    if numpy.array_equal(A.shape, B.shape)==0:
        return False

    r1,c1 = A.nonzero()
    r2,c2 = B.nonzero()

    lidx1 = numpy.ravel_multi_index((r1,c1), A.shape)
    lidx2 = numpy.ravel_multi_index((r2,c2), B.shape)

    sidx1 = lidx1.argsort()
    sidx2 = lidx2.argsort()

    index_match = numpy.array_equal(lidx1[sidx1], lidx2[sidx2])
    if index_match==0:
        return False
    else:  
        v1 = A.data
        v2 = B.data        
        V1 = v1[sidx1]
        V2 = v2[sidx2]        
        return numpy.allclose(V1,V2, atol=atol)
    

def check_01(t):
    n0=0
    n1=0
    for i in t:
        if i==0:
            n0+=1
        elif i==1:
            n1+=1
    return (n0,n1)


def check_spin_conv(t):
    (p,q,r,s)=t
    (i,j,k,l)=(p%2,q%2,r%2,s%2)
    if check_01((i,j)) == check_01((k,l)):
        return True
    else:
        return False


