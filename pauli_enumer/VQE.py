import datetime
import numpy
import math
from scipy.sparse import csc_matrix
from sparse_matrices_helper import *
from function import *
from openfermion import *
from generator import*

class vqe:

    def __init__(self, op_pool, hf_ref_state, sys_information, method='expm', adapt_tol:float=1e-3, adapt_maxiter:int=9999, n_procs=None):
            self.op_pool=op_pool
            self.hf_ref_state=hf_ref_state
            self.hamiltonian_fermOp=sys_information[1]
            self.n_qubits=sys_information[2]
            self.method=method
            self.adapt_tol=adapt_tol
            self.adapt_maxiter=adapt_maxiter
            self.n_procs=n_procs


    def adapt_vqe(self, add_para_num:int=1):
        psi_adapt_array = self.hf_ref_state
        new_operator_pool = []                  
        new_parameters = []                    
        new_operator_pool_idx = []              
        iter_idx = 0
        convergence=False
        err=[]

        hamiltonian_spMat=get_sparse_operator(self.hamiltonian_fermOp, self.n_qubits)
        e_FCI,v_FCI = check_FCI(hamiltonian_spMat, self.hf_ref_state.todense())
        print("FCI energy:"," %.8f"%e_FCI)

        table1="{0:^5}{1:^10}{2:^12}{3:^12}{4:^10}{5:^10}{6:^12}"
        print(table1.format("Iter", "ΔH", "Ei", "Error", "s2", "sz", "opt-time"))

        while convergence==False and (iter_idx <= self.adapt_maxiter): 
            iter_idx += 1

            operator_list_idx=grad_select(self.op_pool[2], psi_adapt_array, hamiltonian_spMat,add_para_num=add_para_num, n_procs=self.n_procs)

            for i in operator_list_idx:
                new_operator_pool.append(self.op_pool[2][i])                  
                new_parameters.append(0.)
                new_operator_pool_idx.append(i)
            new_parameters_array = numpy.array(new_parameters)

            start3=datetime.datetime.now()
            energy_i, new_parameters, psi_adapt_array=VQE_para_optimize(self.method, new_parameters_array, new_operator_pool, self.hf_ref_state, hamiltonian_spMat)
            end3=datetime.datetime.now()
            
            E0=energy_i
            Hv=csc_matrix.dot(hamiltonian_spMat, psi_adapt_array)
            H2=csc_matrix.dot(Hv.T.conj(), Hv)[0,0].real
            E2=(csc_matrix.dot(psi_adapt_array.T.conj(), Hv)[0,0].real)**2
            delta_H=numpy.sqrt(abs(H2-E2))
            convergence=check_convergence(delta_H, self.adapt_tol)

            error=abs(e_FCI-energy_i)
            err.append(float("%.4e"%error))
            [N_spMat, N_spMat_list, s2_spMat, sz_spMat]=phy_operator(n_qubits=self.n_qubits)[1]
            s2=measure_expection(s2_spMat, psi_adapt_array)
            sz=measure_expection(sz_spMat, psi_adapt_array)
            
            print(table1.format(iter_idx, 
                                "%.4f"%delta_H, 
                                "%.6f"%energy_i, 
                                "%.2e"%error,
                                "%.2e"%s2,
                                "%.2e"%sz, 
                                "%.2e"%(end3-start3).total_seconds()))
            
        end_warning('Convergence energy(Ha)', 'Convergence error(exp)', "%.6f"%energy_i, "%.4e"%error)

        return  energy_i, psi_adapt_array, new_operator_pool_idx, new_parameters, error


