from openfermion import*
from function import*
from generator import *
from sparse_matrices_helper import *
import time

print('Generate time:',time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))

class mole_ansatz(object):

    def __init__(self, sys_information, cal_spMat=True, n_procs=8):
        self.nelectron=sys_information[0]
        self.hamiltonian_fermOp=sys_information[1]
        self.n_qubits=sys_information[2]
        self.occ_indices=sys_information[3]
        self.n_procs=n_procs
        self.cal_spMat=cal_spMat

    def UCCSD(self, n_orb_occ):
        n_orb=int(self.n_qubits/2)
        n_orb_vir = n_orb - n_orb_occ
        occ_indices = [i for i in range(n_orb_occ)]
        vir_indices = [i + n_orb_occ for i in range(n_orb_vir)]

        T1_singles = []
        T2_doubles = []
        for p_idx in range(len(vir_indices)):
            p = vir_indices[p_idx]
            for q_idx in range(len(occ_indices)):
                q = occ_indices[q_idx]
                tpq_list = SFEx_generator(p, q)
                for idx in range(len(tpq_list)):
                    tpq = tpq_list[idx]
                    tpq = tpq - openfermion.hermitian_conjugated(tpq)
                    tpq = openfermion.normal_ordered(tpq)
                    if (tpq.many_body_order() > 0):
                        T1_singles.append(tpq)

        for p_idx in range(len(vir_indices)):
            p = vir_indices[p_idx]
            for q_idx in range(p_idx, len(vir_indices)):
                q = vir_indices[q_idx]
                for r_idx in range(len(occ_indices)):
                    r = occ_indices[r_idx]
                    for s_idx in range(r_idx, len(occ_indices)):
                        s = occ_indices[s_idx]

                        tpqrs_list = DFEx_generator(p,q,r,s)
                        for idx in range(len(tpqrs_list)):
                            tpqrs = tpqrs_list[idx]
                            tpqrs = tpqrs -  openfermion.hermitian_conjugated(tpqrs)
                            tpqrs = openfermion.normal_ordered(tpqrs)
                            if (tpqrs.many_body_order() > 0):
                                T2_doubles.append(tpqrs)

        Op_pool=[]
        uccsd_fermOp = T1_singles + T2_doubles
        Op_pool.append(uccsd_fermOp)
        uccsd_QubitOp = [jordan_wigner(op) for op in uccsd_fermOp]
        Op_pool.append(uccsd_QubitOp)
        if self.cal_spMat==True:
            uccsd_spMatOp=qubitOpList_to_spMatList_mp(uccsd_QubitOp, self.n_qubits, n_procs=self.n_procs)
            Op_pool.append(uccsd_spMatOp)

        print("UCCSD pool size", str(len(uccsd_QubitOp)))

        return Op_pool

    def UCCGSD(self):
        n_orb=int(self.n_qubits/2)
        occ_indices = [i for i in range(n_orb)]
        vir_indices = [i for i in range(n_orb)]

        T1_singles = []
        T2_doubles = []

        for p_idx in range(len(vir_indices)):
            p = vir_indices[p_idx]
            for q_idx in range(len(occ_indices)):
                q = occ_indices[q_idx]
                tpq_list = SFEx_generator(p, q)
                for idx in range(len(tpq_list)):
                    tpq = tpq_list[idx]
                    tpq = tpq - hermitian_conjugated(tpq)
                    tpq = normal_ordered(tpq)
                    if (tpq.many_body_order() > 0):
                        T1_singles.append(tpq)

        pq = -1
        for p_idx in range(len(vir_indices)):
            p = vir_indices[p_idx]
            for q_idx in range(p_idx, len(vir_indices)):
                q = vir_indices[q_idx]
                pq += 1
                rs = -1
                for r_idx in range(len(occ_indices)):
                    r = occ_indices[r_idx]
                    for s_idx in range(r_idx, len(occ_indices)):
                        s = occ_indices[s_idx]
                        rs += 1
                        if (pq > rs):
                            continue
                        tpqrs_list = DFEx_generator(p,q,r,s)
                        for idx in range(len(tpqrs_list)):
                            tpqrs = tpqrs_list[idx]
                            tpqrs = tpqrs -hermitian_conjugated(tpqrs)
                            tpqrs = normal_ordered(tpqrs)
                            if (tpqrs.many_body_order() > 0):
                                T2_doubles.append(tpqrs)
        Op_pool=[]
        uccgsd_fermOp = T1_singles + T2_doubles
        Op_pool.append(uccgsd_fermOp)
        uccgsd_QubitOp = [jordan_wigner(op) for op in uccgsd_fermOp]
        Op_pool.append(uccgsd_QubitOp)
        if self.cal_spMat==True:
            uccsd_spMatOp=qubitOpList_to_spMatList_mp(uccgsd_QubitOp, self.n_qubits, n_procs=self.n_procs)
            Op_pool.append(uccsd_spMatOp)

        print("UCCGSD pool size", str(len(uccgsd_QubitOp)))

        return Op_pool
    
    def SA_UCCGSD(self):
        n_orb=int(self.n_qubits/2)
        occ_indices = [i for i in range(n_orb)]
        vir_indices = [i for i in range(n_orb)]

        T1_singles = []
        T2_doubles = []

        for p_idx in range(len(vir_indices)):
            p = vir_indices[p_idx]
            for q_idx in range(len(occ_indices)):
                q = occ_indices[q_idx]
                tpq_list = SA_SFEx_generator(p, q)
                for idx in range(len(tpq_list)):
                    tpq = tpq_list[idx]
                    tpq = tpq - hermitian_conjugated(tpq)
                    tpq = normal_ordered(tpq)
                    if (tpq.many_body_order() > 0):
                        T1_singles.append(tpq)

        pq = -1
        for p_idx in range(len(vir_indices)):
            p = vir_indices[p_idx]
            for q_idx in range(p_idx, len(vir_indices)):
                q = vir_indices[q_idx]
                pq += 1
                rs = -1
                for r_idx in range(len(occ_indices)):
                    r = occ_indices[r_idx]
                    for s_idx in range(r_idx, len(occ_indices)):
                        s = occ_indices[s_idx]
                        rs += 1
                        if (pq > rs):
                            continue
                        tpqrs_list = SA_DFEx_generator([p, q], [r, s])
                        for idx in range(len(tpqrs_list)):
                            tpqrs = tpqrs_list[idx]
                            tpqrs = tpqrs -hermitian_conjugated(tpqrs)
                            tpqrs = normal_ordered(tpqrs)
                            if (tpqrs.many_body_order() > 0):
                                T2_doubles.append(tpqrs)
        
        Op_pool=[]
        uccgsd_fermOp = T1_singles + T2_doubles
        Op_pool.append(uccgsd_fermOp)
        uccgsd_QubitOp = [jordan_wigner(op) for op in uccgsd_fermOp]
        Op_pool.append(uccgsd_QubitOp)
        if self.cal_spMat==True:
            uccsd_spMatOp=qubitOpList_to_spMatList_mp(uccgsd_QubitOp, self.n_qubits, n_procs=self.n_procs)
            Op_pool.append(uccsd_spMatOp)

        print("SA-UCCGSD pool size", str(len(uccgsd_QubitOp)))

        return Op_pool

    def QEB(self):
        orb_idx_list=get_all_spin_conservation_orbitals(self.n_qubits)
        qubit_Op_pool=[]
        for i in orb_idx_list:
            if len(i)==2:
                (p,q)=i
                qubit_Op_pool.append(SQEx_generator(p,q))
            if len(i)==4:
                (p,q,r,s)=i
                if check_spin_conv((p,q,r,s)) ==True:
                    qubit_Op_pool.append(DQEx_generator(p,q,r,s))
                if check_spin_conv((p,r,q,s)) ==True:
                    qubit_Op_pool.append(DQEx_generator(p,r,q,s))
                if check_spin_conv((p,s,q,r)) ==True:
                    qubit_Op_pool.append(DQEx_generator(p,s,q,r))
        op_pool=[]
        op_pool.append([])
        op_pool.append(qubit_Op_pool)
        if self.cal_spMat==True:
            spMatOp=qubitOpList_to_spMatList_mp(qubit_Op_pool, self.n_qubits, n_procs=self.n_procs)
            op_pool.append(spMatOp)
        print("QEB pool size",str(len(qubit_Op_pool)))
        return op_pool

    def RPs(self):
        orb_idx_list=get_all_spin_conservation_orbitals(self.n_qubits)
        qubit_Op_pool=[]
        for i in orb_idx_list:
            if len(i)==2:
                (p,q)=i
                qubit_Op_pool+=SPEx_generator(p,q)
            if len(i)==4:
                (p,q,r,s)=i
                qubit_Op_pool+=DPEx_generator(p,q,r,s)
        op_pool=[]
        op_pool.append([])
        op_pool.append(qubit_Op_pool)
        if self.cal_spMat==True:
            spMatOp=qubitOpList_to_spMatList_mp(qubit_Op_pool, self.n_qubits, n_procs=self.n_procs)
            op_pool.append(spMatOp)
        print("RPs pool size",str(len(qubit_Op_pool)))
        return op_pool


