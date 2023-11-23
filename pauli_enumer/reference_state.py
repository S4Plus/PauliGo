import openfermion
import scipy
class single_reference_state(object):

    def __init__(self, sys_information):
        self.n_qubits=sys_information[2]
        self.occ_indices=sys_information[3]
        
    def HF_state(self):
        hf_ref_states=openfermion.jw_configuration_state(self.occ_indices, self.n_qubits).reshape([-1, 1])
        hf_ref_states=scipy.sparse.csc_matrix(hf_ref_states)
        return hf_ref_states
