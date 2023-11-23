from compiler import PH_compiler, My_compiler
from arch import load_graph
import pickle
from mypauli import pauliString
from functions import print_qc

ac = load_graph('all12')
with open('./benchmark/LiH.pickle', 'rb') as f:
    block_list = pickle.load(f)
blocks = []
for b in block_list:
    blocks.append([pauliString(i[0], coeff=(i[1].real + i[1].imag)) for i in b])

print('n_block: ', len(blocks))
cp1 = My_compiler(blocks, ac)
cp2 = PH_compiler(blocks, ac)
res1, qc1 = cp1.compile()
res2, qc2 = cp2.compile()

print('our opt_level3: ', end="")
qc1 = print_qc(qc1, opt_level = 3)
print('ph opt_level3: ', end="")
qc2 = print_qc(qc2, opt_level = 3)
print('n_parameter: ', qc1.num_parameters)
