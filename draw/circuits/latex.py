import subprocess
import os

import fitz  # PyMuPDF

def convert_pdf_to_svg(pdf_path, svg_path):
    # 打开PDF文件
    doc = fitz.open(pdf_path)
    # 遍历每一页
    for page_num in range(len(doc)):
        page = doc.load_page(page_num)  # 加载页面
        svg = page.get_svg_image()  # 获取SVG格式的页面
        # 保存SVG文件
        with open(f"{svg_path}_{page_num + 1}.svg", "w") as svg_file:
            svg_file.write(svg)

def latex_file_to_svg(latex_file_path, output_svg_path):
    # 获取LaTeX文件的目录和文件名
    tex_dir, tex_file_name = os.path.split(latex_file_path)
    base_name = os.path.splitext(tex_file_name)[0]
    
    # 编译LaTeX文件为PDF
    subprocess.run(["pdflatex", latex_file_path], check=True, cwd=tex_dir)
    
    # # 获取生成的PDF文件路径
    pdf_file_path = os.path.join(tex_dir, f"{base_name}.pdf")

    convert_pdf_to_svg(pdf_file_path, output_svg_path)
    
    # # 将PDF文件转换为SVG
    # subprocess.run(["pdf2svg", pdf_file_path, output_svg_path], check=True)

# 示例用法

def draw_svg(latex_file_path, output_svg_path):
    f = open(latex_file_path, 'r')
    code = f.read()
    code = code.replace('{{}}_', '').replace(':', '')
    f = open(latex_file_path, 'w+')
    f.write(code)
    f.close()
    latex_file_to_svg(latex_file_path, output_svg_path)

class gate:
    def __init__(self, op, q):
        self.op = op
        self.q = q

from qiskit import QuantumCircuit
op_blocks = [[gate('cx', [5, 1]), gate('cx', [4, 0]), gate('cx', [1, 0])],
             [gate('cx', [3, 5]), gate('cx', [2, 4]), gate('cx', [5, 4])],
             [gate('cx', [5, 3]), gate('cx', [4, 2]), gate('cx', [3, 2])],
             [gate('swap', [2, 4]), gate('swap', [3, 5]), gate('cx', [3, 1]), gate('cx', [2, 0]), gate('cx', [1, 0])]
             ]
i = 0
for ops in op_blocks:
    qc = QuantumCircuit(6)
    for op in ops:
        if op.op == 'cx':
            qc.cx(op.q[0], op.q[1])
        if op.op == 'swap':
            qc.swap(op.q[0], op.q[1])
    qc.draw('latex_source', filename='./c_{}.tex'.format(i), idle_wires=False)
    draw_svg('./c_{}.tex'.format(i), './c_{}'.format(i))
    i += 1

exit()

# 将LaTeX文件编译为SVG
for i in range(5):
    latex_file_path = "./cir_{}.tex".format(i)  # 你的LaTeX文件路径
    output_svg_path = "./cir_{}".format(i)  # 你想要输出的SVG文件路径
    draw_svg(latex_file_path, output_svg_path)
