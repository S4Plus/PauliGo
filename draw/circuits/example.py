from qiskit import QuantumCircuit
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

qc = QuantumCircuit(8)
for i in range(7):
    qc.cx(i, i+1)

qc.draw('latex_source', filename='./motiExam1.tex', idle_wires=False)
draw_svg('./motiExam1.tex', './motiExam1')

qc = QuantumCircuit(8)
for i in range(0, 8, 2):
    if i!=4:
        qc.cx(i, i+1)
    else:
        qc.cx(i+1, i)
# for i in range(1, 8, 4):
#     qc.cx(i, i+2)
# qc.cx(3, 7)
for i in [1,3,7]:
    qc.cx(i, 4)

qc.draw('latex_source', filename='./motiExam2.tex', idle_wires=False)
draw_svg('./motiExam2.tex', './motiExam2')