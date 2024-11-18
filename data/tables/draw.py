# import os
# import pandas as pd

# # 指定文件夹路径
# folder_path = '.'

# # 遍历文件夹
# for root, dirs, files in os.walk(folder_path):
#     for file in files:
#         if file.endswith('.xlsx'):
#             # 获取文件的完整路径
#             file_path = os.path.join(root, file)
#             print(file_path)

#             # 读取Excel文件
#             df = pd.read_excel(file_path, sheet_name='graphResult')

#             # 遍历DataFrame的每一行
#             for index, row in df.iterrows():
#                 # 遍历行中的每个单元格
#                 for col_name, cell_value in row.iteritems():
#                     print(f"Row {index}, Column {col_name}, Value: {cell_value}")

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
# 指定文件夹路径
folder_path = Path('.')

# 存储每个工作表的数据
data = []

# 遍历文件夹中的所有.xlsx文件
for file_path in folder_path.glob('*.xlsx'):
    df = pd.read_excel(file_path, sheet_name='graphResult')
    data.append(df.iloc[:, [0, 1]])  # 只取第一列和第二列

# 设置横坐标（假设每个sheet的第一列是相同的）
x_labels = data[0].iloc[:, 0]
def to01(num):
    c = ''
    for i in range(4):
        c = str(num % 10) + c
        num //= 10
    return c
x_labels = [to01(num) for num in x_labels]
print(x_labels)
x = np.array(range(len(x_labels)))
width = 0.4

# 创建柱状图
plt.figure(figsize=(10, 6))

# 为每个工作表绘制柱状图
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']  # 定义颜色列表
labels = ['PauliForest', 'Paulihedral']
styles = [# {'edgecolor' : 'black', 'hatch' : '//', 'fill' : False}, 
#               {'edgecolor' : 'black', 'hatch' : '\\\\', 'fill' : False},
              {'edgecolor' : 'black', 'fill' : False}, 
              {'color' : 'black'}]
for i, dataset in enumerate(data):
    plt.bar(x - width / 2 + i * width, dataset.iloc[:, 1], width, label=labels[i], **styles[i])
    print(max(dataset.iloc[:, 1]))

# 添加图例
plt.legend(fontsize=14)

# 设置横坐标标签
plt.xticks(x, x_labels)

# 添加标题和轴标签
# plt.title('Distribution of measurement results')
# plt.xlabel('X Axis Label')
plt.ylabel('Proportion', fontsize=14)
plt.yticks(fontsize=14)

# 显示图形
plt.savefig('./prob.png')