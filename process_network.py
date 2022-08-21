#encoding:utf-8
import numpy as np
import pickle
import pandas as pd
import time
import re

# with open('E:\迅雷下载\文件\9606.protein.links.detailed.v11.5.txt\\9606.protein.links.detailed.v11.5.txt.gz.inproc', 'r') as f:
#     while True:
#                 data = f.readline()
#                 if not data:
#                     break
#                 str_data = (data.split(' ')[0] + '\t' + data.split(' ')[1] + '\t' + data.split(' ')[5] + '\n')
#                 print(str_data)
#                 time.sleep(10)
#                 #fw.write(str_data)

with open('E:\\shiyan_data\\BRCA\\String\\9606.protein.gene_textmining.txt', 'r') as f:
    gene_combined = f.readlines()

with open('E:\\shiyan_data\\BRCA\\String\\9606.protein.gene_experiment.txt', 'r') as f:
    gene_experiment = f.readlines()

with open('E:\\shiyan_data\\GBM\\lab\\Common_gene_diff_exp.bin', 'rb') as f:
    common_gene = pickle.load(f)

select_gene_combined = []
select_gene_combined.append(gene_combined[0])
select_gene_experiment = []
select_gene_experiment.append(gene_combined[0])
# print(select_gene_combined)
#     # time.sleep(10)

start = time.time()
for i in range(1, len(gene_combined)):
    if gene_combined[i].split('\t')[0] in common_gene and gene_combined[i].split('\t')[1] in common_gene:
        select_gene_combined.append(gene_combined[i])
    if gene_experiment[i].split('\t')[0] in common_gene and gene_experiment[i].split('\t')[1] in common_gene:
        select_gene_experiment.append(gene_experiment[i])
end = time.time()
print("running time:", end - start)

with open('E:\\shiyan_data\\GBM\\lab\\String\\select.protein.gene.texting.txt', 'w') as f:
    for i in range(len(select_gene_combined)):
        f.write(select_gene_combined[i])

with open('E:\\shiyan_data\\GBM\\lab\\String\\select.protein.gene.experiment.txt', 'w') as f:
    for i in range(len(select_gene_experiment)):
        f.write(select_gene_experiment[i])


with open('E:\\shiyan_data\\GBM\\lab\\String\\select.protein.gene.texting.txt', 'r') as f:
    gene_texting = f.readlines()

with open('E:\\shiyan_data\\GBM\\lab\\String\\select.protein.gene.experiment.txt', 'r') as f:
    gene_experiment = f.readlines()

replace_gene = []
replace_gene.append(gene_experiment[0])
for i in range(1, len(gene_texting)):
    if (gene_texting[i].split('\t')[0] == gene_experiment[i].split('\t')[0]) and (gene_texting[i].split('\t')[1] == gene_experiment[i].split('\t')[1]):
        if int(gene_experiment[i].split('\t')[2].split('\n')[0]) != 0:
            gene_num = max(gene_texting[i].split('\t')[2].split('\n')[0], gene_experiment[i].split('\t')[2].split('\n')[0])
        else:
            gene_num = 0
        replace_gene.append(gene_texting[i].split('\t')[0] + '\t' + gene_texting[i].split('\t')[1] + '\t' + str(gene_num) + '\n')


with open('E:\\shiyan_data\\GBM\\lab\\String\\select.protein.gene.all.txt', 'w') as f:
    for i in range(len(replace_gene)):
        f.write(replace_gene[i])


with open('E:\\shiyan_data\\GBM\\lab\\Common_gene_diff_exp.bin', 'rb') as f:
    common_gene_2 = pickle.load(f)

with open('E:\\shiyan_data\\GBM\\lab\\Common_gene_diff_exp.bin', 'rb') as f:
    common_gene = pickle.load(f)

f1 = open('E:\\shiyan_data\\GBM\\lab\\String\\select.protein.gene.all.txt', 'r')
select_gene_combined = f1.read()
f1.close()

n = len(common_gene)
print(len(common_gene))

start = time.time()
while len(common_gene_2):
    res = max(common_gene_2, key=len, default='')
    #print(res)
    num = common_gene.index(res)
    select_gene_combined = select_gene_combined.replace(str(res), str(num))
    common_gene_2.remove(res)
end = time.time()
print('Running time:', end - start)

with open('E:\\shiyan_data\\GBM\\lab\\String\\select_gene_index.txt', 'w') as f:
    f.write(select_gene_combined)

#将string网络处理为一个矩阵
with open('E:\\shiyan_data\\GBM\\lab\\String\\select_gene_index.txt', 'r') as f:
    stringV10_sub = f.readlines()

n = len(common_gene)
num = np.ndarray(shape=(n, n), dtype=float)
for i in range(1, len(stringV10_sub)):
    row = int(stringV10_sub[i].split('\t')[0])
    col = int(stringV10_sub[i].split('\t')[1])
    translate = float(stringV10_sub[i].split('\t')[2])/1000
    # if float(stringV10_sub[i].split('\t')[2]) > 0.0:
    #     translate = 1
    # else:
    #     translate = 0
    num[row][col] = translate
    num[col][row] = translate

with open('E:\\shiyan_data\\GBM\\lab\\String\\string_sub_matrix.txt', 'w') as f:
    for i in range(len(num)):
        f.write(str(num[i]))
        f.write('\n')

with open('E:\\shiyan_data\\GBM\\lab\\String\\string_sub_matrix.bin', 'wb') as f:
    pickle.dump(num, f, protocol=1)
