#encoding:utf-8
from numpy import *
import numpy as np
import pandas as pd
import pickle
import time

def Get_gene():
    with open('E:\shiyan_data\GBM\lab\A.txt', 'r') as f:
        GBM_mua = f.readlines()

    print(GBM_mua[0].split('\t'))
    print(len(GBM_mua[0].split('\t')))

    GBM_mua_gene_sub = []
    GBM_mua_sub = []
    GBM_mua_sub_total_col = []
    for i in range(1, len(GBM_mua[0].split('\t'))):
        GBM_mua_col = [col.split('\t')[i] for col in GBM_mua]
        if np.sum(np.array(GBM_mua_col[1:], dtype=int)) >= 0.4:
            GBM_mua_gene_sub.append(GBM_mua_col[0])
            GBM_mua_sub.append(GBM_mua_col[1:])
            GBM_mua_sub_total_col.append(np.sum(np.array(GBM_mua_col[1:], dtype=int)))

    print(GBM_mua_gene_sub)
    print(len(GBM_mua_gene_sub))

    with open('E:\shiyan_data\GBM\lab\Common_gene_diff_exp.bin', 'wb') as f:
        pickle.dump(GBM_mua_gene_sub, f, protocol=1)

    with open('E:\shiyan_data\GBM\lab\Common_gene_diff_exp.txt', 'w') as f:
        f.write(str(GBM_mua_gene_sub))

    with open('E:\shiyan_data\GBM\lab\GBM_mua_sub.bin', 'wb') as f:
        pickle.dump(GBM_mua_sub, f, protocol=1)

    with open('E:\shiyan_data\GBM\lab\GBM_mua_sub.txt', 'w') as f:
        f.write(str(GBM_mua_sub))

    with open('E:\shiyan_data\GBM\lab\exchange_num_col_total_diff_exp.bin', 'wb') as f:
        pickle.dump(GBM_mua_sub_total_col, f, protocol=1)

    return GBM_mua_gene_sub, GBM_mua_sub

def Get_string():
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
        if (gene_texting[i].split('\t')[0] == gene_experiment[i].split('\t')[0]) and (
                gene_texting[i].split('\t')[1] == gene_experiment[i].split('\t')[1]):
            if int(gene_experiment[i].split('\t')[2].split('\n')[0]) != 0:
                gene_num = max(gene_texting[i].split('\t')[2].split('\n')[0],
                               gene_experiment[i].split('\t')[2].split('\n')[0])
            else:
                gene_num = 0
            replace_gene.append(
                gene_texting[i].split('\t')[0] + '\t' + gene_texting[i].split('\t')[1] + '\t' + str(gene_num) + '\n')

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
        # print(res)
        num = common_gene.index(res)
        select_gene_combined = select_gene_combined.replace(str(res), str(num))
        common_gene_2.remove(res)
    end = time.time()
    print('Running time:', end - start)

    with open('E:\\shiyan_data\\GBM\\lab\\String\\select_gene_index.txt', 'w') as f:
        f.write(select_gene_combined)

    # 将string网络处理为一个矩阵
    with open('E:\\shiyan_data\\GBM\\lab\\String\\select_gene_index.txt', 'r') as f:
        stringV10_sub = f.readlines()

    n = len(common_gene)
    num = np.ndarray(shape=(n, n), dtype=float)
    for i in range(1, len(stringV10_sub)):
        row = int(stringV10_sub[i].split('\t')[0])
        col = int(stringV10_sub[i].split('\t')[1])
        translate = float(stringV10_sub[i].split('\t')[2]) / 1000
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

def get_Gene_num(GBM_mua_gene_sub, num_1, num_2):

    geng_num_1 = GBM_mua_gene_sub[num_1]
    geng_num_2 = GBM_mua_gene_sub[num_2]

    return geng_num_1, geng_num_2

#计算海明距离
def Hamming_distance(s1, s2):
    sum_1 = 0
    sum_2 = 0
    num_1 = 0
    num_2 = 0
    cover = 0
    assert len(s1) == len(s2)  #判断两列是否长度相等
    for i in range(len(s1)):
        if s1[i] == '1' and s2[i] == '0':
            sum_1 += 1
        if s1[i] == '0' and s2[i] == '1':
            sum_2 += 1
        if s1[i] == '1':
            num_1 += 1
        if s2[i] == '1':
            num_2 += 1

    cover = num_1 + num_2 - (num_1 - sum_1)
    return float(((sum_1/num_1) + (sum_2/num_2))/2), cover

#获取每个基因与其他基因的海明距离
def Count_HM(GBM_mua_gene_sub, GBM_mua_sub):
    dict_par = {}
    gene_network = []
    Num_gene = len(GBM_mua_gene_sub)
    hm_mat = np.ndarray(shape=(Num_gene, Num_gene))    #9483*9483
    cover_mat = np.ndarray(shape=(Num_gene, Num_gene))

    for i in range(Num_gene):
        print(i)
        exchange_line_0 = GBM_mua_sub[i]

        for j in range(Num_gene):
            exchange_line_1 = GBM_mua_sub[j]
            #计算每两列之间的海明距离
            sum_hm, cover = Hamming_distance(exchange_line_0, exchange_line_1)
            #print(sum_hm)
            hm_mat[i][j] = sum_hm  #将海明距离写入矩阵
            cover_mat[i][j] = cover
            gene_num1, gene_num2 = get_Gene_num(GBM_mua_gene_sub, i, j)
            dict_par[str(gene_num1) + ',' + str(gene_num2)] = hm_mat[i][j]
            gene_network_list = str(gene_num1) + '\t' + str(gene_num2) + '\t' + str(hm_mat[i][j]) + '\n'
            gene_network.append(gene_network_list)
    key_sorted = sorted(zip(dict_par.values(), dict_par.keys()), reverse=True)

    #输出海明距离矩阵、二进制格式
    with open('E:\shiyan_data\GBM\lab\hm_data_example.bin', 'wb') as f:
        pickle.dump(hm_mat, f, protocol=1)

    with open('E:\shiyan_data\GBM\lab\hm_data_example.txt', 'w') as f:
        for i in range(len(hm_mat)):
            f.write(str(hm_mat[i]))
            f.write('\n')

    #输出覆盖度矩阵、二进制格式
    # with open('E:\shiyan_data\BRCA\cover_data_example.bin', 'wb') as f:
    #     pickle.dump(cover_mat, f, protocol=1)
    #
    # with open('E:\shiyan_data\BRCA\cover_data_example.txt', 'w') as f:
    #     for i in range(len(cover_mat)):
    #         f.write(str(cover_mat[i]))
    #         f.write('\n')

    #将获得的海明距离写入文件
    with open('E:\shiyan_data\GBM\lab\Hm_gene_example.txt', 'w') as f:
        f.write(str(key_sorted))

    # with open('E:\shiyan_data\BRCA\Gene_network_mutex.txt', 'w') as f:
    #     for i in range(len(gene_network)):
    #         f.write(gene_network[i])

if __name__ == '__main__':
    # GBM_mua_gene_sub, GBM_mua_sub = Get_gene()
    # Count_HM(GBM_mua_gene_sub, GBM_mua_sub)
    Get_string()
