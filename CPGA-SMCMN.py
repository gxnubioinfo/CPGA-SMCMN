import numpy as np
import pickle
import time
import random

SELECT1 = 0.7
SELECT2 = 0.5
N_GENERATIONS = 1000
MUTATION_RATE = 0.03

def Clustering(hm_data, string_sub_matrix, select_num1=0.9, select_num2=0.8):

    #print('select_num1:', select_num1)
    cluster_index = []
    for i in range(len(hm_data)):
        neighbours = np.where((hm_data[:, i] > select_num1))[0]
        #neighbours = np.insert(neighbours, 0, i)
        neighbours = neighbours.tolist()
        cluster_index.append(neighbours)

    string_cluster = []
    for j in range(len(string_sub_matrix)):
        string_neighbours = np.where((string_sub_matrix[:, j] > select_num2))[0]
        #string_neighbours = np.insert(string_neighbours, 0, j)
        string_neighbours = string_neighbours.tolist()
        string_cluster.append(string_neighbours)


    return cluster_index, string_cluster

def Get_random_index_2(common_gene_num_total, cluster_index_sub):
    common_gene_num_total_sub_list = []
    for i in range(len(cluster_index_sub)):
        num = common_gene_num_total[cluster_index_sub[i]]
        common_gene_num_total_sub_list.append(num)

    common_gene_num_total_sub = np.array(common_gene_num_total_sub_list)
    data_index = np.random.choice(cluster_index_sub, size=1, replace=False, p=(common_gene_num_total_sub) / (common_gene_num_total_sub.sum()))[0]

    return data_index

def Select_gene_index(cluster_index, string_cluster_index, string_sub_matrix, exchange_num_col_total_diff_exp, K=3, select_num = 0.7):
    #print(K)
    n = len(cluster_index)

    order_list = []

    isProcess = np.full((n,), -1)

    cluster_index_total = []
    for i in range(len(cluster_index)):
        num = len(cluster_index[i])
        cluster_index_total.append(num)

    cluster_index_total = np.array(cluster_index_total)
    gene_cluster_num = np.random.choice(np.arange(len(cluster_index)), size=1, replace=False, p=(cluster_index_total)/(cluster_index_total.sum()))[0]

    #  = random.randint(0, int(len(cluster_index) - 1))
    order_list.append(gene_cluster_num)
    isProcess[gene_cluster_num] = 1

    flag = len(cluster_index[gene_cluster_num])
    while flag > 0:
        #gene_index_num = Get_big_index(common_gene_num_total, cluster_index[gene_cluster_num])
        gene_num = Get_random_index_2(exchange_num_col_total_diff_exp, cluster_index[gene_cluster_num])
        if string_sub_matrix[gene_cluster_num][gene_num] > select_num:
            order_list.append(gene_num)
            isProcess[gene_num] = 1
            break
        flag -= 1

    if len(order_list) == K:
        #print('?????????????????????')
        return order_list

    #union_list = list(set(cluster_index[gene_cluster_num]).union(set(cluster_index[gene_num])))
    string_uion_list = list(set(string_cluster_index[gene_cluster_num]).union(set(string_cluster_index[gene_num])))
    intersection_list = [i for i in cluster_index[gene_cluster_num] if i in cluster_index[gene_num]]
    intersection_list = [i for i in string_uion_list if i in intersection_list]

    if len(intersection_list) == 0:
        return []

    gene_num = Get_random_index_2(exchange_num_col_total_diff_exp, intersection_list)

    stop_m = 0
    while (np.sum(isProcess) != n):
        flag = len(order_list)

        if set(order_list) < set(cluster_index[gene_num]) and isProcess[gene_num] == -1:
            distance_gene_list = []
            for i in range(len(order_list)):
                distance_gene = string_sub_matrix[order_list[i]][gene_num]
                distance_gene_list.append(distance_gene)
            #print('distance_gene_list:', distance_gene_list)
            if (sum(distance_gene_list)/len(distance_gene_list)) > select_num: # any  all
                isProcess[gene_num] = 1
                order_list.append(gene_num)
                string_uion_list = list(set(string_uion_list).union(set(string_cluster_index[gene_num])))
                intersection_list = [i for i in intersection_list if i in cluster_index[gene_num]]
                intersection_list = [i for i in string_uion_list if i in intersection_list]
                if len(intersection_list) == 0:
                    return []
                gene_num = Get_random_index_2(exchange_num_col_total_diff_exp, intersection_list)
                #gene_index_num = Get_big_index(common_gene_num_total, union_list)
                #gene_index_num = random.randint(0, int(len(union_list) - 1))
            else:
                gene_num = Get_random_index_2(exchange_num_col_total_diff_exp, intersection_list)
                #gene_index_num = Get_big_index(common_gene_num_total, union_list)
                #gene_index_num = random.randint(0, int(len(union_list) - 1))

        else:
            isProcess[gene_num] = 1
            gene_num = Get_random_index_2(exchange_num_col_total_diff_exp, intersection_list)
            #gene_index_num = Get_big_index(common_gene_num_total, union_list)
            #gene_index_num = random.randint(0, int(len(union_list) - 1))

        if len(order_list) == K:
            #print('?????????????????????')
            return order_list

        if len(order_list) == flag:
            stop_m += 1
        else:
            stop_m = 0

        if stop_m >= 920:
            break

    return []

def F(x, y, z):
    return x + y + z
    #return x + y           # ?????????????????????

def get_fitness(pop, gene_patient_indices_select, hm_data, string_sub_matrix):
    x, y, z = translateDNA(pop, gene_patient_indices_select, hm_data, string_sub_matrix)
    pred = F(x, y, z)
    return pred

# pop???????????????????????????????????????DNA?????????????????????????????????
def translateDNA(pop, gene_patient_indices_select, hm_data, string_sub_matrix):
# ??????????????????x???y???
    X = []
    Y = []
    Z = []
    for i in pop:
        y = 0
        z = 0
        union_list = []

        # ?????????????????????????????????
        pop_i = len(i)
        i_n = 0
        j_n = 1
        while True:
            if j_n == pop_i:
                i_n = i_n + 1
                j_n = i_n + 1
            if i_n == pop_i - 1:
                break
            z = z + (string_sub_matrix[i[i_n]][i[j_n]] + string_sub_matrix[i[j_n]][i[i_n]])/2
            y = y + (hm_data[i[i_n]][i[j_n]] + hm_data[i[j_n]][i[i_n]])/2
            j_n = j_n + 1

        y = (y/((pop_i * (pop_i - 1))/2))
        z = z/((pop_i * (pop_i - 1))/2)

        # ?????????????????????????????????
        for j in i:
            gene_patient_indices_select_1 = gene_patient_indices_select[j].split('\n')
            gene_patient_indices_select_list = gene_patient_indices_select_1[0].split('\t')[1:]
            # print(gene_patient_indices_select_list)
            # time.sleep(20)
            union_list = list(set(gene_patient_indices_select_list).union(set(union_list)))

        x = (len(union_list)/48) # ?????????48????????????8??????75
        X.append(x)
        Y.append(y)
        Z.append(z)

    #edge = np.array(edge)
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)

    return X, Y, Z

def mutation(pop, cluster_index, gene_patient_indices_select, string_sub_matrix, hm_data):
    new_pop = []

    # ????????????????????????????????????????????????????????????
    for father in pop:
        # ????????????????????????????????????
        child = father
        #print('????????????', child)
        # ??????????????????????????????????????????
        flag = np.random.randint(0, 2)
        # if np.random.rand() < MUTATION_RATE:
        if flag == 0:
            mutation_1(child, cluster_index, gene_patient_indices_select, string_sub_matrix, hm_data)
        else:
            mutation_2(child, cluster_index, gene_patient_indices_select, string_sub_matrix, hm_data)

        #print('????????????', child)
        new_pop.append(child)


    return new_pop

def mutation_1(child, cluster_index, gene_patient_indices_select, string_sub_matrix, hm_data):
    # ???child????????????????????????????????????????????????????????????????????????
    mother = random.sample(child, 1)[0]
    child.remove(mother)
    alternate = cluster_index[child[0]]
    for i in child:
        alternate = list(set(alternate).intersection(cluster_index[i]))
    alternate.remove(mother)

    # union_list = []
    # for i in child:
    #     gene_patient_indices_select_1 = gene_patient_indices_select[i].split('\n')
    #     gene_patient_indices_select_list = gene_patient_indices_select_1[0].split('\t')[1:]
    #     union_list = list(set(gene_patient_indices_select_list).union(set(union_list)))

    sum_total = []
    for i in alternate:
        gene_patient_indices_select_1 = gene_patient_indices_select[i].split('\n')
        gene_patient_indices_select_list = gene_patient_indices_select_1[0].split('\t')[1:]
        sum_total.append(len(gene_patient_indices_select_list))

    mother_mutation_point = sum_total.index(max(sum_total))
    while True:
        # mother_sample = gene_patient_indices_select[alternate[mother_mutation_point]]
        # mother_sample = mother_sample.split('\n')
        # mother_sample = mother_sample[0].split('\t')[1:]
        # mother_sample_sub = [i for i in mother_sample if i not in union_list]
        # union_list_sub = [j for j in union_list if j not in mother_sample]
        # #print('mutation_1:', len(mother_sample), len(union_list))
        # EX = ((len(mother_sample_sub)/len(mother_sample)) + (len(union_list_sub)/len(union_list)))/2

        EX = []
        for i in child:
            ME = hm_data[i][alternate[mother_mutation_point]]
            EX.append(ME)

        string_dictance = []
        for i in child:
            distance = string_sub_matrix[i][alternate[mother_mutation_point]]
            string_dictance.append(distance)
        #print('string_dictance:', string_dictance)

        if alternate[mother_mutation_point] not in child and (sum(EX)/len(EX)) > SELECT1 and (sum(string_dictance)/len(string_dictance)) > SELECT2:
            child.append(alternate[mother_mutation_point])
            break
        else:
            sum_total[mother_mutation_point] = 1
            mother_mutation_point = sum_total.index(max(sum_total))

        if (np.array(sum_total) == 1).all():
            child.append(mother)
            break

def mutation_2(child, cluster_index, gene_patient_indices_select, string_sub_matrix, hm_data):
    # ???child????????????????????????????????????????????????????????????????????????
    sum_total = []
    for i in child:
        gene_patient_indices_select_1 = gene_patient_indices_select[i].split('\n')
        gene_patient_indices_select_list = gene_patient_indices_select_1[0].split('\t')[1:]
        sum_total.append(len(gene_patient_indices_select_list))
    child_mutation_point = sum_total.index(min(sum_total))
    mother = child[child_mutation_point]
    child.remove(mother)

    alternate = cluster_index[child[0]]
    for i in child:
        alternate = list(set(alternate).intersection(cluster_index[i]))
    alternate.remove(mother)

    # union_list = []
    # for i in child:
    #     gene_patient_indices_select_1 = gene_patient_indices_select[i].split('\n')
    #     gene_patient_indices_select_list = gene_patient_indices_select_1[0].split('\t')[1:]
    #     union_list = list(set(gene_patient_indices_select_list).union(set(union_list)))

    sum_total = []
    for i in alternate:
        gene_patient_indices_select_1 = gene_patient_indices_select[i].split('\n')
        gene_patient_indices_select_list = gene_patient_indices_select_1[0].split('\t')[1:]
        sum_total.append(len(gene_patient_indices_select_list))

    mother_mutation_point = sum_total.index(max(sum_total))
    while True:
        # mother_sample = gene_patient_indices_select[alternate[mother_mutation_point]]
        # mother_sample = mother_sample.split('\n')
        # mother_sample = mother_sample[0].split('\t')[1:]
        # mother_sample_sub = [i for i in mother_sample if i not in union_list]
        # union_list_sub = [j for j in union_list if j not in mother_sample]
        # #print('mutation_2:', len(mother_sample), len(union_list))
        # EX = ((len(mother_sample_sub)/len(mother_sample)) + (len(union_list_sub)/len(union_list)))/2

        EX = []
        for i in child:
            ME = hm_data[i][alternate[mother_mutation_point]]
            EX.append(ME)

        string_dictance = []
        for i in child:
            distance = string_sub_matrix[i][alternate[mother_mutation_point]]
            string_dictance.append(distance)
        #print('string_dictance:', string_dictance)
        #time.sleep(1)

        if alternate[mother_mutation_point] not in child and (sum(EX)/len(EX)) > SELECT1 and (sum(string_dictance)/len(string_dictance)) > SELECT2:
            child.append(alternate[mother_mutation_point])
            break
        else:
            sum_total[mother_mutation_point] = 1
            mother_mutation_point = sum_total.index(max(sum_total))

        if (np.array(sum_total) == 1).all():
            child.append(mother)
            break

def select(pop, fitness, pop_size):
    """????????????????????????????????????????????????????????????????????????????????????????????????????????????"""
    #print('?????????', pop)
    idx = np.random.choice(np.arange(pop_size), size=pop_size - 1, replace=True, p=(fitness) / (fitness.sum()))
    idx = idx.tolist()
    max_fitness_index = np.argmax(fitness)
    idx.append(max_fitness_index)
    idx = np.array(idx)
    #print('?????????', pop[idx])
    return pop[idx]

def print_info(pop, gene_patient_indices_select, hm_data, string_sub_matrix):
    global best_population
    # x, y, z = translateDNA(pop, gene_patient_indices_select, hm_data, string_sub_matrix)
    Fitness = get_fitness(pop, gene_patient_indices_select, hm_data, string_sub_matrix)

    # max_fitness_index = np.where(x == np.max(x))[0]
    max_fitness_index = np.where(Fitness == np.max(Fitness))[0]
    max_fitness_gene = []

    # print('?????????????????????', x[max_fitness_index[0]])
    print('?????????????????????', Fitness[max_fitness_index[0]])
    for i in range(len(max_fitness_index)):
        max_fitness_gene.append(pop[max_fitness_index[i]])

    for i in range(len(max_fitness_gene)):
        if i >= len(max_fitness_gene) - 1:
            break
        j = i + 1
        while j < len(max_fitness_gene):
            if set(max_fitness_gene[i]) == set(max_fitness_gene[j]):
                del max_fitness_gene[j]
                j = j - 1
            j += 1

    fitness_gene = []
    for i in max_fitness_gene:
        for j in i:
            fitness_gene.append(common_gene_diff_exp[j])

    print("?????????????????????????????????", fitness_gene)
    fitness[max_fitness_index] = 0

def Read_data():
    with open('E:\\shiyan_data\\GBM\\lab\\hm_data_example.bin', 'rb') as f:
        hm_data = pickle.load(f)

    with open('E:\\shiyan_data\\GBM\\lab\\Common_gene_diff_exp.bin', 'rb') as f:
        common_gene_diff_exp = pickle.load(f)

    with open('E:\\shiyan_data\\GBM\\lab\\String\\string_sub_matrix.bin', 'rb') as f:
        string_sub_matrix = pickle.load(f)

    with open('E:\\shiyan_data\\GBM\\lab\\exchange_num_col_total_diff_exp.bin', 'rb') as f:
        exchange_num_col_total_diff_exp = pickle.load(f)

    with open('E:\\shiyan_data\\GBM\\lab\\gene_patient_indices_select.txt', 'r') as f:
        gene_patient_indices_select = f.readlines()

    return hm_data, common_gene_diff_exp, string_sub_matrix, exchange_num_col_total_diff_exp, gene_patient_indices_select
if __name__ == '__main__':

    hm_data, common_gene_diff_exp, string_sub_matrix, exchange_num_col_total_diff_exp, gene_patient_indices_select = Read_data()

    cluster_index, string_cluster_index = Clustering(hm_data, string_sub_matrix, SELECT1, SELECT2)

    flag = True
    path_list_or = []
    n = 0
    path = 0
    start_time = time.time()
    for K in range(6, 9):
        if K <= 6:
            SELECT1 = 0.7

        best_population = {}

        # ????????????
        start_time = time.time()
        #for i in range(n):
        while path < int(len(hm_data)/4):
            # ?????????????????????
            order_list = Select_gene_index(cluster_index, string_cluster_index, string_sub_matrix, exchange_num_col_total_diff_exp, K, SELECT2)
            n += 1
            # ???????????????????????????
            # filter_gene = []
            # order_list_1 = Select_gene_index(cluster_index, string_cluster_index, string_sub_matrix, common_gene_num_total, filter_gene, K, 0.7)
            # order_list_2 = Select_gene_index(cluster_index, string_cluster_index, string_sub_matrix, common_gene_num_total, order_list_1, K, 0.7)
            # order_list = [order_list_1, order_list_2]
            if len(order_list) >= 2:
                path_list_or.append(order_list)
                path += 1

        path = 0

        max_fitness_pop = []
        break_circle = 0

        for i in range(10):
            for j in range(N_GENERATIONS):  # ??????N???
                path_list = np.array(mutation(path_list_or, cluster_index, gene_patient_indices_select, string_sub_matrix, hm_data))
                fitness = get_fitness(path_list, gene_patient_indices_select, hm_data, string_sub_matrix)
                max_fitness_index = np.where(fitness == np.max(fitness))[0]
                if set(max_fitness_pop) == set(path_list[max_fitness_index[0]]):
                    break_circle += 1
                else:
                    max_fitness_pop = path_list[max_fitness_index[0]]
                    break_circle = 0

                if break_circle == 10:
                    print_info(path_list, gene_patient_indices_select, hm_data, string_sub_matrix)
                    break
                # ????????????????????????
                pop_size = len(path_list)
                path_list = select(path_list, fitness, pop_size)
                # path_list.append
                path_list = path_list.tolist()
                #print("i, j:", i, j)
            end_time = time.time()
            #N_GENERATIONS += 200
            # ??????????????????
            print('???' + str(K) + '??????????????????????????????????????????', end_time - start_time)
            print_info(path_list, gene_patient_indices_select, hm_data, string_sub_matrix)

        path_list = []
