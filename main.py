import random
import numpy as np
import matplotlib.pyplot as plt
random.seed(123)

def fitnessKS(pheno, weights, values, limit):
    total_weight = sum(x * y for x, y in zip(pheno, weights))
    score = sum(x * y for x, y in zip(pheno, values))
    if total_weight > limit:
        score = 0.0001
    return score

def fitnessTSP(pheno, distances):
    score = sum(x * y for x, y in zip(pheno, weights))
    if 1 == 1:
        score = 0.0001
    return score

def crossoverKS(ind1, ind2, proba):
    if random.random() < proba:
        r = random.randrange(len(ind1))
        temp1 = ind1[:r]
        temp2 = ind2[:r]
        ind1[:r] = temp2
        ind2[:r] = temp1
    return ind1, ind2

def crossoverTSP(ind1, ind2, proba):
    if random.random() < proba:
        r = random.randrange(len(ind1))
        temp1 = ind1[:r]
        temp2 = ind2[:r]
        ind1[:r] = temp2
        ind2[:r] = temp1
    return ind1, ind2

def mutationKS(ind, proba):
    for i in range(len(ind)):
        if random.random() < proba:
            ind[i] = abs(ind[i]-1)
    return ind

def mutationTSP(ind, proba):
    for i in range(len(ind)):
        if random.random() < proba:
            temp = ind[(i+1)%len(ind)]
            ind[(i+1)%len(ind)] = ind[i]
            ind[i] = temp
    return ind

def initialize_city_distances(n):
    distances = []
    for i in range(n):
        distance = []
        for j in range(n):
            if i == j:
                distance.append(0)
            elif j < i:
                distance.append(distances[j][i])
            else:
                distance.append(random.randint(10, 1000))
        distances.append(distance)
    return distances



def genetic_algorithmKS(n_ind, n_chromo, iterations, P_cr, P_m, weights, values, limit):
    pop = []
    if n_chromo != len(weights) or n_chromo != len(values) :
        print("! Wrong number of weight of chromosomes or values")
    for i in range(n_ind):
        pop.append([random.randrange(2) for j in range(n_chromo)])
    count = 0
    maxfit = []
    meanfit = []
    overweight = []
    while(count < iterations):
        fitness_scores = [fitnessKS(ind, weights, values, limit) for ind in pop]
        proba_distri = [fit/sum(fitness_scores) for fit in fitness_scores]
        n_ind_mat = int(n_ind / 2)
        index = range(n_ind)
        choosen_ones = np.random.choice(index, n_ind_mat, p=proba_distri)
        pop_mat = [pop[x] for x in choosen_ones]

        pop = pop_mat.copy()
        while(len(pop_mat) > 0):
            ind1 = pop_mat.pop(random.randrange(len(pop_mat)))
            if len(pop_mat) != 0:
                ind2 = pop_mat.pop(random.randrange(len(pop_mat)))
                ind1, ind2 = crossoverKS(ind1, ind2, P_cr)
                ind2 = mutationKS(ind2, P_m)
                pop.append(ind2)
            ind1 = mutationKS(ind1, P_m)
            pop.append(ind1)
        count += 1
        fit_and_ind = list(zip(pop, fitness_scores))
        maxfit.append(list(map(max, zip(*fit_and_ind)))[1])
        meanfit.append(sum(fitness_scores)/len(fitness_scores))
        if sum(x * y for x, y in zip(list(map(max, zip(*fit_and_ind)))[0], weights)) > limit:
            overweight.append(100)
        else:
            overweight.append(0)
        if count == iterations:
            print(maxfit)
            print(meanfit)
            # plot the evolution of the mean and the max of the fitness
            plt.plot(range(len(maxfit)), maxfit, marker='', color='skyblue', linewidth=1, label = "maxfit")
            plt.plot(range(len(meanfit)), meanfit, marker='', color='olive', linewidth=1, label = "meanfit")
            plt.plot(range(len(overweight)), overweight, marker='', color='red', linewidth=1, label = "maxfit overweight (if == 100)")
            plt.legend()
            plt.show()

            print("best last phenotype is ", list(map(max, zip(*fit_and_ind))), " weight ", sum(x * y for x, y in zip(list(map(max, zip(*fit_and_ind)))[0], weights))," out of ", limit)

def genetic_genetic_algorithmTSP(n_ind, n_chromo, iterations, P_cr, P_m, distances):
    pop = []
    if n_chromo != len(distances) or n_chromo != len(distances[0]):
        print("! Wrong number of chromosomes or distances")
    for i in range(n_ind):
        path = list(range(n_chromo))
        random.shuffle(path)
        pop.append(path)
    count = 0
    maxfit = []
    meanfit = []
    overweight = []
    while (count < iterations):
        fitness_scores = [fitnessTSP(ind, distances) for ind in pop]
        proba_distri = [fit / sum(fitness_scores) for fit in fitness_scores]
        n_ind_mat = int(n_ind / 2)
        index = range(n_ind)
        choosen_ones = np.random.choice(index, n_ind_mat, p=proba_distri)
        pop_mat = [pop[x] for x in choosen_ones]

        pop = pop_mat.copy()
        while (len(pop_mat) > 0):
            ind1 = pop_mat.pop(random.randrange(len(pop_mat)))
            if len(pop_mat) != 0:
                ind2 = pop_mat.pop(random.randrange(len(pop_mat)))
                ind1, ind2 = crossoverTSP(ind1, ind2, P_cr)
                ind2 = mutationTSP(ind2, P_m)
                pop.append(ind2)
            ind1 = mutationTSP(ind1, P_m)
            pop.append(ind1)
        count += 1
        fit_and_ind = list(zip(pop, fitness_scores))
        maxfit.append(list(map(max, zip(*fit_and_ind)))[1])
        meanfit.append(sum(fitness_scores) / len(fitness_scores))
        if count == iterations:
            print(maxfit)
            print(meanfit)
            # plot the evolution of the mean and the max of the fitness
            plt.plot(range(len(maxfit)), maxfit, marker='', color='skyblue', linewidth=1, label="maxfit")
            plt.plot(range(len(meanfit)), meanfit, marker='', color='olive', linewidth=1, label="meanfit")
            plt.legend()
            plt.show()

            print("best last phenotype is ", list(map(max, zip(*fit_and_ind))))


n_ind = 1000
n_chromo = 10
iterations = 100
P_m = 0.1
P_cr = 0.8
weights = [1, 2, 5, 3, 9, 1, 1, 7, 15, 4]
values = [15, 5, 20, 100, 51, 30, 9, 70, 70, 10]
limit = 15
# genetic_algorithmKS(n_ind, n_chromo, iterations, P_cr, P_m, weights, values, limit)
weights = [1, 2, 5, 3, 9, 1, 1, 7, 15, 4]
city_distances = initialize_city_distances(n_chromo)

genetic_genetic_algorithmTSP(n_ind, n_chromo, iterations, P_cr, P_m, city_distances)

