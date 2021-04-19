from numpy import zeros
import random
from math import log10


def find_when_zero(N, m, initial_A, p):
    """
    Function returns an m x len(A) array with 1's where the recessive allele
    frequency is zero. For example, if return_value[1][2] == 1, then that means
    on generation 2, gene 3, there were no recessive alleles.
    """
    total_alleles = p*N
    A = initial_A.copy()
    when_zero = zeros((m, len(A)))

    for generation in range(0, m):
        # calculate all the alleles for the next generation
        alleles = zeros((len(A), total_alleles))
        for i in range(len(A)):  # each gene
            prob_for_recessive = A[i] / total_alleles
            for j in range(total_alleles):
                # if the rng appears lower than prob_for_recessive
                if prob_for_recessive > random.random():
                    alleles[i][j] = 1  # set to recessive
                else:
                    alleles[i][j] = 0  # set to dominant

        for k in range(len(A)):  # set the new A
            # set num of recessive to the appropriate value in A
            A[k] = sum(alleles[k])
            if A[k] == 0:  # check if the allele freq is zero
                when_zero[generation][k] = 1

    return when_zero


def allele_frequency_sim(N, m, A, p):
    """
    Simulates each allele over multiple generations.
    :param N: population size
    :param m: number of generations
    :param A: recessive alleles per gene
    :param p: ploidy
    :return: log10 of the chances the recessive alleles will disappear
    """
    iterations = 10000
    when_zero_sum = zeros((m, len(A)))

    for iter in range(iterations):
        # generate an array mapping out where the allele freq's are
        # zero in terms of generation and gene
        when_zero = find_when_zero(N, m, A, p)

        # do matrix addition to collect the sum, will avg later
        for x in range(m):
            for y in range(len(A)):
                when_zero_sum[x][y] += when_zero[x][y]

    # return the average of all the when_zero matrices
    for i in range(m):
        for j in range(len(A)):
            when_zero_sum[i][j] = log10(when_zero_sum[i][j]/iterations)
    return when_zero_sum


def experiment_1():
    N = 4
    m = 3
    A = [0, 1, 2]
    output = allele_frequency_sim(N, m, A, 2)
    for i in output:
        print(i)


if __name__ == '__main__':
    experiment_1()