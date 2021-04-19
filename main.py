"""
File for implementation of wright fisher founder effect and extension
    Bennet Montgomery
    Evelyn Yach
    Daniel Oh
"""

from math import log10
from scipy import stats as sp
from WFMD import WFMD
from numpy import zeros

def wright_fisher_founder(N, m, A):
    """
    Implementation of the wright fisher founder effect algorithm using a
    solution to the WFMD problem. A solution is built by applying the WFMD
    algorithm to each element of A m times, once for each generation value.
    The outputs of these runs are used to generate the columns for B.

    :param N: number of individuals
    :param m: number of generations to calculate
    :param A: list of incidences of an allele
    :return: 2d list of log10 of probabilities that an allele is eliminated at
    or before that generation
    """
    B = zeros((m, len(A)))

    for i in range(0, len(A)):
        for g in range(1, m+1):
            prob = WFMD(N, (2*N) - A[i], g, 1)
            B[g-1][i] = log10(1 - prob)

    return B

def wright_fisher_founder_polyploid(N, m, A, p):
    """
    Extension of the wright fisher algorithm to polyploidy. Instead of direct
    calls to WFMD, which only works on diploids, this function implements a
    version of the algorithm used on WFMD but with the binomial distribution
    recalibrated to account for arbitrary ploidy.

    :param N: number of individuals
    :param m: number of generations to calculate
    :param A: list of indices of an allele
    :param p: ploidy of the organism
    :return: 2d ixj list of log10 of probabilities that j allele is eliminated
    in or before generation i.
    """
    B = zeros((m, len(A)))
    for i in range(0, len(A)):
        prev_distribution = []

        for j in range(0, (p*N) + 1):
            prev_distribution.append(sp.binom.pmf(j, p*N, A[i]/(p*N)))
        B[0][i] = log10(prev_distribution[0])

        for j in range(1, m):
            new_distr = []
            for x in range(0, (p*N) + 1):
                new_gen = []
                new_distr.append(0)
                for y in range(0, (p*N)+1):
                    new_gen.append(sp.binom.pmf(x, p*N, y/(p*N)))
                for y in range(0, p*N):
                    new_distr[-1] += new_gen[y]*prev_distribution[y]

            prev_distribution = new_distr
            B[j][i] = log10(prev_distribution[0])
    return B


if __name__ == '__main__':
    # N = 4
    # m = 3
    # A = [0, 1, 2]

    # print(wright_fisher_founder(N,m,A))
    # print(wright_fisher_founder_polyploid(N, m, A, 3))
    # unittest.main()
    print('test')