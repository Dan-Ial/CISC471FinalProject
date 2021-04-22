from math import log10
from scipy import stats as sp
from numpy import zeros

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