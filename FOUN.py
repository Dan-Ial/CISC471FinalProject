"""
File for implementation of wright fisher founder effect
    Bennet Montgomery
    Evelyn Yach
    Daniel Oh
"""

from math import log10
from scipy import stats as sp
from WFMD import WFMD
from numpy import zeros
from FOUN_test import *

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

if __name__ == '__main__':
    unittest.main()