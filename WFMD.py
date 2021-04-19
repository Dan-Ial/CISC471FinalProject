"""
File for implementation of solution to the WFMD predecessor to FOUN problem
    Bennet Montgomery
    Daniel Oh
    Evelyn Yach
"""

from scipy import stats as sp


def WFMD(N, m, g, k):
    """
    Implementation of a solution for the WFMD problem. The algorithm takes an
    organism count and allele incidence and then applies wright fisher with
    all related assumptions to these initial conditions. WF is extended across
    multiple generations by taking the expected probabilities of each generation
    and using the distribution of probability distributions to calculate each
    generation.

    :param N: number of organisms (constant across generations)
    :param m: allele frequency
    :param g: number of generations
    :param k: target incidence
    :return: probability of the allele having an incidence of at least k in the
    population between 0 and 1.
    """
    prev_distribution = []

    for j in range(0, (2 * N) + 1):
        prev_distribution.append(sp.binom.pmf(j, 2 * N, ((2*N) - m) / (2 * N)))

    for j in range(1, g):
        new_distr = []
        for x in range(0, (2 * N) + 1):
            new_gen = []
            new_distr.append(0)
            for y in range(0, (2 * N) + 1):
                new_gen.append(sp.binom.pmf(x, 2 * N, y / (2 * N)))
            for y in range(0, 2 * N):
                new_distr[-1] += new_gen[y] * prev_distribution[y]

        prev_distribution = new_distr

    result = 0

    for i in range(0, k):
        result += prev_distribution[i]

    return 1 - result

if __name__ == '__main__':
    print(WFMD(4, 6, 2, 1))