from math import log10
from scipy import stats as sp
from numpy import zeros

def wright_fisher(N, m, A):
    B = zeros((m, len(A)))
    # iterating through alleles in A
    for i in range(0, len(A)):
        # expected number of times allele appears in the next generation
        # for each potential value of the allele frequency in the previous gen
        prev_distribution = []

        # calculate first generation
        for j in range(0, (2*N) + 1):
            prev_distribution.append(sp.binom.pmf(j, 2*N, A[i]/(2*N)))
        B[0][i] = log10(prev_distribution[0])

        # calculate subsequent generations, using prev_distribution as p in
        # Wright-Fisher formula (2N choose k) * p^k * (p-1)^(2N-k)
        for j in range(1, m):
            new_distr = []
            for x in range(0, (2*N) + 1):
                new_gen = []
                new_distr.append(0)
                # generate list of probabilities of an allele showing up x times
                # if it showed up y times last generation
                for y in range(0, (2*N)+1):
                    new_gen.append(sp.binom.pmf(x, 2*N, y/(2*N)))
                # apply previous generation's probability distribution
                for y in range(0, 2*N):
                    new_distr[-1] += new_gen[y]*prev_distribution[y]

            prev_distribution = new_distr
            B[j][i] = log10(prev_distribution[0])
    return B

# p = ploidy
def wright_fisher_polyploid(N, m, A, p):
    B = zeros((m, len(A)))
    # iterating through alleles in A
    for i in range(0, len(A)):
        # expected number of times allele appears in the next generation
        # for each potential value of the allele frequency in the previous gen
        prev_distribution = []

        # calculate first generation
        for j in range(0, (p*N) + 1):
            prev_distribution.append(sp.binom.pmf(j, p*N, A[i]/(p*N)))
        B[0][i] = log10(prev_distribution[0])

        # calculate subsequent generations, using prev_distribution as p in
        # Wright-Fisher formula (ploid*N choose k) * p^k * (p-1)^(ploid*N-k)
        for j in range(1, m):
            new_distr = []
            for x in range(0, (p*N) + 1):
                new_gen = []
                new_distr.append(0)
                # generate list of probabilities of an allele showing up x times
                # if it showed up y times last generation
                for y in range(0, (p*N)+1):
                    new_gen.append(sp.binom.pmf(x, p*N, y/(p*N)))
                # apply previous generation's probability distribution
                for y in range(0, p*N):
                    new_distr[-1] += new_gen[y]*prev_distribution[y]

            prev_distribution = new_distr
            B[j][i] = log10(prev_distribution[0])
    return B


if __name__ == '__main__':
    N = 4
    m = 3
    A = [0, 1, 2]

    print(wright_fisher(N,m,A))
    print(wright_fisher_polyploid(N, m, A, 3))
