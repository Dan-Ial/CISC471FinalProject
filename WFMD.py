from scipy import stats as sp

# WFMD function
def WFMD(N, m, g, k):
    # expected number of times allele appears in the next generation
    # for each potential value of the allele frequency in the previous gen
    prev_distribution = []

    # calculate first generation
    for j in range(0, (2 * N) + 1):
        prev_distribution.append(sp.binom.pmf(j, 2 * N, ((2*N) - m) / (2 * N)))

    # calculate subsequent generations, using prev_distribution as p in
    # Wright-Fisher formula (2N choose k) * p^k * (p-1)^(2N-k)
    for j in range(1, g):
        new_distr = []
        for x in range(0, (2 * N) + 1):
            new_gen = []
            new_distr.append(0)
            # generate list of probabilities of an allele showing up x times
            # if it showed up y times last generation
            for y in range(0, (2 * N) + 1):
                new_gen.append(sp.binom.pmf(x, 2 * N, y / (2 * N)))
            # apply previous generation's probability distribution
            for y in range(0, 2 * N):
                new_distr[-1] += new_gen[y] * prev_distribution[y]

        prev_distribution = new_distr

    result = 0

    for i in range(0, k):
        result += prev_distribution[i]

    return 1 - result

if __name__ == '__main__':
    print(WFMD(4, 6, 2, 1))