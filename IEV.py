from rosalind_reader import read_file

def iev(couples):
    """
    the following are the chances for 1 dominant offspring
    AA-AA: 100%
    AA-Aa: 100%
    AA-aa: 100%
    Aa-Aa: 75%
    Aa-aa: 50%
    aa-aa: 0%
    :param couples: Six nonnegative integers. The integers correspond to the
        number of couples in a population possessing each genotype pairing, as
        listed above
    :return: The expected number of offspring displaying the dominant phenotype
        in the next generation, under the assumption that every couple has
        exactly two offspring.
    """
    # doesnt consider couples[5] since it's 0% every time
    return couples[0]*2 + couples[1]*2 + couples[2]*2 + \
           3/4*couples[3]*2 + 1/2*couples[4]*2


if __name__ == '__main__':
    data = [int(num) for num in read_file('iev.txt')]
    print(iev(data))