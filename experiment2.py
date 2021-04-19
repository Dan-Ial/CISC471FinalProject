from main import wright_fisher_founder
from main import wright_fisher_founder_polyploid

#mark 80 characters v
#1234567890123456789012345678901234567890123456789012345678901234567890123456789

def run_experiment2():
    """
        This is the main function for running the second experiment

        There are two options for running the series of tests that compromise,
        the second experiment. The first option (default) will run all of the
        currently defined tests. The second option will run one specific tests.
    """
    multipleRuns = True  # run multiple B populations

    if multipleRuns:  # run all the tests
        numRuns = 1  # this number can be between 0 and _
        for i in range(0, numRuns+1):
            popA, popB = load_tests(i)
            resultA = wright_fisher_founder(popA[0], popA[1], popA[2])
            resultB = wright_fisher_founder(popB[0], popB[1], popB[2])
            display_results(resultA, resultB, popA, popB, i)

    else:  # run one specific one (debugging/data entry)
        testNum = 0  # change this to the B pop. you want to run
        popA, popB = load_tests(testNum)
        resultA = wright_fisher_founder(popA[0], popA[1], popA[2])
        resultB = wright_fisher_founder(popB[0], popB[1], popB[2])
        display_results(resultA, resultB, popA, popB, testNum)

def display_results(rA, rB, pA, pB, whichB):
    """
        This function handles displaying the output of each test on the console

        Begins by displaying which B population test parameters are being used
        (denoted by an integer), next display the relevant info for population
        A and then display the info for population B.

        :param rA: The results of population A
        :param rB: The results of population B
        :param pA: The initial parameters used for population A
        :param pB: The initial parameters used for population B
        :param whichB: which parameters for B population
        :type rA: list
        :type rB: list
        :type pA: list
        :type pB: list
        :type whichB: int
    """

    # print out initial statement
    print("Using population B number:", whichB)
    # print info pertaining to the control population (A)
    print("Results from the control population:")
    print("parameters: N =", pA[0], " m =", pA[1], " A =", pA[2])
    print(rA[0], "\n", rA[1], "\n", rA[2])
    # print info pertaining to the varying population (B)
    print("Results from the varying population:")
    print("parameters: N =", pB[0], " m =", pB[1], " A =", pB[2])
    print(rB[0], "\n", rB[1], "\n", rB[2])
    # newline character for formatting
    print("\n")

def load_tests(whichB):
    """
        This function handles setting up all the test parameters for the
        experiment

        Each population below is stored as an array of parameters:
            population = [N, m, A]
        where N is the # of individuals in a population, m is the # of
        generations to calculate, and A an array containing the # of
        recessive alleles for each factor in the population

        :param whichB: which parameters for B population
        :type whichB: int
        :return: The list of parameters for population A and B
        :rtype: list
    """

    # define CONTROL, or population A
    # this population will model an average population
    popA = [5, 3, [5]]

    # define VARIANCE, or population B
    if (whichB == 0):
        # this population will model a with no recessive alleles
        popB = [5, 3, [0]]
    elif (whichB == 1):
        # this population will model an high recessive alleles
        popB = [5, 3, [9]]

    return popA, popB


if __name__ == '__main__':
    run_experiment2()
