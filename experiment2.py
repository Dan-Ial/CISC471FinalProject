"""
File for implementation of experiment 2
    Bennet Montgomery
    Evelyn Yach
    Daniel Oh
"""

from main import wright_fisher_founder

def run_experiment2():
    """
        This is the main function for running the second experiment
    """
    # vary the population size in population B, compare to control population A
    print("Vary Population Size")
    run_populationsize()
    print("Vary Number of Generations")
    run_numberofgenerations()

def run_populationsize():
    """
        Runs the part of the experiment where we vary pop size There are
        two options for running the series of tests that compromise, the 2nd
        experiment. The first option (default) will run all of the currently
        defined tests. The second option will run one specific tests.
    """
    multipleRuns = True

    if multipleRuns:  # run all the tests
        numRuns = 5  # this number can be between 0 and 4
        for i in range(1, numRuns+1):
            popA, popB = test_populationsize(i)
            resultA = wright_fisher_founder(popA[0], popA[1], popA[2])
            resultB = wright_fisher_founder(popB[0], popB[1], popB[2])
            display_results(resultA, resultB, popA, popB, i)

    else:  # run one specific one (debugging/data entry)
        testNum = 0  # change this to the B pop. you want to run
        popA, popB = test_populationsize(testNum)
        resultA = wright_fisher_founder(popA[0], popA[1], popA[2])
        resultB = wright_fisher_founder(popB[0], popB[1], popB[2])
        display_results(resultA, resultB, popA, popB, testNum)

def test_populationsize(whichB):
    """
        This function handles setting up all the parameters for testing
        population size. Each population below is stored as an array of
        parameters: population = [N, m, A]
        where N is the # of individuals in a population, m is the # of
        generations to calculate, and A an array containing the # of
        recessive alleles for each factor in the population

        :param whichB: which parameters for B population
        :type whichB: int
        :return: The list of parameters for population A and B
        :rtype: list
    """
    # control pop
    popA = [5, 3, [0, 1, 2, 3, 4, 5]]
    # varying pop
    if (whichB == 1):
        popB = [10, 3, [0, 1, 2, 3, 4, 5]]
    elif (whichB == 2):
        popB = [20, 3, [0, 1, 2, 3, 4, 5]]
    elif (whichB == 3):
        popB = [40, 3, [0, 1, 2, 3, 4, 5]]
    elif (whichB == 4):
        popB = [80, 3, [0, 1, 2, 3, 4, 5]]
    elif (whichB == 5):
        popB = [160, 3, [0, 1, 2, 3, 4, 5]]

    return popA, popB

def run_numberofgenerations():
    """
        Runs the part of the experiment where we increase the number of
        generations on a fixed population size. Each population below is
        stored as an array of parameters: population = [N, m, A]
        where N is the # of individuals in a population, m is the # of
        generations to calculate, and A an array containing the # of
        recessive alleles for each factor in the population
    """
    # control pop
    popA = [5, 10, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]]

    # calculate the results
    resultA = wright_fisher_founder(popA[0], popA[1], popA[2])

    # display results
    display_results(resultA, [], popA, [], "N/A")

def display_results(rA, rB, pA, pB, whichB):
    """
        This function handles displaying the output of each test on the console

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
    if type(whichB) == int:
        print("Using population B number:", whichB)
    # print info pertaining to the control population (A)
    print("Results from the control population:")
    print("parameters: N =", pA[0], " m =", pA[1], " A =", pA[2])
    print(rA)
    # print info pertaining to the varying population (B)
    if type(whichB) == int:
        print("Results from the varying population:")
        print("parameters: N =", pB[0], " m =", pB[1], " A =", pB[2])
        print(rB)
    # newline character for formatting
    print("\n")

if __name__ == '__main__':
    run_experiment2()
