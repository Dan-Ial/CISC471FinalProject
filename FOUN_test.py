import unittest
from FOUN import *

THRESHOLD = 0.001

class FOUN_tests(unittest.TestCase):
    def test_wright_fisher_founder_Rosalind_default(self):
        """
        First unit test for the Founder effect, using the default founder
        effect dataset from the FOUN page on Rosalind.
        """

        N = 4
        m = 3
        A = [0, 1, 2]

        expctd = [[0.0, -0.463935575821, -0.999509892866],
                    [0.0, -0.301424998891, -0.641668367342],
                    [0.0, -0.229066698008, -0.485798552456]]

        real = wright_fisher_founder(N, m, A)

        for i in range(len(expctd)):
            for j in range(len(expctd[i])):
                self.assertAlmostEqual(expctd[i][j],real[i][j],delta=THRESHOLD)

    def test_wright_fisher_founder_handcalculated(self):
        """
        Second unit test for fisher founder. Based on results handcalculated by
        Bennet Montgomery that were also used to generate the hand calculated
        unit test for WFMD. The result was around but not exactly
        -0.4998, -0.3342, -0.2607 for generations 1, 2, and 3 respectively.
        """

        N = 2
        m = 3
        A = [1]

        expctd = [[-0.4998], [-0.3342], [-0.2607]]
        real = wright_fisher_founder(N, m, A)

        for i in range(len(expctd)):
            for j in range(len(expctd[i])):
                self.assertAlmostEqual(expctd[i][j],real[i][j],delta=THRESHOLD)

    def test_wright_fisher_founder_zero_output(self):
        """
        Third unit test for fisher founder. Output should be 0.0 in the case
        where an allele starts eliminated.
        """

        N = 3
        m = 1
        A = [0]

        expctd = [[0]]
        real = wright_fisher_founder(N, m, A)

        for i in range(len(expctd)):
            for j in range(len(expctd[i])):
                self.assertAlmostEqual(expctd[i][j],real[i][j],delta=THRESHOLD)

    def test_wright_fisher_polyploid_Rosalind_default(self):
        """
        First unit test for the polyploid implementation. Functionally the same
        as the original algorithm, but instead of the values of A indicating
        the recessive allele, they just indicate the proportion of allele in the
        nploid set being measured. Results needed to be hand calculated, as the
        distributions will not be the same as in a diploid. Due to the need to
        calculate by hand, only the 3rd column of rosalind data was considered.
        Test was done assuming a triploid organism.
        """
        N = 3
        m = 1
        A = [2]
        p = 3

        expctd = [[-0.98230]]

        real = wright_fisher_founder_polyploid(N, m, A, p)

        for i in range(len(expctd)):
            for j in range(len(expctd[i])):
                self.assertAlmostEqual(expctd[i][j],real[i][j],delta=THRESHOLD)

    def test_wright_fisher_polyploid_sanity(self):
        """
        Second unit test for polyploid implementation. In this test, the full
        Rosalind default dataset is passed and ploidy set to 2. The polyploid
        function should work on any ploidy including 2, and therefore should
        return a result that is equivalent to the original unextended algorithm.
        """
        N = 4
        m = 3
        A = [0, 1, 2]

        expctd = wright_fisher_founder(N, m, A)
        real = wright_fisher_founder_polyploid(N, m, A, 2)

        for i in range(len(expctd)):
            for j in range(len(expctd[i])):
                self.assertAlmostEqual(expctd[i][j],real[i][j],delta=THRESHOLD)