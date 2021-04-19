"""
Unit test file for WFMD implementation.
"""

import unittest
from WFMD import WFMD

THRESHOLD = 0.005 # equivalent to a tolerance of 1/2 a percent

class WFMD_Tests(unittest.TestCase):
    def test_WFMD_Rosalind_default(self):
        """
        First unit test for WFMD, uses the example data that appears on the
        rosalind page for WFMD. The expected output is 0.772. The test succeeds
        if the WFMD algorithm spits out a probability between 0.772 - THRESHOLD
        and 0.772 + THRESHOLD, THRESHOLD being some very small variance
        tolerance.
        """
        N = 4
        m = 6
        g = 2
        k = 1

        expected = 0.772
        real = WFMD(N, m, g, k)

        self.assertAlmostEqual(expected, real, delta=THRESHOLD)

    def test_WFMD_hand_calculated(self):
        """
        Second unit test for WFMD, uses data that Bennet Montgomery calculated
        by hand on starting situation 2 individuals with an allele incidence of
        3. The outcome of the calculation was around but not exactly 0.302.
        """
        N = 2
        m = 3
        g = 3
        k = 2

        expected = 0.302
        real = WFMD(N, m, g, k)

        self.assertAlmostEqual(expected, real, delta=THRESHOLD)

    def test_WFMD_zero_output(self):
        """
        Third unit test for WFMD, using an outlier case where no generations
        are calculated. The probabiltiy of the allele having the same incidence
        or greater after 0 generations should be 1.
        """
        N = 1
        m = 0
        g = 0
        k = 0

        self.assertEqual(WFMD(N,m,g,k), 1)