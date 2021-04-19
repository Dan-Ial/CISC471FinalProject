import unittest
from WFMD import WFMD

THRESHOLD = 0.005 # equivalent to a tolerance of 1/2 a percent

class WFMD_Tests(unittest.TestCase):
    def test_WFMD_Rosalind_default(self):
        N = 4
        m = 6
        g = 2
        k = 1

        expected = 0.772
        real = WFMD(N, m, g, k)

        self.assertAlmostEqual(expected, real, delta=THRESHOLD)

    def test_WFMD_hand_calculated(self):
        N = 2
        m = 3
        g = 3
        k = 2

        expected = 0.302
        real = WFMD(N, m, g, k)

        self.assertAlmostEqual(expected, real, delta=THRESHOLD)

    def test_WFMD_zero_output(self):
        N = 1
        m = 0
        g = 0
        k = 0

        self.assertEqual(WFMD(N,m,g,k), 1)