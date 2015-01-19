"""
Tests the properties of the hydrogen gas
"""

import unittest

from .. import h2_mu

class H2Test(unittest.TestCase):

    """The test case for H2"""

    def test_mu(self):
	"""Tests the chemical potential of H2 at normal condition"""

	mu = h2_mu(300.0, 101325.0)
	self.assertAlmostEqual(mu, -0.301376045993850)

