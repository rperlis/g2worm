import unittest
from calculator_worm import *

class DefaultVariablesTest(unittest.TestCase): 

	def test_defaults(self): 

		risk = make_pathways(100)

		print(risk)

if __name__ == '__main__':
	unittest.main()
