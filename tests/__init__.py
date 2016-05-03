#!/usr/bin/env python
import unittest

class InitializationTests(unittest.TestCase):

	def test_initialization(self):
		#Check the test suite runs by affirming 1+1=2
		self.assertEqual(1+1,2)
	def test_import(self):
		#Ensure the test suite can import the module
		try: 
 			import PDS
		except ImportError: 
 			self.fail("PDS could not be imported.")

 	def test_graph_tool(self):
 		# checks if graph-tool is installed
 		try:
 			import graph_tool
 		except ImportError:
 			self.fail("Graph-tool is not or incorrectly installed.")