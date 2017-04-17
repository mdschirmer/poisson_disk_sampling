import sys

try:
	import graph_tool
except ImportError:
	print("Graph-tool is not or incorrectly installed.")
	sys.exit()
