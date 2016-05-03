try:
	import graph_tool
except ImportError:
	self.fail("Graph-tool is not or incorrectly installed.")