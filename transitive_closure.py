#!/usr/bin/env python

from copy import deepcopy
from api.graph import Graph, read_data, write_data, make_symmetric, clean_graph

def find_paths(G):
	"""
	Extracts all paths between two nodes in a graph if these are connected by single
	intermediate node.
	"""

	temp_dic = dict()

	for initial_node in sorted(G.graph):
		target_nodes = set(G.graph[initial_node].keys())

		for target_node in target_nodes:

			transitive_nodes = target_nodes.intersection(set(G.graph[target_node].keys()))

			if initial_node in transitive_nodes:
				raise SelfEdgeError("Graph seems to contain self edges, please check the input!.")

			if (target_node, initial_node) not in temp_dic and len(transitive_nodes) > 0:
				temp_dic[(initial_node, target_node)] = transitive_nodes

	paths = set()

	for path in temp_dic:
		initial_node = path[0]
		terminal_node = path[1]

		for intermediate_node in temp_dic[path]:
			paths.add((initial_node, intermediate_node, terminal_node))

	return paths

# Seems to be buggy

# def find_paths(G):
# 	"""Extracts all paths between two nodes in a graph if these are connected by single
# 	intermediate node."""

# 	paths = set()

# 	for initial_node in G.graph:
# 		target_nodes = set(G.graph[initial_node].keys())

# 		for intermediate_node in G.graph:
# 			if initial_node == intermediate_node: continue
	
# 			terminal_nodes = set(G.graph[intermediate_node])
			
# 			for target_node in target_nodes:
# 				if target_node in terminal_nodes:
# 					if (target_node, intermediate_node, initial_node) not in paths:
# 						paths.add((initial_node, intermediate_node, target_node))

# 	return paths

def linear_transitive_closure(G, paths, verbose):

	P = deepcopy(G)
	updated = False
	calls = 0

	for path in paths:
		calls += 1

		initial_node, other_node, terminal_node = path

		forward_prob = G.get_edge(initial_node, terminal_node)
		backward_prob = G.get_edge(terminal_node, initial_node)

		if forward_prob != backward_prob:
			raise ValueError("Graph does not seem to be undirected!")

		new_prob = G.get_edge(initial_node, other_node) * G.get_edge(other_node, terminal_node) 

		if new_prob > forward_prob and new_prob > backward_prob:
			if verbose: 
				print ("Updating egde P(%s,%s) = P(%s,%s) = %1.2f with P(%s,%s) * P(%s,%s) = %1.2f * %1.2f = %1.2f" % (initial_node, 
				terminal_node, terminal_node, initial_node, backward_prob, initial_node, other_node, other_node, 
				terminal_node, G.get_edge(initial_node, other_node), G.get_edge(other_node, terminal_node), new_prob))

			P.add_edge(initial_node, terminal_node, new_prob)
			P.add_edge(terminal_node, initial_node, new_prob)
			updated = True

	if verbose: 
		print ("Calls of transitive_closure %s." % calls)

	return updated, P

def linear_iterative_transitive_closure(G, verbose):
	"""Repeat transitive_closure until no changes in the graph occur!"""
	
	P = deepcopy(G) 
	updated = True
	paths = find_paths(G)
	counter = 1

	
	while updated:
		
		if verbose: 
                     print("\nIteration %d \n" % counter)
		updated, P = linear_transitive_closure(P, paths, verbose)
		counter += 1
		
	return P


def opt():
	import optparse
	description = """Performs a probablistic transitive closure on a similarity graph. May filter the input graph by coverage, 
	or edge weights before or after transitve closure. Consider that the probablistic transitve closure may 
	increase the weight of a edge connecting two nodes, thus the prefiltering weights may yield different results as compared to 
	postfiltering. 
	"""


    # Initiate a OptionParser Class
	parser = optparse.OptionParser(description=description)

    # Call add_options to the parser
	parser.add_option('-i', help='Input data.', dest='input')
	parser.add_option('-o', help='Output data.', dest='output')
	parser.add_option('-c', help='Coverage cut-off (default: no cutoff = 0).', dest='cut_off', type='float', default=0)
	parser.add_option('-v', help='Verbose mode (default: False).', dest='verbose', action="store_true", default=False)
	
	parser.add_option('-a', help='Pre filter cut-off (filters edges below a cut-off before transtive closure, default: no cutoff = 0).', dest='pre', type='float', default=0)
	parser.add_option('-b', help='Post filter cut-off (filters edges below a cut-off after transitive closure, default: no cutoff = 0).', dest='post', type='float', default=0)

	return parser


def main():
	parser = opt()
	(options, argv) = parser.parse_args()

	if options.verbose: 
		print ("Reading data ...\n")
	
	data = read_data(options.input, options.cut_off)

	if options.verbose: 
		print ("Making graph symmetric and performing transitive closure ...")
	
	symmetric_graph = make_symmetric(data)

	if options.pre != 0:
		if options.verbose: 
			print ("Pre cut-off: Cleaning edges from graph with cutoff below %s ..." % options.pre) 
		
		symmetric_graph = clean_graph(symmetric_graph, options.pre, options.verbose)
	
	transitive_graph = linear_iterative_transitive_closure(symmetric_graph, options.verbose)

	if options.post != 0:
		if options.verbose: 
			print ("Post cut-off: Cleaning edges from graph with cutoff below %s ..." % options.post)
		
		cleaned_graph = clean_graph(transitive_graph, options.post, options.verbose)
		
		if options.verbose: 
			print ("Writing data ...")
		
		write_data(cleaned_graph, options.output)
	
	else:
		if options.verbose: 
			print ("Writing data ...")
		
		write_data(transitive_graph, options.output)
	

if __name__ == "__main__":
	main()

