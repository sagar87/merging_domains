#!/usr/bin/env python

from collections import defaultdict

def rewrite_graph(mapping, in_file, out_file):

	in_fh = open(in_file)
	out_fh = open(out_file, 'w')

	for line in in_fh:
		initial_node, terminal_node, prob = line.split("\t")
		int_prob = float(prob) * 10000
		int_prob = int(int_prob)
		
		w_line = "{in_node}\t{out_node}\t{num}\n".format(in_node = mapping[initial_node], 
			out_node = mapping[terminal_node], num=int_prob)
		out_fh.write(w_line)

	in_fh.close()
	out_fh.close()

def read_mapping(in_file):
	
	mapping = defaultdict()
	
	with open(in_file) as fh:
		for line in fh:
			idx, idx_int = line.split("\t")

			mapping[idx] = int(idx_int)

	return mapping

def arg():
    import argparse
    description = """A description"""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('mapping', help='mapping file') #  nargs = '+' for more one argument
    parser.add_argument('graph', help='graph file') #  nargs = '+' for more one argument
    parser.add_argument('output', help='output graph file') #  nargs = '+' for more one argument

    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    mapping = read_mapping(args.mapping)
    rewrite_graph(mapping, args.graph, args.output)

if __name__ == "__main__":
    main()