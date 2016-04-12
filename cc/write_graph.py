#!/usr/bin/env python

from collections import defaultdict
import os

def read_components(in_file):

    components = defaultdict()

    with open(in_file) as in_fh:
        for i, line in enumerate(in_fh):
            if i == 0:
                total_comp = int(line.split()[0])
            else:
                node, component = map(int, line.split())
                components[node] = component

    return total_comp, components

def read_mapping(in_file):

    mapping = defaultdict()
    
    with open(in_file) as in_fh:
        for line in in_fh:
            identifier = line.split()[0]
            int_idx = int(line.split()[1])

            mapping[int_idx] = identifier

    return mapping

def rewrite_graph(in_file, out_file, mapping, component):
    """
    Writes each component into a separate file.
    """

    graph = open(in_file)
    out_fh = open(out_file, 'w')

    for line in graph:
        initial_node, terminal_node, prob = map(int, line.split())
        prob = prob / float(10000)
        
        # verify that initial node terminal node are found in the same component
        if component[initial_node] != component[terminal_node]:
            print (initial_node, component[initial_node])
            print (terminal_node, component[terminal_node])
            raise ValueError("Something went wrong whilst identifying connected components.")

        cluster = component[initial_node]

        out_str = "{c}\t{i}\t{o}\t{p}\n".format(c=cluster, i=mapping[initial_node], 
            o=mapping[terminal_node], p=prob)
        out_fh.write(out_str)

    graph.close()
    out_fh.close()


def arg():
    import argparse
    description = """Takes a mapping, a graph, its componets and writes all sub graphs into seperate files."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('mapping', help='Input mapping') #  nargs = '+' for more one argument
    parser.add_argument('components', help='Input components') #  nargs = '+' for more one argument
    parser.add_argument('graph', help='Input graph') #  nargs = '+' for more one argument
    parser.add_argument('output', help='Output folder') #  nargs = '+' for more one argument

    return parser

def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    mapping = read_mapping(args.mapping)
    total_nodes, components = read_components(args.components)
    rewrite_graph(args.graph, args.output, mapping, components)



if __name__ == "__main__":
    main()