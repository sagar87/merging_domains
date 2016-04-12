#!/usr/bin/env python

from copy import deepcopy

class Graph():

    def __init__(self):
        self.graph = {}

        self.nodes = 0
        self.edges = 0

    def __str__(self):
        """Displays the graph as a table."""

        header = "\t"
        table = ""

        for terminal_node in self.graph:
            header += "{:.4}\t".format(terminal_node)

        header += "\n" + "--------" * (len(self.graph)+ 1) + "\n"

        for initial_node in self.graph:
            row = "{:.4}\t".format(initial_node)

            for terminal_node in self.graph.keys():
                if terminal_node in self.graph[initial_node]:
                    row += "%1.2f\t" % self.graph[initial_node][terminal_node]
                else:
                    row += "\t"

            table += row + "\n"

        return header + table

    def __iter__(self):
        for node in self.graph.keys():
            yield node

    def __getitem__(self, node):
        if node not in self.graph.keys():
            raise IndexError("Accessed node is not present in the graph.")

        return self.graph[node]

    def get_edge(self, initial_node, terminal_node):
        """Returns the weight of the edge pointing from initial_node to 
        terminal_node if it exists, 0 otherwise"""

        try:
            return self.graph[initial_node][terminal_node]
        except KeyError:
            return 0

    def total_nodes(self):
        """Returns the number of nodes in the graph."""
        return self.nodes

    def total_edges(self):
        return self.edges

    def add_node(self, node):
        if node not in self.graph:
            self.graph[node] = dict()
            self.nodes += 1

    def add_edge(self, initial_node, terminal_node, weight):
        """Adds or modifies a edge to the graph."""
        
        if initial_node not in self.graph:
            raise ValueError("Initial Node is not present in similarity graph.")
        if terminal_node not in self.graph:
            raise ValueError("Terminal node is not present in similarity graph.")

        if terminal_node not in self.graph[initial_node]:
            self.edges += 1
        
        self.graph[initial_node][terminal_node] = weight
     
    def delete_edge(self, initial_node, terminal_node):
        """Deletes a edge."""

        del self.graph[initial_node][terminal_node]
        self.edges -= 1

def read_data(in_file, cut_off):
    """Reads data into a graph. Does not allow self edges!"""

    G = Graph()

    with open(in_file) as data:
        for num, line in enumerate(data):
            if num%1000 == 0:
                print("Reading line {0}...      ".format(num), end="")
            # There is always an intial node
            initial_node = line.split()[0]
            G.add_node(initial_node)

            try:
                terminal_node = line.split()[1]
            except IndexError:
                terminal_node = None

            if terminal_node != None:
                G.add_node(terminal_node)

            try:
                weight = float(line.split()[2])
            except IndexError:
                weight = None
            
            try:
                cov = float(line.split()[3])
            except IndexError:
                cov = None # some old datasets do not have coverage data -> we use all data

            if initial_node != terminal_node and weight != None: # no self edges
                if cov == None: # only if coverage satisfies criteria
                    G.add_edge(initial_node, terminal_node, weight) 
                elif cov > cut_off:
                    G.add_edge(initial_node, terminal_node, weight) 
            
    return G

def write_data(G, out_file):

    with open(out_file, 'w+') as f_out:

        for initial_node in G.graph.keys():
            if len(G.graph[initial_node]) == 0:
                f_out.write("%s\n" % initial_node)
            else:
                for terminal_node in G.graph[initial_node]:
                    f_out.write("%s\t%s\t%s\n" % (initial_node, terminal_node, G.graph[initial_node][terminal_node]))


def make_symmetric(G):
    """Converts graph G to a symmetric, undirected graph S."""

    S = deepcopy(G)

    for initial_node in G:
        for terminal_node in G[initial_node]:
            # select the higher probability
            max_prob = max([G.get_edge(initial_node, terminal_node), G.get_edge(terminal_node, initial_node)])

            S.add_edge(initial_node, terminal_node, max_prob)
            S.add_edge(terminal_node, initial_node, max_prob)
    
    return S

def clean_graph(G, cut_off, verbose):
    """Cleans the graph from edges that are below"""
    C = deepcopy(G)
    
    for node in G:
        for edge in G[node]:
            if G[node][edge] < cut_off:
                
                if verbose:
                    print ("Deleting %s -> %s (%s)." % (node, edge, G[node][edge]))
                C.delete_edge(node, edge)

    return C

def arg():
    import argparse
    description = """Prints number of nodes and edges to stdout."""
    epilog= '"It is not the strongest of the species that survive, nor the most intelligent, but the one most responsive to change." -Charles Darwin'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help='Input graph') #  nargs = '+' for more one argument
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode', default=False)
    parser.add_argument('-c', type=float, help='Coverage cut-off (default = 0)', default=0)
    parser.add_argument('-a', type=float, help='Pre-filter (default = 1)', default=1)


    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    similarity_graph = read_data(args.input, args.c)
    
    print ("Total nodes: %s." % similarity_graph.total_nodes())
    print ("Total edges: %s." % similarity_graph.total_edges())

    if args.a != 1:
        similarity_graph = clean_graph(similarity_graph, args.a, args.verbose)
        print ("Total after filtering edges:", similarity_graph.total_edges())
    
    symmetric = make_symmetric(similarity_graph)
    print ("Edges symmetric graph:", symmetric.total_edges())


if __name__ == "__main__":
    main()
