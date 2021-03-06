helper/mapping.py                                                                                   000755  000765  000024  00000002645 12654636200 014560  0                                                                                                    ustar 00Sagar                           staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python
from collections import defaultdict

def map_identifier(input_file):

	key_to_idx = defaultdict()
	# graph = defaultdict(list)
	idx = 0
	line_counter = 0
	with open(input_file) as in_fh:

		for line in in_fh:
			line_counter += 1

			if line_counter % 1000 == 0:
				print ("Processing Line {}...".format(line_counter))
			
			query, template, prob = line.split("\t")

			if query not in key_to_idx:
				key_to_idx[query] = idx 
				idx += 1

			if template not in key_to_idx:
				key_to_idx[template] = idx
				idx += 1

	return key_to_idx

def write_mapping(key_to_idx, out_file):

	with open(out_file, 'w') as out_fh:

		for k, v in key_to_idx.items():
				out_fh.write("{str}\t{idx}\n".format(str=k, idx=v))


def arg():
    import argparse
    description = """Converts huge graphs in hopefully smaller ones."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help='Input file') #  nargs = '+' for more one argument
    parser.add_argument('map', help='Output mapping') #  nargs = '+' for more one argument

    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    key_to_idx = map_identifier(args.input)
    write_mapping(key_to_idx, args.map)


if __name__ == "__main__":
    main()                                                                                           helper/quick_find.py                                                                                000755  000765  000024  00000004223 12654603510 015231  0                                                                                                    ustar 00Sagar                           staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python

class QuickFind():

    def __init__(self, N):
        self._count = N
        self._idx = [ i for i in range(N) ]
        self._size = [ 1 for dummy_size in range(N) ]

    def count(self):
        return self._count

    def connected(self, p, q):
        return self.find(p) == self.find(q)

    def find(self, p):
        while p != self._idx[p]:
            # this is path compression
            self._idx[p] = self._idx[self._idx[p]]
            p = self._idx[p]

        return p

    def union(self, p, q):
        i = self.find(p)
        j = self.find(q)
        
        if (i == j): return

        # weightening -> add smaller root to the larger one
        if (self._size[i] < self._size[j]):
            self._idx[i] = j
            self._size[j] += self._size[i]  
        else:
            self._idx[j] = i
            self._size[i] += self._size[j]

        self._count -= 1

def arg():
    import argparse
    description = """A connected components algorithm based on weighted quick union (with path compression)."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help='Input graph') #  nargs = '+' for more one argument
    parser.add_argument('output', help='Connected components') #  nargs = '+' for more one argument
    parser.add_argument('N', type=int, help='Size of the QuickFind Array.')

    return parser

def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    QF = QuickFind(args.N+1)

    print ("Created Quick Find Array")

    with open(args.input) as in_fh:
        for i, line in enumerate(in_fh):
            if i % 10000 == 0:
                print ("Processing Edge {0}.".format(i))
            p, q = map(int, line.split()[:2])
            QF.union(p, q)

    with open(args.output, 'w') as out_fh:
        out_fh.write(str(QF._count) + "\n")
        for node in range(args.N+1):
            line = str(node) + "\t" + str(QF.find(node)) + "\n"
            out_fh.write(line)

if __name__ == "__main__":
    main()                                                                                                                                                                                                                                                                                                                                                                             helper/reformat.py                                                                                  000755  000765  000024  00000002665 12654603403 014745  0                                                                                                    ustar 00Sagar                           staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python

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
    main()                                                                           helper/write_graph.py                                                                               000755  000765  000024  00000005135 12654634362 015444  0                                                                                                    ustar 00Sagar                           staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python

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
    main()                                                                                                                                                                                                                                                                                                                                                                                                                                   cc_paths.py                                                                                         000644  000765  000024  00000005205 12654636470 013433  0                                                                                                    ustar 00Sagar                           staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python

import os, subprocess, tempfile, shutil

# Enter the path to the subscripts
MAPPING = ''
REFORMAT = ''
QUICKFIND = ''
WRITEGRAPH = ''


def wrapper(in_file, out_folder):

    # set up a tempfolder
    tmp_dir = tempfile.mkdtemp()
    tmp = os.path.join(tmp_dir, tmp_dir.split("/")[-1])
    # print tmp
    # print os.path.join(os.getcwd(), "helper/mapping.py")
    
    # creates a file for the mapping
    tmp_map = tmp + ".map"
    tmp_grp = tmp + ".grp"
    tmp_cc = tmp + ".cc"

    print ("Mapping identifier to integers...")

    mapping = subprocess.Popen(" ".join(['python3', MAPPING, in_file, tmp_map]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = mapping.communicate()

    print ("Identifying total number of nodes...")
    # Get the number of total nodes
    nodes = subprocess.Popen(" ".join(["wc", "-l", tmp_map]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = nodes.communicate()

    num_nodes = int(output.split()[0]) - 1

    print ("Total number of nodes {0}. Reformatting graph...".format(num_nodes))

    reformat = subprocess.Popen(" ".join(['python3', REFORMAT, tmp_map, in_file, tmp_grp]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = reformat.communicate()

    print ("Identifying connected components.")

    cc = subprocess.Popen(" ".join(['python3', QUICKFIND, tmp_grp, tmp_cc, str(num_nodes)]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = cc.communicate()    

    print ("Writing connected components...")

    wc = subprocess.Popen(" ".join(['python3', WRITEGRAPH, tmp_map, tmp_cc, tmp_grp, out_folder]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = wc.communicate() 

    print ("All done!")

    shutil.rmtree(tmp_dir)

def arg():
    import argparse
    description = """A pipeline for connected components."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help='Input graph') #  nargs = '+' for more one argument
    parser.add_argument('output', help='Output Connected components folder') #  nargs = '+' for more one argument


    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    wrapper(args.input, args.output)

if __name__ == "__main__":
    main()                                                                                                                                                                                                                                                                                                                                                                                           cc.py                                                                                               000755  000765  000024  00000005233 12654635410 012231  0                                                                                                    ustar 00Sagar                           staff                           000000  000000                                                                                                                                                                         #!/usr/bin/env python

import os, subprocess, tempfile, shutil




def wrapper(in_file, out_folder):

    # set up a tempfolder
    tmp_dir = tempfile.mkdtemp()
    tmp = os.path.join(tmp_dir, tmp_dir.split("/")[-1])
    # print tmp
    # print os.path.join(os.getcwd(), "helper/mapping.py")
    
    # creates a file for the mapping
    tmp_map = tmp + ".map"
    tmp_grp = tmp + ".grp"
    tmp_cc = tmp + ".cc"

    print ("Mapping identifier to integers...")

    mapping = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "helper/mapping.py"), in_file, tmp_map]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = mapping.communicate()

    print ("Identifying total number of nodes...")
    # Get the number of total nodes
    nodes = subprocess.Popen(" ".join(["wc", "-l", tmp_map]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = nodes.communicate()

    num_nodes = int(output.split()[0]) - 1

    print ("Total number of nodes {0}. Reformatting graph...".format(num_nodes))

    reformat = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "helper/reformat.py"), tmp_map, in_file, tmp_grp]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = reformat.communicate()

    print ("Identifying connected components.")

    cc = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "helper/quick_find.py"), tmp_grp, tmp_cc, str(num_nodes)]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = cc.communicate()    

    print ("Writing connected components...")

    wc = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "helper/write_graph.py"), tmp_map, tmp_cc, tmp_grp, out_folder]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = wc.communicate() 

    print ("All done!")

    shutil.rmtree(tmp_dir)

def arg():
    import argparse
    description = """A pipeline for connected components."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('input', help='Input graph') #  nargs = '+' for more one argument
    parser.add_argument('output', help='Output Connected components folder') #  nargs = '+' for more one argument


    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    wrapper(args.input, args.output)

if __name__ == "__main__":
    main()                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     