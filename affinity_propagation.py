#!/usr/bin/env python

import numpy as np
import numpy.matlib
from itertools import groupby
import os, subprocess, tempfile, shutil
from collections import defaultdict

def read_file(data, input_preference):
    """ 
    Reads in a graph.
    """
    
    sub_graphs = list()

    with open(data) as fin:
        
        for component, sub_graph in groupby(fin, lambda line: line.split()[0]):
            
            idx = 0
            node_idx = {} # maps each node to its index in the matrix
            edges = {} # maps node i,j tuples to the corresponding weights

            for edge in sub_graph:
                component, initial_node, terminal_node, prob = edge.strip().split()

                if initial_node not in node_idx:
                    node_idx[initial_node] = idx
                    idx += 1

                if terminal_node not in node_idx:
                    node_idx[terminal_node] = idx
                    idx += 1    

                edges[(node_idx[initial_node], node_idx[terminal_node])] = float(prob)

            N = max(node_idx.values())+1

            # print ("Found %s nodes." % N)

            S = np.zeros([N, N])
            
            for i, j in edges.keys():
                S[i, j] = edges[(i, j)]

            
            if input_preference == 0:
                M = np.median(S)
                # print ("Setting input preference (Median) to %s." % M)
                np.fill_diagonal(S, M)
            else:
                # print ("Setting input preference to %s." % input_preference)
                np.fill_diagonal(S, input_preference)

            sub_graphs.append((S, node_idx))
    
    return sub_graphs

def affinity_propagation(S, threshold, lam):
    
    dim = S.shape

    if dim[0] != dim[1]:
        raise ValueError("S has to be a square matrix !")

    m, n = dim 
    # Initialize matrices
    N = np.zeros(dim)
    A = np.zeros(dim)
    R = np.zeros(dim)

    S = S + 1e-12 * np.random.randn(m, n) * (np.max(S) - np.min(S)) # Remove Degeneracies
    
    exit = False
    it = 0

    while it < threshold:

        Rold = np.copy(R)
        AS = A + S

        Y = np.max(AS, axis=1) # an array with the maximum value along each column
        Y = np.array(Y)[np.newaxis]
        I = np.argmax(AS, axis=1) # an array which contains all row indices of the maximum values are found

        for i in range(m): 
            AS[i, I[i]] = - float('inf') # sets the highest value of matrix AS to - infinity

        Y2 = np.max(AS, axis=1) 
        I2 = np.argmax(AS, axis=1) 

        R = S - np.matlib.repmat(Y.T, 1, n) # subtracts from the similarty matrix the highest values of each row

        for i in range(m):
            R[i, I[i]] = S[i, I[i]] - Y2[i]

        R = (1-lam) * R + lam * Rold

        # compute in A
        
        Aold = np.copy(A)
        Rp = np.maximum(R, 0) # replaces all values in R that are < 0 with 0

        for k in range(m):
            Rp[k, k] = R[k, k]

        A = np.matlib.repmat(np.sum(Rp, axis = 0), m, 1) - Rp
        dA = np.diag(A)
        A = np.minimum(A, 0)

        for k in range(m):
            A[k, k] = dA[k]

        A = (1-lam) * A + lam * Aold

        # Check if values of matrices have changed
        Adiff = A - Aold
        Rdiff = R - Rold

        Amax = Adiff.max()
        Rmax = Rdiff.max()

        # print ("Iteration %s (Rmax = %e, Amax = %e, max Iteration = %s)." % (it, Amax, Rmax, threshold))

        it += 1

    E = R + A
    
    I = np.where(np.diag(E)>0)
    I = I[0]
    
    K = max(I.shape)

    tmp = np.max(S[:, I], 1)
    c = np.argmax(S[:, I], 1)

    c[I] = np.array(range(K))
    idx = I[c]

    return idx      

def read_predictions(in_file):

    predictions = set()

    with open(in_file) as in_fh:
        for line in in_fh:
            domain = line.strip().split()[1]
            predictions.add(domain)

    return predictions

def read_cluster(in_file):

    domains = set()
    clusters = defaultdict(list)

    with open(in_file) as in_fh:
        for line in in_fh:
            exemplar, candidate = line.strip().split()
            domains.add(exemplar)
            domains.add(candidate)
            clusters[exemplar].append(candidate)

    return domains, clusters

def write_to_file(idx, node_idx, mode, outfile):

    idx_node = { value:key for key, value in node_idx.items() }
    clusters = { cluster:list() for cluster in np.unique(idx) }

    for node, cluster in enumerate(idx):
        clusters[cluster].append(node)

    with open(outfile, mode) as fout:
        for cluster in clusters:
            for node in clusters[cluster]:
                string = idx_node[cluster] + "\t" + idx_node[node] + "\n"
                fout.write(string)

def opt():
    import optparse
    usage = "usage: affinity_propagation.py -n number of nodes (int) -p edge probability (float) -i number of iterations (int) "
    description = """Performs affinity propagation on a random graph."""
    # Initiate a OptionParser Class
    parser = optparse.OptionParser(usage=usage, description=description)
    
    # Call add_options to the parser
    parser.add_option('-i', help='Input data',  dest='data')
    parser.add_option('-p', help='Predictions', dest='pred')
    
    parser.add_option('-s', help='Self responsibility is a float between 0 and 1 (default = median of the similarty matrix)', 
        default=0, type='float', dest='preference')
    parser.add_option('-d', help='Dampening factor (default = 0.5)', 
        default=0.5, type='float', dest='lam')
    parser.add_option('-t', help='Number of iterations (default = 500)', 
        default=500, type='int', dest='threshold')
    parser.add_option("-o", help="Output file.", dest="output") 

    return parser


def main():
    parser = opt()
    (options, argv) = parser.parse_args()

    tmp_dir = tempfile.mkdtemp()
    tmp = os.path.join(tmp_dir, tmp_dir.split("/")[-1])
    tmp_clu = tmp + ".clu"

    subgraphs = read_file(options.data, options.preference)
    print ("Total number of subgraphs %s." % len(subgraphs))

    for i, subgraph in enumerate(subgraphs):
        S, node_idx = subgraph
        idx = affinity_propagation(S, options.threshold, options.lam)

        mode = 'w' if i == 0 else 'a'
        
        write_to_file(idx, node_idx, mode, tmp_clu)

    domains, clusters = read_cluster(tmp_clu)
    predictions = read_predictions(options.pred)

    missing_fragments = predictions-domains

    with open(options.output, 'w') as out_fh:
        print ("Writing graph ...")
        for cluster in clusters:
            for fragment in clusters[cluster]:
                out_fh.write('{cluster}\t{frag}\n'.format(cluster = cluster, frag = fragment))

        print ("Writing Singletons ...")
        for fragment in missing_fragments:
            out_fh.write('{frag}\t{frag}\n'.format(frag = fragment))

    print ("All done.")

    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
    main()