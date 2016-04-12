#!/usr/bin/env python

from alignment_matrices import *
from ffindex import read_data, read_index, read_entry_data
from collections import defaultdict
from math import ceil
from multiprocessing import Pool
import os, subprocess, tempfile, shutil

DEBUG = False

class AlignmentMatrix():

    def __init__(self, query, temp, prob, pred, fmat, bmat, post):
        """
        Lousy name for a fancy object.
        """
        self._query = query
        self._temp = temp
        self._prob = prob

        self._pred = pred
        self._fmat = fmat # holds the forward matrix
        self._bmat = bmat # holds the backward matrix
        self._post = self._convert_matrix(post)

    def __repr__(self):
        return "{query} -> {template}".format(query=self._query, template=self._temp)

    def _convert_matrix(self, sparse_matrix):

        max_query = max(map(lambda x: x[0], sparse_matrix))
        max_temp = max(map(lambda x: x[1], sparse_matrix))

        mat = np.zeros([max_temp, max_query])

        for entry in sparse_matrix:
            mat[entry[1]-1, entry[0]-1] = entry[2]

        return mat

    def display_post(self):
        """
        Displays the posterior matrix in a convinient way :).
        """

        nrow, ncol = self._post.shape
        
        ans = ""
        
        for i in range(nrow):
            for j in range(ncol):                
                if self._post[i, j] == 0:
                    formated = '0'
                else:
                    formated = "{0:.1f}".format(self._post[i, j])
                
                while len(formated) < 4:
                    formated += ' '
                    
                ans += formated                 
            ans += "\n"
        
        return ans

    def qstart(self, start, end, n, offset):
        """
        Computes the probability qstart for query fragment n.
        """
        
        qstart = 0

        for entry in sorted(self._bmat, key=lambda x: x[0]):
            qidx, tidx, prob = entry 
            
            if qidx > start + offset:
                break

            qstart += prob
        
        if DEBUG:
            print ('|-P(query fragment {n} starts <= {pend}) = {qs:.4f} (dm = {dm}, {s} -> {e} ({l})).'.format(n = n,
                s = start, e = end, l = end - start, dm = offset, pend = start + offset, qs = qstart))
        
        return qstart

    def qend(self, start, end, n, offset):
        """
        Computes the probability qend for a query fragment n.
        """

        qend = 0

        for entry in sorted(self._fmat, key=lambda x: x[0]):
            qidx, tidx, prob = entry 
                    
            if qidx >= end - offset:
                qend += prob

        if DEBUG:
            print ('|-P(query fragment {n} ends >= {pend}) = {qe:.4f} (dm = {dm}, {s} -> {e} ({l})).'.format(n = n,
                s = start, e = end, l = end - start, dm = offset, pend = end - offset, qe = qend))

        return qend

    def tstart(self, start, end, n, offset):
        """
        Computes the probability tstart for template fragment n.
        """
        
        tstart = 0

        for entry in sorted(self._bmat, key=lambda x: x[1]):
            qidx, tidx, prob = entry 
            
            if tidx > start + offset:
                break
            
            if tidx <= start + offset:
                tstart += prob
        
        if DEBUG:
            print ('  |-P(template fragment {n} starts <= {pend}) = {qe:.4f} (dm = {dm}, {s} -> {e} ({l})).'.format(n = n,
                s = start, e = end, l = end - start, dm = offset, pend = start + offset, qe = tstart))
        
        return tstart

    def tend(self, start, end, n, offset):
        """
        Computes the probability tend for a template fragment n.
        """

        tend = 0

        for entry in sorted(self._fmat, key=lambda x: x[1]):
            qidx, tidx, prob = entry 
            
            if tidx >= end - offset:
                tend += prob
        
        if DEBUG:
            print ('  |-P(template fragment {n} ends >= {pend}) = {qe:.4f} (dm = {dm}, {s} -> {e} ({l})).'.format(n = n,
                s = start, e = end, l = end - start, dm = offset, pend = end - offset, qe = tend))

        return tend

    def ali_starts_after(self, qstart, tstart, offset):

        if offset < 1: 
            raise ValueError('Offset has to be greater or equal to 1.')

        # extract matrix boundaries and set the col/row offset to a max value
        nrow, ncol = self._post.shape
        col_offset, row_offset = sys.maxsize, sys.maxsize

        # correct indexing
        qstart = qstart - 1
        tstart = tstart - 1

        # check whether the start positions make sense
        if qstart >= ncol or tstart >= nrow:
            if DEBUG:
                print ('  |- Warning ! Check domain boundaries qstart = {qs} and template start {ts}.'.format(
                    qs = qstart, ts = tstart, dm = offset, re = max(terms)))
            return 0

        # correct for matrix dimensions
        if qstart + offset > ncol:
            col_offset = offset + (ncol - (qstart + offset)) 
        if tstart + offset > nrow:
            row_offset = offset + (nrow - (tstart + offset))

        if col_offset != sys.maxsize or row_offset != sys.maxsize:
            offset = min(row_offset, col_offset)
            
            if DEBUG:
                print ('Corrected offset {0}'.format(offset))

        tmp_mat = self._post[tstart : tstart + offset, qstart : qstart + offset]
        
        terms = list()

        for i in range(tmp_mat.shape[0]):
            term = np.sum(tmp_mat[i, 0 : i]) + np.sum(tmp_mat[0 : i, i]) + tmp_mat[i, i]
            terms.append(term)

        if DEBUG:
            print ('  |-P(ali starts after query = {qs}, template = {ts}) = {re:.4f} (dm = {dm}).'.format(
                qs = qstart, ts = tstart, dm = offset, re = max(terms)))
        if max(terms) > 1:
            import pdb; pdb.set_trace()

        return max(terms)

    def ali_starts_before(self, qstart, tstart, offset):

        if offset < 1: 
            raise ValueError('Offset has to be greater or equal to 1.')

        # extract matrix boundaries and set the col/row offset to a max value
        nrow, ncol = self._post.shape
        col_offset, row_offset = offset, offset

        # correct for matrix dimensions
        if qstart - offset < 0:
            col_offset += qstart - offset 
        if tstart - offset < 0:
            row_offset += tstart - offset

        if col_offset != offset or row_offset != offset:
            offset = min(col_offset, row_offset)
            
            if DEBUG:
                print ('Corrected offset {0}'.format(offset))

        # extract the corresponding submatrix
        tmp_mat = self._post[tstart - offset : tstart, qstart - offset : qstart]
        terms = list()

        for i in range(tmp_mat.shape[0]):
            term = np.sum(tmp_mat[i, i+1:]) + np.sum(tmp_mat[i+1:, i]) + tmp_mat[i, i]
            terms.append(term)

        return max(terms)

    def ali_ends_before(self, qend, tend, offset):

        if offset < 1: 
            raise ValueError("Offset has to be greater or equal to 1.")

        # extract matrix boundaries and set the col/row offset to a max value
        nrow, ncol = self._post.shape
        col_offset, row_offset = offset, offset

        # check whether the start positions make sense
        if qend <= 0 or tend <= 0:
            return 0

        if qend >= ncol or tend >= nrow:
            if DEBUG:
                print ('  |-Warning ! Query end = {qs} and template end = {ts} exceed dimensions of posterior matrix.'.format(
                    qs = qend, ts = tend))
            return 0

        # correct for matrix dimensions
        if qend - offset < 0:
            col_offset += qend - offset 
        if tend - offset < 0:
            row_offset += tend - offset

        if col_offset != offset or row_offset != offset:
            offset = min(col_offset, row_offset)
            
            if DEBUG:
                print ('Corrected offset {0}'.format(offset))

        # extract the corresponding submatrix
        tmp_mat = self._post[tend - offset : tend, qend - offset : qend]
        terms = list()

        for i in range(tmp_mat.shape[0]):
            term = np.sum(tmp_mat[i, i+1:]) + np.sum(tmp_mat[i+1:, i]) + tmp_mat[i, i]
            terms.append(term)

        if DEBUG:
            print ('  |-P(ali ends before query = {qs}, template = {ts}) = {re:.4f} (dm = {dm}).'.format(
                qs = qend, ts = tend, dm = offset, re = max(terms)))
        if max(terms) > 1:
            import pdb; pdb.set_trace() 
        return max(terms)

    def ali_ends_after(self, qend, tend, offset):

        if offset < 1: 
            raise ValueError("Offset has to be greater or equal to 1.")

        # extract matrix boundaries and set the col/row offset to a max value
        nrow, ncol = self._post.shape
        col_offset, row_offset = sys.maxsize, sys.maxsize

        # correct indexing
        qend = qend - 1
        tend = tend - 1 

        # check whether the start positions make sense
        if qend >= ncol or tend >= nrow:
            raise ValueError("Domain boundaries do not match matrix dimensions.")

        # correct for matrix dimensions
        if qend + offset > ncol:
            col_offset = offset + (ncol - (qend + offset)) 
        if tend + offset > nrow:
            row_offset = offset + (nrow - (tend + offset))

        if col_offset != sys.maxsize or row_offset != sys.maxsize:
            offset = min(row_offset, col_offset)
            
            if DEBUG:
                print ("Corrected offset {0}".format(offset))

        tmp_mat = self._post[tend : tend + offset, qend : qend + offset]
        terms = list()

        for i in range(tmp_mat.shape[0]):
            term = np.sum(tmp_mat[i, 0 : i]) + np.sum(tmp_mat[0 : i, i]) + tmp_mat[i, i]
            terms.append(term)

        return max(terms)

    
def get_all_templates(query, data, prob):
    """ 
    Takes a query and returns all of the templates in a set that exceed the
    specified threshold prbability.
    """
    
    templates = set()
    
    for domain in sorted(query):

        lines = read_lines(query[domain], data)
        hhr = parse_result(lines)
        
        for tmp in hhr:
            if tmp.probability > prob:
                temp, temp_idx = tmp.template_id.split("_")
                templates.add(str(temp))
    
    return list(templates)
                
def read_predictions(in_file):
    """ 
    Reads in a file which contains predictions. Returns a dictionary which 
    maps each query sequence to its template sequence.
    """
    predictions = defaultdict(list)
    
    with open(in_file) as fh:
        for line in fh:
            strip_line = line.strip()
            domain = strip_line.split()[0]
            number = int(strip_line.split()[1].split("_")[1])

            positions = strip_line.split()[2].split(",")
            fragments = list()

            for pos in positions:
                start = int(pos.split("-")[0])
                # Added 6. Jan 2016 to handle cases where a fragements consist of only one residue
                try:
                    end  = int(pos.split("-")[1])
                except IndexError:
                    print("! Warning: This prediction seems to be corrupt {0}.".format(strip_line))
                    continue
                
                fragments.append((start, end, number))

            for frag in fragments:
                predictions[domain].append(frag)

    return predictions  

def wrapper(data):

    alignment = data[0]
    P_low = data[1]
    cov_min = data[2]

    if DEBUG:
        print ("\nProcessing Alignment: {a} with P(query, template) = {p}.\n".format(a=alignment, p=alignment._prob))

    # this are the results
    export = list()

    query = alignment._query
    templ = alignment._temp
    
    qseen = list()
    
    for qfrag in sorted(alignment._pred[query], key=lambda x: x[1]):
        # find the right boundaries 
        qbound = sorted([ p for p in alignment._pred[query] if p[2] == qfrag[2] ])
        
        # old qlen was the length of all query fragments taken together 
        qlen = sum([ p[1] - p[0] + 1 for p in qbound ])
        
        offset = min([ abs(ceil((qlen - 10)/float(2))), ceil(qlen * (1 - cov_min)) ])

        if offset < 1:
            offset = 1
        
        qs = qbound[0][0]
        qe = qbound[-1][1]
        qn = qfrag[2]

        if qn in qseen:
            continue
        else:
            qseen.append(qn)

        qend = alignment.qend(qs, qe, qn, offset)          
        prob_mn = alignment._prob * qend

        if prob_mn < P_low:
            if DEBUG:
                print ('|_P(q, t) * qend = {mn:.4f} < {pl}. Skipping whole query.'.format(mn = prob_mn, pl = P_low))
            break

        qstart = alignment.qstart(qs, qe, qn, offset)     
        prob_mn *= qstart

        if alignment._prob * qstart * qend < P_low:
            if DEBUG:
                print ('|_P(q, t) * qend * qstart = {mn:.4f} < {pl}. Switching to next query fragment.'.format(mn = prob_mn, 
                    pl = P_low))
            continue

        if DEBUG:
            print ('|_Prob_mn = {mn:.4f}. Comparing to templates.'.format(q = query, 
                f = qfrag, qs = qstart, qe = qend, mn = prob_mn))
                
        tseen = list()

        for tfrag in sorted(alignment._pred[templ], key=lambda x: x[1]):

            # find the right boundaries 
            tbound = sorted([ p for p in alignment._pred[templ] if p[2] == tfrag[2] ])
        
            ts = tbound[0][0]
            te = tbound[-1][1]
            tn = tfrag[2]

            if tn in tseen:
                continue
            else:
                tseen.append(tn)
        
            tend = alignment.tend(ts, te, tn, offset)
            tmp = prob_mn * tend

            if tmp < P_low:
                if DEBUG:
                    print ('  |_P(q, t) * qend * qstart * tend = {mn:.4f} < {pl}. Switching to next query fragment.'.format(mn = tmp, 
                        pl = P_low))
                break

            tstart = alignment.tstart(ts, te, tn, offset)
            tmp *= tstart

            if tmp < P_low:
                if DEBUG:
                    print ('  |_P(q, t) * qend * qstart * tend * tstart = {mn:.4f} < {pl}. Switching to next template fragment.'.format(mn = tmp, 
                        pl = P_low))
                continue
            
            ali_sa = alignment.ali_starts_after(qs, ts, offset)
            ali_eb = alignment.ali_ends_before(qe, te, offset)
            
            tmp *= ali_sa
            tmp *= ali_eb

            if DEBUG:
                print ('  |_Prob_mn = {mn:.4f}.'.format(qs = tstart, qe=tend, mn=tmp))
            
            query_frag = query + "_" + str(qn)
            templ_frag = templ + "_" + str(tn)
            #if tmp > 1:
            #    import pdb; pdb.set_trace() 
            if tmp != 0:
                export.append((query_frag, templ_frag, tmp))

    if len(export) != 0:
        return export
    else:
        return None

def write_results(output, results, mode):
    
    with open(output, mode) as out_fh:
        for query in results:
            for entry in query:
                out_fh.write("{q}\t{t}\t{p:.4f}\n".format(q=entry[0], t=entry[1], p=entry[2]))

def read_graph(in_file):
    
    graph = defaultdict(dict)
    
    with open(in_file) as in_fh:
        for line in in_fh:
            strip_line = line.strip()
            initial_node = strip_line.split()[0]
            terminal_node = strip_line.split()[1]
            probability = float(strip_line.split()[2])
            
            if initial_node == terminal_node:
                print ('Removing seledge! {inn} -> {tn} with prob = {prob}.'.format(inn = initial_node, tn = terminal_node, prob = probability))
                continue
            
            if initial_node not in graph:
                graph[initial_node][terminal_node] = probability
                graph[terminal_node][initial_node] = probability
           
            if initial_node in graph:
                if terminal_node not in graph[initial_node]:
                    graph[initial_node][terminal_node] = probability
                    graph[terminal_node][initial_node] = probability
                if terminal_node in graph[initial_node]:
                    if probability > graph[initial_node][terminal_node]:
                        print ('Updating {inn} -> {tn} with prob = {prob}.'.format(inn = initial_node, tn = terminal_node, prob = probability))
                        graph[initial_node][terminal_node] = probability
                        graph[terminal_node][initial_node] = probability
                    else:
                        pass
    return graph

def arg():
    import argparse
    description = """Sim_graph creates a similarity graph from domain predictions."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('mat_data', help='Matrix ffdata input file') 
    parser.add_argument('mat_index', help='Matrix ffindex input file')
    parser.add_argument('predictions', help='File containing predictions') 
    parser.add_argument('output', help='Output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode')
    parser.add_argument('-p', '--p_low', help='P_low probability (default = 0.1)', type=float, default=0.1)
    parser.add_argument('-c', '--cov', help='Cov min (default = 0.8)', type=float, default=0.8)
    parser.add_argument('-m', '--cpu', help='Number of cores (default = 8).', type=int, default=8)
    parser.add_argument('-b', '--block', help='Blocksize (default = 32).', type=int, default=32)

    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])
    
    tmp_dir = tempfile.mkdtemp()
    tmp = os.path.join(tmp_dir, tmp_dir.split("/")[-1])
    tmp_graph = tmp + ".graph"

    if args.verbose:
        DEBUG = True

    entries = read_index(args.mat_index)
    data = read_data(args.mat_data)
    
    # load predictions 
    predictions = read_predictions(args.predictions)
    # process list (input for map)

    block_size = args.block
    current_block = 1
    pool = Pool(args.cpu)

    for block_start in range(0, len(entries), block_size):
        block_end = min(len(entries), block_start + block_size)
        block = entries[block_start:block_end]

        print ("Processing Block {cb} ({bs}-{be})".format(cb=current_block, bs=block_start, be=block_end))

        plist = list()
      
        for entry in block:
            query = entry.name.split(".")[0]

            if (entry.length == 1):
                print ("skip: "+ entry.name)
                continue

            entry_data = read_entry_data(entry, data)
            alis = read_alignment_matrices(entry.length, entry_data)

            for ali in alis.alignments:
                template = ali.template
                #print (query + '->' + template)
                
                if query == template:
                    print ("Skipping {0} (Self alignmnent).".format(template))
                    break
                
                prob = ali.alignment_probability / float(100)

                if prob < args.p_low:
                    if DEBUG:
                        print ("Not adding {q} - {t}. Prob(query, template) < P_low.".format(q=query, t=template))
                    break
                
                #if query == 'corrupt_mat' and template =='WEWKIWABA':
                #    import pdb; pdb.set_trace()
               

                pred = dict()
                pred[query] = predictions[query]
                pred[template] = predictions[template]
                
                # check whether all matrices contain values
                if (len(ali.alignment_start_matrix) == 0) or (len(ali.alignment_end_matrix) == 0) or (len(ali.alignment_posterior_matrix) == 0):
                    print ('! Warning detected empty matrices! Skipping '+ query + '->' + template)
                else:                
                    plist.append((AlignmentMatrix(query, template, prob, pred, ali.alignment_end_matrix, ali.alignment_start_matrix, ali.alignment_posterior_matrix), args.p_low, args.cov))
        
        # Uncomment to debug! 
        
        if args.verbose:
            processed_data = list()
            for process in plist:
                chunck = wrapper(process)
                processed_data.append(chunck)
        else:
            processed_data = pool.map(wrapper, plist)
        
        results = list(filter(lambda x: x != None, processed_data))
        
        mode = "a" if block_start > 0 else "w"
        print ("Writing Block {cb}.".format(cb=current_block))
        write_results(tmp_graph, results, mode)
        current_block += 1

    graph = read_graph(tmp_graph)

    with open(args.output, 'w') as out_fh:
        for initial_node in graph:
            for terminal_node in graph[initial_node]:
                out_fh.write('{inn}\t{on}\t{prob}\n'.format(inn = initial_node, on = terminal_node, prob = graph[initial_node][terminal_node]))
   

if __name__ == "__main__":
    main()
