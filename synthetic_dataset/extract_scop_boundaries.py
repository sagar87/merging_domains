#!/usr/bin/env python

from hh_reader import parse_result
from ffindex import read_data, read_index, read_lines
from re import match, search, findall
from collections import defaultdict

def extract_uniprot_id(template_info):

	old_pattern = findall(r'[OPQ][0-9][A-Z0-9]{3}[0-9]', str(template_info))
	new_pattern = findall(r'([A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})', str(template_info))

	uniprot_id = list()
	if len(new_pattern) > 0:
		for np in new_pattern:
			uniprot_id.append(np[0])

	if len(old_pattern) > 0:
		for op in old_pattern:
			uniprot_id.append(op)

	if uniprot_id:
		return ' '.join(uniprot_id)
	else:
		return ""


def extract_data(data, index, ev, cov_min, sim, max_len):

	domains = defaultdict(dict)

	for num, idx in enumerate(index):
		
		if num % 100 == 0:
			print ('Processing {num}/{total}'.format(num = num, total = len(index)))
		
		lines = read_lines(idx, data)
		hhr_data = parse_result(lines)

		for ali in hhr_data:

			cov = ali.aligned_cols / float(ali.query_length)

			# calculate whether the template is much larger than the query SCOP
			temp_len = ali.end[1] - ali.start[1] + 1 
			rel_size = (temp_len / float(ali.query_length)) -1
			
			if (cov > cov_min) and (ali.evalue < ev) and (ali.similarity > sim) and (rel_size < max_len):
				
				uniprot_id = extract_uniprot_id(ali.template_info)
				
				# print ('SCOP: {scop_domain} UNIPROT: {up} Qneff: {qneff} Tneff: {tneff} P: {prob} EV: {eval} Sc: {score}, ACols: {acol}, Id: {ident}, Sim: {sim}, SumPr: {sumProb}').format(
				# 	scop_domain = idx.name,
				# 	up = uniprot_id,
				# 	qneff = ali.query_neff,
				# 	tneff = ali.template_neff,
				# 	prob = ali.probability,
				# 	eval = ali.evalue,
				# 	score = ali.score,
				# 	acol = ali.aligned_cols,
				# 	ident = ali.identity,
				# 	sim = ali.similarity,
				# 	sumProb = ali.sum_probs)


				# print ('QL: {query_len} TL: {temp_len} Cov: {cov} Temp: {start}-{end} max_len: {max_len}'.format(
				# 	query_len = ali.query_length, 
				# 	temp_len = str(ali.end[1] - ali.start[1] + 1),
				# 	cov = cov,
				# 	start = ali.start[1],
				# 	end = ali.end[1],
				# 	max_len = rel_size ))
				
				# #print (idx.name, ali.aligned_cols, ali.similarity, ali.sum_probs)
				# #print (ali.template_info)
				
				scop = ali.query_id
				cluster = ali.template_id.split("|")[1]
				boundaries = (ali.start[1], ali.end[1])
				
				domains[scop][cluster] = boundaries

	return domains


def write_results(domains, out_file):

	# write mapping
	with open(out_file + '_map.dat', 'w') as out_fh:
		for scop in domains:
			for cluster in domains[scop]:
				out_fh.write('{scop}\t{cluster}\t{start}-{end}\n'.format(
					scop = scop,
					cluster = cluster,
					start = str(domains[scop][cluster][0]),
					end = str(domains[scop][cluster][1])))
	# write predictions, they need unique identifier
	seen_predictions = defaultdict(int)

	with open(out_file + '_pred.dat', 'w') as out_fh:
		for scop in domains:
			for cluster in domains[scop]:
				seen_predictions[cluster] += 1
				out_fh.write('{cluster}\t{cluster}_{idx}\t{start}-{end}\n'.format(
					cluster = cluster,
					idx = seen_predictions[cluster],
					start = str(domains[scop][cluster][0]),
					end = str(domains[scop][cluster][1])))



def arg():
    import argparse
    description = """Designed to extract SCOP domains from an HHblits search."""
    epilog= '"A nice quote" -'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('data', help = 'ffdata results') 
    parser.add_argument('index', help = 'ffindex results') #  nargs = '+' for more one argument
    parser.add_argument('out_file', help = 'results file')

    parser.add_argument('-ev', type = float, help = 'E value (default = 10**-3)', default = 10**-3)
    parser.add_argument('-cov', type = float, help = 'Cov (default = 0.9)', default = 0.9)
    parser.add_argument('-sim', type = float, help = 'Similarity (default = 0.6)', default = 0.6)
    parser.add_argument('-max_len', type = float, help = 'Maximum Overlength (default = 0.1)', default = 0.1)
    parser.add_argument('-pr', type = float, help = 'Threshold probability (default = 99.0)', default = 99)

    parser.add_argument('-v', '--verbose', action = 'store_true', help = 'verbose mode')

    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    data = read_data(args.data)
    index = read_index(args.index)

    dataset = extract_data(data, index, args.ev, args.cov, args.sim, args.max_len)
    write_results(dataset, args.out_file)
    

if __name__ == "__main__":
    main()