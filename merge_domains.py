#!/usr/bin/env python

from copy import deepcopy
from collections import defaultdict
from itertools import groupby

class NotStringError(ValueError):
	pass

class NotIntegerError(ValueError):
	pass

class WrongBoundaryError(ValueError):
	pass

## The cool thing about this script

class Domain(object):
    """A doubly linked list implementation of protein domain."""

    def __init__(self, protein, domain, start, end):
        """ Creates a domain object. """

        if type(protein) != str:
            raise NotStringError("A protein is a string!")
        if type(domain) != int:
            raise NotIntegerError("A domain is an integer!")
        if type(start) != int or type(end) != int:
            raise NotIntegerError("Start and end boundaries have to be integers!")
        if start > end:   
        	raise WrongBoundaryError("Domain boundaries for %s seem to be wrong: Start (%s) < End (%s)." % (protein, start, end))

        self._protein = protein
        self._domain = domain
        self._start = start
        self._end = end
        
        self._upstream = None
        self._downstream = None

    def __str__(self):
        """ Returns a string representation of domain and all its downstream domains. """
        ans = self._protein + "\t" + str(self._start) + "-" + str(self._end) + ","
        
        if self._downstream == None:
            return ans.strip(",")
        else:
            temp = self._downstream
            while temp != None:
                ans += str(temp._start) + "-" + str(temp._end) + ","
                temp = temp._downstream

        return ans.strip(",")

    def __iter__(self):

    	here = self
    	while here:
    		yield here
    		here = here._downstream

    def _verify_protein(self, other):
    	""" Simply checks if proteins are the same. """
    	return self._protein == other._protein

    def _verify_domain(self, other):
    	""" A invariant that verfies that other domain does not interfere with domain boundaries. """
    	
    	if other._start >= self._start and other._start <= self._end: # Start of other domain lies between the boundary of this domain
    		return False
    	elif other._end >= self._start and other._end <= self._end: # End position of other domain lies between the boundary of this domain
    		return False
    	else: # no interference
    		return True 

    def _is_upstream(self, other):
        """ Returns True iff this domain is directly upstream from the other domain. """
        # Before we do anything check invariants
        assert self._verify_protein(other), "Proteins have to the same."
        assert self._verify_domain(other), "Other domain conflicts domain boundaries."

        if self._end < other._start:
        	return True
        else:
        	return False

    def _is_downstream(self, other):
        """ Returns True iff this domain is directly downstream as compared to the other domain. """ 
        assert self._verify_protein(other), "Proteins have to the same."
        assert self._verify_domain(other), "Other domain conflicts domain boundaries."

        if other._end < self._start:
            return True
        else:
            return False
    
    def protein(self):
    	""" Returns the name of the protein. """
    	return self._protein

    def reference(self):
    	""" This information is used in the reference to retrieve the positional information of the 
    	given protein. """
    	return (self._domain, self._start, self._end)

    def first(self):
        """ Moves upwards to the most upstream domain present in the protein. """
        # Base case
        if self._upstream == None:
            return self
        else:
            # traverse until we find the first domain
            temp = self._upstream
            while temp._upstream != None:
                temp = temp._upstream

        return temp

    def last(self):
        """ Moves downwards to the most downstream domain present in the protein. """
        # Base case
        if self._downstream == None:
            return self
        else:
            # traverse until the very end
            temp = self._downstream
            while temp._downstream != None:
                temp = temp._downstream

        return temp

    def find_neighbors(self, other):
    	""" Returns the neighboring domains (left and right) of other domain. """
        # Before we do anything check invariants
        assert self._verify_protein(other), "Proteins have to the same."
        assert self._verify_domain(other), "Other domain conflicts domain boundaries."


        other_first = other.first()
        other_last = other.last()

        # get the neighbors
    	prev = self._upstream
    	next = self._downstream

    	# Other domain is inbetween the right domains
    	if other_first._is_downstream(self) and (next == None or other_last._is_upstream(next)):
    		return self, next

    	# Other domain is downstream, however, the next domain is not upstream from it 
    	elif other_first._is_downstream(self) and (next != None or not other_last._is_upstream(next)):
    		return next.find_neighbors(other)
    	
    	# Other domain is between previous domain and this domain
    	elif other_last._is_upstream(self) and (prev == None or other_first._is_downstream(prev)):
    		return prev, self
    	
    	# Other domains is upstream from the previous domain, hence we traverse backwards
    	elif other_last._is_upstream(self) and (prev != None or not other_first._is_downstream(next)):
    		return prev.find_neighbors(other)

    	# This is the case when the other domain encloses the this domain
    	elif not other_first._is_downstream(self) and not other_last._is_upstream(self):
    		return None, None

    def merge_domains(self, other):
        """ Merges two domains iff other domain is either directly downstream or upstream of the 
        protein. """
        left, right = self.find_neighbors(other)

        if left != None and right != None:
        	# link other to the upstream domain
        	left._downstream = other
        	other.first()._upstream = left
        	# link other to the downstream domain
        	other.last()._downstream = right
        	right._upstream = other

        # when other is the new first domain      
        elif left == None and right != None: 
        	other.last()._downstream = right
        	right._upstream = other

        # when other is the new last domain
        elif right == None and left != None:  
            left._downstream = other
            other.first()._upstream = left 

    def copy(self):
    	""" Returns a copy of the domain. """
    	here = self.first()
    	return deepcopy(here)

# The "Core" functions

def create_reference(reference):
	"""
	Yields a dictionary which maps for each proteins and all of its domains the position to the
	position in the protein.

	Example:
	a	a_1	1-20,200-240
	a	a_2	60-80
	a	a_3	100-150

	The structure of this protein "a" is schematically: 
	|1 domain 1 20|-|60 domain 2 80|-|100 domain 3 150|-|200 domain 1 240|

	Thus create_reference creates a the following dictionary:

	{'a': {(1, 1, 20): 0, (2, 60, 80): 1, (3, 100, 150): 2, (1, 200, 240): 3}}
	"""

	reference_dic = defaultdict(dict)

	for protein in reference:

		idx = 1 # is the integer that specifies the position within the protein

		for domain in reference[protein]:
			# does the magic
			reference_dic[domain._protein][(domain._domain, domain._start, domain._end)] = idx
			idx += 1

	return reference_dic

def is_included(left, domain, right, reference):
	""" Returns True iff domain is between left and right. """
	if left == None or right == None:
		return False

	pos_left = reference[left.protein()][left.reference()]
	pos_domain_left = reference[domain.first().protein()][domain.first().reference()]
	
	pos_domain_right = reference[domain.last().protein()][domain.last().reference()]
	pos_right = reference[right.protein()][right.reference()]

	return pos_left + 1 == pos_domain_left and pos_domain_right == pos_right - 1

def is_downstream(left, domain, right, reference):
	""" Returns True iff domain is downstream of left. """
	if left == None or right != None:
		return False

	pos_left = reference[left.protein()][left.reference()]
	pos_domain = reference[domain.first().protein()][domain.first().reference()]

	return pos_left + 1 == pos_domain 

def is_upstream(left, domain, right, reference):
	""" Returns True iff domain is upstream of right. """
	if right == None or left != None:
		return False

	pos_right = reference[right.protein()][right.reference()]
	pos_domain = reference[domain.last().protein()][domain.last().reference()]

	return pos_right - 1 == pos_domain

def merge_cluster(cluster_a, cluster_b, reference, relation):
	""" Takes two clusters a and b, a reference and a realation (0-3) and returns the upstream 
	probability between of the two clusters. """
	temp = list()
	merge = list()
	
	fcn = { 0:is_upstream, 1:is_downstream, 2:is_included, 3:is_included }

	# identify the domains which are directly upstream of domain and copy the to temp
	for domain_a in cluster_a:
		for domain_b in cluster_b:
			if domain_a.protein() == domain_b.protein():

				# treat the enclosed case specially
				if relation == 3:
					left, right = domain_a.find_neighbors(domain_b)

					if fcn[relation](left, domain_b, right, reference):
						temp.append(domain_a)
						temp.append(domain_b)

						copy_a = domain_a.copy()
						copy_b = domain_b.copy()
						# merge the domains
						copy_a.merge_domains(copy_b)
						# append the marged domains
						merge.append(copy_a.first())

				# if domain is directly upstream, downstream or included we'll proceed in the following way
				else:
					left, right = domain_b.find_neighbors(domain_a)
					
					if fcn[relation](left, domain_a, right, reference):
						temp.append(domain_a)
						temp.append(domain_b)

						copy_a = domain_a.copy()
						copy_b = domain_b.copy()
						# merge the domains
						copy_b.merge_domains(copy_a)
						# append the marged domains
						merge.append(copy_b.first())

	# create two new clusters
	new_a = list()
	new_b = list()

	for domain in cluster_a:   
		if domain not in temp:
			new_a.append(domain)

	for domain in cluster_b:
		if domain not in temp:
			new_b.append(domain)

	return new_a, new_b, merge

def probability(cluster_a, cluster_b, reference, verbose = False):
	""" Takes two clusters a and b, and a reference and returns the upstream probability between
	of the two clusters. """
	# upstream, downstream, included, enclosed
	candidates = [0, 0, 0, 0]
	
	for domain_a in cluster_a:

		for domain_b in cluster_b:
			# if we have the same protein in the cluster
			if domain_a.protein() == domain_b.protein():

				left, right = domain_b.find_neighbors(domain_a)
				
				# Check the four cases !
				if left == None and right == None:
					# try to find out whether the domain b encloses domain a
					left, right = domain_a.find_neighbors(domain_b)

					if is_included(left, domain_b, right, reference):
						if verbose: print "Domain %s encloses domain %s." % (domain_a, domain_b)
						candidates[3] += 1
					else:
						if verbose: print "Domain %s does NOT encloses domain %s." % (domain_a, domain_b)

				elif left == None and right != None:
					if is_upstream(left, domain_a, right, reference):
						if verbose: print "Domain %s is directly upstream of %s." % (domain_a, domain_b)
						candidates[0] += 1
					else:
						if verbose: print "Domain %s is NOT upstream of %s." % (domain_a, domain_b)

				elif left != None and right == None:
					if is_downstream(left, domain_a, right, reference):
						if verbose: print "Domain %s is directly downstream of %s" % (domain_a, domain_b)
						candidates[1] += 1
					else:
						if verbose: print "Domains %s is NOT directly downstream of %s." % (domain_a, domain_b)

				elif left != None and right != None:
					if is_included(left, domain_a, right, reference):
						if verbose: print  "Domain %s is inbetween of %s and %s" % (domain_a, left, right)
						candidates[2] += 1
					else:
						if verbose: print  "Domain %s is NOT inbetween of %s and %s" % (domain_a, left, right)

	probabilities = list()
	
	for candidate in candidates:
		if candidate != 0:
			probabilities.append(candidate/float(len(cluster_b)))
		else:
			probabilities.append(0)

	return probabilities

def cluster_intersect(cluster_a, cluster_b):
	""" Returns True iff the clusters have at least one protein in common, otherwise False. """

	proteins_a = set([domain.protein() for domain in cluster_a])
	proteins_b = set([domain.protein() for domain in cluster_b])

	intersection = proteins_a.intersection(proteins_b)
	
	if len(intersection) > 0:
		return True
	else:
		return False

def merge_iteration(clusters, reference, cut_off, verbose = False):
	""" Iterates one time over all clusters. """

	merged_domains = dict() # store new clusters here
	merge_idx = 0


	for cl_a in sorted(clusters):
		for cl_b in sorted(clusters):
			if cl_a != cl_b and cluster_intersect(clusters[cl_a], clusters[cl_b]):
				
				if verbose:
					print "#### Comparing Clusters %s and %s ####" % (cl_a, cl_b)
					print "Old A:", print_cluster(clusters[cl_a])
					print "Old B:", print_cluster(clusters[cl_b])
					print
					print "#### Calculating P(%s|%s) ####" % (cl_a, cl_b)
				
				# Calculate probabilities
				prob = probability(clusters[cl_a], clusters[cl_b], reference, verbose) 

				if verbose:
					print "\n#### Calculated probabilities ####"
					print "P_up = N(%s directly upstream in %s)/N(%s) = %s" % (cl_a, cl_b, cl_b, prob[0])
					print "P_down = N(%s directly downstream in %s)/N(%s) = %s" % (cl_a, cl_b, cl_b, prob[1])
					print "P_inc = N(%s included in %s)/N(%s) = %s" % (cl_a, cl_b, cl_b, prob[2])
					print "P_enc = N(%s encloses %s)/N(%s) = %s\n" % (cl_a, cl_b, cl_b, prob[3])

				max_prob = max(prob)
				relation = prob.index(max(prob))

				if max_prob > cut_off:
					new_a, new_b, merge = merge_cluster(clusters[cl_a], clusters[cl_b], reference, relation)				

					if verbose:
						print "#### Newly formed cluster ####"
						print "New A:", print_cluster(new_a)
						print "New B:", print_cluster(new_b)
						print "Merge:", print_cluster(merge)
						print

					clusters[cl_a] = new_a
					clusters[cl_b] = new_b

					# If we merged some domains add a new cluster which is stored in another dicitonary
					if len(merge) != 0:
						merged_domains[merge_idx] = merge
						merge_idx += 1

	# remove empty clusters
	empty_cluster = [ key for key in clusters if len(clusters[key]) == 0 ]

	for key in empty_cluster:
		del clusters[key]

	# merge old clusters with merged cluster
	clusters.update(merged_domains)
	
	return merge_idx, clusters

def wrapper(clusters, reference, cut_off, verbose = False):
	""" Calls merge_iteration until no new merges occur. """

	merges = 1
	iterations = 0


	while merges != 0:
		print "Working on iteration %s" % iterations
		merges, clusters = merge_iteration(clusters, reference, cut_off, verbose)
		iterations += 1

	return iterations, clusters 

## Parser functions

def merge_fragments(fragments):
	""" A helper function that merges a single domain which consists of several sequence fragments. 
	Only required to correctly parse data. """
	
	if len(fragments) == 1:
		return fragments.pop()

	domain = fragments.pop(0)

	while len(fragments) != 0:
		domain.merge_domains(fragments.pop(0))

	return domain.first()

def parse_domains(in_file):
	""" Parses predictions such that for every domain prediction a domain object is created which 
	can be accesssed in a dictionary of the form { <protein>:{<idx>:<domain>} }. """
	domains = defaultdict(dict)

	with open(in_file) as in_fh:
		
		for line in in_fh:
			strip_line = line.strip()
			protein, idx = strip_line.split()[1].split("_")
			fragments = strip_line.split()[2].split(",")

			temp_list = list()
			
			for fragment in fragments: 
				start, end = fragment.split("-")
				temp_list.append(Domain(protein, int(idx), int(start), int(end))) 

			domains[protein][int(idx)] = merge_fragments(temp_list)

	return domains

def parse_cluster(in_file, domains):

	clusters = defaultdict(list)
	with open(in_file) as in_fh:

		for line in in_fh:
			striped_line = line.strip()
			cluster = striped_line.split()[0]
			protein, idx = striped_line.split()[1].split("_")
			domain = domains[protein][int(idx)]
			clusters[cluster].append(domain)

	return clusters


def parse_predictions(in_file):

	proteins = defaultdict(list)

	with open(in_file) as in_fh:
		for cluster, members in groupby(in_fh, lambda line: line.split()[0]):
			
			for member in members:
				split_line = member.split()
				domain = split_line[1]
				name = domain.split("_")[0]
				idx = int(domain.split("_")[1])
				boundaries = split_line[2].split(",")

				fragments = list()

				for fragment in boundaries:
					start = int(fragment.split("-")[0])
					end = int(fragment.split("-")[1])
					fragments.append(Domain(name, idx, start, end))
					proteins[name].append(Domain(name, idx, start, end))

	predictions = dict()
	
	for protein in proteins:
		predictions[protein] = merge_fragments(proteins[protein])

	return predictions

def print_cluster(cluster):
	ans = ""
	for domain in cluster:
		ans += str(domain.first()).replace("\t", " ") + "\n"
	return ans.strip(", ")

def write_merged_cluster(out_file, merged_cluster):

	with open(out_file, 'w+') as out_fh:

		for cluster in sorted(merged_cluster):
			result = str(cluster) + "\t" + print_cluster(merged_cluster[cluster])
			# import pdb; pdb.set_trace()
			out_fh.write(result)

def sort_proteins(merged_cluster):
	"""
	Takes the merged clusters as input and sorts domains.
	"""
	proteins = defaultdict(dict)

	cluster_num = 1
	for cluster in merged_cluster:
		for domain in merged_cluster[cluster]:
			if domain._protein in proteins.keys():
				domain_number = max(proteins[domain._protein].keys()) + 1
				proteins[domain._protein][domain_number] = domain
			else:
				proteins[domain._protein][1] = domain

	return proteins

def write_results(out_file, sorted_proteins):

	with open(out_file, 'w') as out_fh:
		for protein in sorted_proteins:
			for domain in sorted_proteins[protein]:
				identifier = protein + '_' + str(domain) 
				boundaries = str(sorted_proteins[protein][domain].first()).split("\t")[1]
				print boundaries
				out_fh.write('{protein}\t{identifier}\t{boundaries}\n'.format(protein = protein,
					identifier = identifier, boundaries = boundaries))



def arg():
    import argparse
    description = """Merge domains uses the domain references and the clusters and merges 
    domain fragments until no more merges can occur."""
    epilog= '"These pointers in C were like fucking pain in the ass." - Donghan Lee'
    # Initiate a ArgumentParser Class
    parser = argparse.ArgumentParser(description = description, epilog = epilog)
    
    # Call add_options to the parser
    parser.add_argument('cluster', help='AP clusters.') 
    parser.add_argument('predictions', help='Predictions.')
    parser.add_argument('out_file', help='Output file.')

    parser.add_argument('-c', help='Cut off probability for cluster merge (default = 0.98).', type=float, default=0.98)
    parser.add_argument('-v', '--verbose', action='store_true', help='verbose mode', default=False)

    return parser


def main():
    import sys
    parser = arg()
    args = parser.parse_args(sys.argv[1:])

    reference = create_reference(parse_predictions(args.predictions))
    print len(reference)
    domains = parse_domains(args.predictions)
    clusters  = parse_cluster(args.cluster, domains)
    
    # print "#### Clusters BEFORE Merging domains ####"    
    # for cl in clusters:
    # 	print cl + " -> " +print_cluster(clusters[cl])
    # print

    merges, new_cluster = wrapper(clusters, reference, args.c, args.verbose)
    sorted_domain = sort_proteins(new_cluster)
    write_results(args.out_file, sorted_domain)

    write_merged_cluster(args.out_file + 'old_output.dat', new_cluster)

    # print "#### Clusters AFTER Merging domains (Iterations: %s) ####" % merges
    # for cl in new_cluster:
    # 	print str(cl) + " -> " + print_cluster(new_cluster[cl])

if __name__ == "__main__":
    main()

