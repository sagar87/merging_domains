#!/usr/bin/env python
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
    main()