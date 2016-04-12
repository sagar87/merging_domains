#!/usr/bin/env python

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
    tmp_sg = tmp + ".sg"

    print ("Mapping identifier to integers...")

    mapping = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "cc/mapping.py"), in_file, tmp_map]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = mapping.communicate()

    print ("Identifying total number of nodes...")
    # Get the number of total nodes
    nodes = subprocess.Popen(" ".join(["wc", "-l", tmp_map]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = nodes.communicate()

    num_nodes = int(output.split()[0]) - 1

    print ("Total number of nodes {0}. Reformatting graph...".format(num_nodes))

    reformat = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "cc/reformat.py"), tmp_map, in_file, tmp_grp]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = reformat.communicate()

    print ("Identifying connected components.")

    cc = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "cc/quick_find.py"), tmp_grp, tmp_cc, str(num_nodes)]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = cc.communicate()    

    print ("Writing connected components...")

    wc = subprocess.Popen(" ".join([os.path.join(os.getcwd(), "cc/write_graph.py"), tmp_map, tmp_cc, tmp_grp, tmp_sg]), 
        cwd=os.getcwd(), stderr=subprocess.STDOUT, stdout=subprocess.PIPE, shell=True)
    
    output, error = wc.communicate() 

    print ("Sorting edges by component...")
    
    s = subprocess.Popen(" ".join(["sort", "-k", "1", "-n", tmp_sg, ">", out_folder]),     
        cwd=os.getcwd(), shell=True)
    
    s.communicate()
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
    main()
