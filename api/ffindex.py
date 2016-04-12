#!/usr/bin/env python

'''
Created on Apr 30, 2014

@author: meiermark
'''


import sys
import mmap
import os
from collections import namedtuple

FFindexEntry = namedtuple("FFindexEntry", "name, offset, length")


def read_index(ffindex_filename):
    entries = []
    
    fh = open(ffindex_filename)
    for line in fh:
        tokens = line.split("\t")
        entries.append(FFindexEntry(tokens[0], int(tokens[1]), int(tokens[2])))
    fh.close()
    
    return entries


def read_data(ffdata_filename):
    fh = open(ffdata_filename, "r")
    data = mmap.mmap(fh.fileno(), 0, prot=mmap.PROT_READ)
    fh.close()
    return data


def get_entry_by_name(name, index):
    #TODO: bsearch
    for entry in index:
        if(name == entry.name):
            return entry
    return None


def read_lines(entry, data):

    try:
        lines = data[entry.offset:entry.offset + entry.length - 1].decode("utf-8").split("\n")
    except UnicodeDecodeError: 
        # Added 6.1.2016: To prevent that the script crashes !
        return None

    return lines


def read_entry_data(entry, data):
    return data[entry.offset:entry.offset + entry.length - 1]


def write_entry(entries, data_fh, entry_name, offset, data):
    data_fh.write(data[:-1])
    data_fh.write(bytearray(1))

    entry = FFindexEntry(entry_name, offset, len(data))
    entries.append(entry)

    return offset + len(data)


def write_entry_with_file(entries, data_fh, entry_name, offset, file_name):
    with open(file_name, "rb") as fh:
        data = bytearray(fh.read())
        return write_entry(entries, data_fh, entry_name, offset, data)


def finish_db(entries, ffindex_filename, data_fh):
    data_fh.close()
    write_entries_to_db(entries, ffindex_filename)


def write_entries_to_db(entries, ffindex_filename):
    entries = sorted(entries, key=lambda x: x.name)
    index_fh = open(ffindex_filename, "w")

    for entry in entries:
      index_fh.write("{name:.64}\t{offset}\t{length}\n".format(name=entry.name, offset=entry.offset, length=entry.length))

    index_fh.close()


def write_entry_to_file(entry, data, file):
    lines = read_lines(entry, data)
    
    fh = open(file, "w")
    for line in lines:
        fh.write(line+"\n")
    fh.close()


def opt():
    import optparse
    usage = "usage: .py -arg1 "
    description = """A description"""
    # Initiate a OptionParser Class
    parser = optparse.OptionParser(usage=usage, description=description)
    
    # Call add_options to the parser
    parser.add_option('-f', help='ffdata file.', dest='data')
    parser.add_option('-i', help='ffindex file.', dest='index')


    return parser


def main():
    parser = opt()
    (options, argv) = parser.parse_args()

    data = read_data(options.data)
    index = read_index(options.index)
    test = read_lines(index[1], data)
    print(len(index))
    print(index[1])
    for l in test:
        print(l)

if __name__ == "__main__":
    main()