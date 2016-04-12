#!/usr/bin/env python

import numpy as np
import sys
from ctypes import *
from collections import namedtuple
import ffindex

AllAlignments = namedtuple("AllAlignments", "query, query_length, alignments")
AlignmentData = namedtuple("AlignmentData", "index, template, template_length, alignment_probability, alignment_similarity, alignment_start_profile, alignment_end_profile, alignment_start_matrix, alignment_end_matrix, alignment_posterior_matrix")


def read_unsigned_short_int(data, index):
    num = data[index+1] | data[index] << 8
    return (num, index + 2)


def read_signed_short_int(data, index):
    num = data[index+1] | data[index] << 8
    num = c_short(num).value
    return (num, index + 2)


def bit_16_to_float(input):
    normalization_63 = 536870912

    exp_shift = 13
    mant_shift = 13

    m = input & 2046
    m = m << mant_shift

    e = input & 129024
    e = e << exp_shift
    e += normalization_63
    
    res = e | m
    res = cast(pointer(c_int32(res)), POINTER(c_float)).contents.value
    
    return res


def bit_8_to_float(input):
    normalization_111 = 939524096
    exp_shift = 19
    mant_shift = 19

    m = input & 15
    m = m << mant_shift

    e = input & 240
    e = e << exp_shift
    e += normalization_111

    res = e | m
    res = cast(pointer(c_int32(res)), POINTER(c_float)).contents.value
    
    return res


def read_name(data, index):
    start_index = index
    
    while(data[index] != 0):
        index += 1
    end_index = index
    
    name = data[start_index:end_index].decode("utf-8")
    
    return (name, end_index+1)


def process_profile(query_length, matrix):
    profile = np.zeros(query_length+1, dtype=np.float)
    
    for tuple in matrix:
        profile[tuple[0]] += tuple[2]

    profile /= sum(profile)
        
    return profile


def read_matrix(data, index):
    posteriors = []
    while True:
        (query_index, index) = read_unsigned_short_int(data, index)

        if(query_index == 0):
            break
        
        (template_index, index) = read_unsigned_short_int(data, index)

        while True:
            c = data[index]
            index += 1

            if(c == 0):
                break

            probability = bit_8_to_float(c)
            assert(probability >= 0.0 and probability <= 1.0)

            posteriors.append([query_index, template_index, probability])
            template_index += 1
        
    return (posteriors, index)
  

def read_alignment_matrices(length, entry_data):
    index = 0
    (query_name, index) = read_name(entry_data, index)
    # print ("Query:", query_name)

    if query_name.find("|") != -1:
        query_name = query_name.split("|")[1]

    
    (query_length, index) = read_unsigned_short_int(entry_data, index)
    
    actual_ali_index = -1
    
    alignments = []
    
    while(index < length - 2):
        actual_ali_index += 1

        try:
          (template_name, index) = read_name(entry_data, index)
          # print (template_name)
        
          (template_length, index) = read_unsigned_short_int(entry_data, index)
          
          alignment_probability = int(entry_data[index])
          index += 1
  
          (alignment_similarity_raw, index) = read_signed_short_int(entry_data, index)
          alignment_similarity = alignment_similarity_raw / 10.0
          
          (backward_matrix, index) = read_matrix(entry_data, index)
          (forward_matrix, index) = read_matrix(entry_data, index)
          (posteriors, index) = read_matrix(entry_data, index)
          
          if template_name.find("|") != -1:
            template_name = template_name.split("|")[1]
            
#           start_i = [x[0] for x in backward_matrix]
#           start_j = [x[1] for x in backward_matrix]
#           start_value = [x[2] for x in backward_matrix]
#           backwrd_matrix = (start_i, start_j, start_value)
          
          alignments.append(AlignmentData(actual_ali_index, template_name, template_length, alignment_probability, alignment_similarity, process_profile(query_length, backward_matrix).tolist(), process_profile(query_length, forward_matrix).tolist(), backward_matrix, forward_matrix, posteriors))
        except:
          break
        
    return AllAlignments(query_name, query_length, alignments)


def read_alignment_matrices_from_file(file_path):
  with open(file_path, "rb") as fh:
    data = fh.read()
    return read_alignment_matrices(len(data), data)


def main():
    data = ffindex.read_data(sys.argv[1])
    entries = ffindex.read_index(sys.argv[2])
    e = ffindex.get_entry_by_name('BAPNUNABA.a3m', entries)
    for entry in [e]:
        print(entry.name)
        if(entry.length == 1):
            print("skip: "+entry.name)
            continue
        entry_data = ffindex.read_entry_data(entry, data)
        alis = read_alignment_matrices(entry.length, entry_data)
        for ali in alis.alignments:
            print(ali.alignment_start_matrix)


if __name__ == '__main__':
    main()
