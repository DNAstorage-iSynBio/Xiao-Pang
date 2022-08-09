#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import xiaopang as xp

data_path = os.path.dirname(os.path.realpath(__file__))
filename = os.path.join(data_path, 'data','Natural selection and central dogma of molecular biology.txt')
with open(filename, 'r') as f1: data =''.join(f1.readlines())

# Quaternary Huffman direct encoding
quaternary_huffman=xp.Quaternary_huffman(data=data)
[code_initial,symbol_composition,seat]=quaternary_huffman.codegenerate()

# Code table addition 
code_table=xp.Code_table_addition(symbol_composition=symbol_composition,seat=seat)
code_table_add=code_table.f_huffman_table()+code_initial     

# Bases extension mapping
extension_mapping=xp.Extension_mapping(code_initial=code_table_add)
code_final=extension_mapping.f_extension_mapping_twice()

# Sequence division and index addition
segmentation_index=xp.Segmentation_index(sequences_initial=code_final,n_bases_payload=994)
dna_sequences=segmentation_index.f_index_and_num2base()
