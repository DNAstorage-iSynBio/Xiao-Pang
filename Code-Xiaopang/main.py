#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import xiaopang as xp

data_path = os.path.dirname(os.path.realpath(__file__))
filename = os.path.join(data_path, 'data','Natural selection and central dogma of molecular biology.txt')
with open(filename, 'r') as f1: data =''.join(f1.readlines())

# Encoding -- quaternary Huffman direct encoding
quaternary_huffman=xp.Quaternary_huffman(data=data)
[encode_initial,symbol_composition,seat]=quaternary_huffman.codegenerate()

# Encoding -- code table addition 
code_table=xp.Code_table_addition(symbol_composition=symbol_composition,seat=seat)
encode_table_add=code_table.f_huffman_table()+encode_initial     

# Encoding -- bases extension mapping
extension_mapping=xp.Extension_mapping(code_initial=encode_table_add)
encode_final=extension_mapping.f_extension_mapping_twice()

# Sequence division and index addition
segmentation_index=xp.Segmentation_index(sequences_initial=encode_final,n_bases_payload=994)
dna_sequences=segmentation_index.f_index_and_num2base()

# Decoding -- the inverse process of bases extension mapping
mapping_inverse=xp.Extension_mapping_inverse(sequences_silicon=encode_final)
decode_initial=mapping_inverse.mapping_inverse_twice()

# Decoding -- reconstruct the correspondence between codewords and characters
code_table_reading=xp.Code_table_reading(data=decode_initial)
[decode_payloads,decode_codewords,symbol_composition]=code_table_reading.f_code_table_reading()

# Decoding -- decode into original text or picture data
decode_result=xp.decode_txt_pic(decode_payloads,decode_codewords,symbol_composition)




