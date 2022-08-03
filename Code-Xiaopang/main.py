#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import xiaopang as xp

filename='C:/Users/isynbio/Desktop/python_program/txt/E.txt'
with open(filename, 'r') as f1: data =''.join(f1.readlines())

[symbol_composition,symbol_frequency]=xp.p_generate(data) 
seat,codewords=xp.seat_and_codewords(symbol_frequency)
code_initial=xp.codegenerate(codewords,symbol_composition,data) 

huf_chart1=xp.enchart1_word(symbol_composition)  
huf_chart2=xp.enchart2_loc(seat)        
huf_chart01=xp.enchart0(symbol_composition)     
huf_chart02=xp.enchart0(seat)           
code_table_add=huf_chart01+huf_chart02+huf_chart1+huf_chart2+code_initial     

[code5,code6]=xp.c5c6()   
code_first_extension_mapping=xp.code_5to6(code_table_add,code5,code6)    
code_final=xp.encode_homopolymer(code_first_extension_mapping)   

dna_sequences=xp.f_index_and_num2base(code_final,994)
