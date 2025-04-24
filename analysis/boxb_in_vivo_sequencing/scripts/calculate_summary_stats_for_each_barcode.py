#!/usr/bin/env python
# coding: utf-8

# # Count number of each nucleotide within reads and the number of each nucleotide at each position for distinct inserts

# ## Import Libraries

# In[1]:


import pandas as pd
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools as it
import sys


# ## Define Paths to Files

# In[ ]:


# barcode = 'CCGCACCGCC'
# annotations_file = "../annotations/barcode_annotations_test.csv"
# fastq_file = f"../data/split_fastq/i79_p23/{barcode}_R1.fastq"
# fastq_umi = f"../data/split_fastq/i79_p23/{barcode}_R2.fastq"
# output_file = f"../data/summary_stats/testing/{barcode}.csv"

umi_start = 0
umi_length = 7
target_edit_loc = [1,2,4,7,9,12,13,15]
target_sequence = 'ATTATGGTGTAATTCTA'
target_5_var_loc = [8,10,11,14,16]
target_3_var_loc = [0,3,5,6,8]
primer_sequence= "AGCACAACA" #this is the beginning of the primer for amplifying, it is present in some truncated reads that need to be filtered out

annotations_file = sys.argv[1]
fastq_file=sys.argv[2]
fastq_umi=sys.argv[3]
output_file= sys.argv[4]
barcode = sys.argv[5]


# ## Load barcode annotations and subset to barcode that we want to process

# In[42]:


annotations = pd.read_csv(annotations_file).set_index('reverse_complement')

annotations_barcode = annotations.loc[barcode]


# ## Iterate through reads and parse the summary stats we want to collect

# In[ ]:


R1_file = FastqGeneralIterator(fastq_file)
R2_file = FastqGeneralIterator(fastq_umi)
count_table = dict()
umi_table = dict()
proper_upstream=annotations_barcode['upstream_seq']
proper_downstream=annotations_barcode['downstream_seq']


if pd.isna(proper_downstream): # ensures edge cases handle correctly
    proper_downstream = "" 

for read1, read2  in zip(R1_file, R2_file):
    umi = read2[umi_start:umi_length][1]
    target = read1[1][annotations_barcode['target_start']:(annotations_barcode['target_start'] + annotations_barcode['target_length'])]
    upstream_seq=read1[1][annotations_barcode['target_start']-10:annotations_barcode['target_start']]
    downstream_seq=read1[1][(annotations_barcode['target_start']+17):(annotations_barcode['target_start']+27)]
    


    if primer_sequence not in read1 and re.fullmatch(proper_upstream, upstream_seq) and (
    (downstream_seq == '' and proper_downstream == '') or 
    re.fullmatch(proper_downstream, downstream_seq)):

        
        insert = read1[1][annotations_barcode['variable_start']:(annotations_barcode['variable_start'] + annotations_barcode['variable_length'])]

        if insert not in count_table:
            count_table[insert] = {'umi_counts': 0}

        if insert not in umi_table:
            umi_table[insert] = dict()

        if umi not in umi_table[insert]:
            umi_table[insert][umi] = 1
            count_table[insert]['umi_counts'] += 1
        else:
            umi_table[insert][umi] += 1
            continue

        variable_target = ''.join(target[loc] for loc in target_edit_loc)

        for num,nt in it.product(range(len(variable_target) + 1), ['A','G','T','C','N']):
            key = f'num_{num}_{nt}'
            if  key not in count_table[insert]:
                count_table[insert][key] = 0
        
        for nt in ['A','G','T','C','N']:
            count = variable_target.count(nt)
            count_table[insert][f"num_{count}_{nt}"] += 1


        for num,nt in it.product(range(len(variable_target)), ['A','G','T','C','N']):
            key = f'pos_{num}_{nt}'
            if  key not in count_table[insert]:
                count_table[insert][key] = 0

        for pos,char in enumerate(variable_target):
            count_table[insert][f'pos_{pos}_{char}'] += 1
    else:
        pass


# ## Write summary stats to output file

# In[36]:


count_table_df = pd.DataFrame.from_dict(count_table, orient='index')
count_table_df.to_csv(output_file, index_label='insert')

