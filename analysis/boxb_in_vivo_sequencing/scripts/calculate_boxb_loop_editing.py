#!/usr/bin/env python
# coding: utf-8

# ## Import Libraries

# In[ ]:


import re
import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools as it
import sys


# ## Define Paths to Files

# In[ ]:


# barcode = 'AAGACGCCCA'
# annotations_file = "../annotations/barcode_annotations.csv"
# fastq_file = f"../data/split_fastq/i79_p1/{barcode}_R1.fastq"
# fastq_umi = f"../data/split_fastq/i79_p1/{barcode}_R2.fastq"
# output_file = f"../data/summary_stats/i79_p1/{barcode}.csv"

umi_start = 0
umi_length = 7
target_edit_loc = [1,2,4,7,9,12,13,15]
target_sequence = 'ATTATGGTGTAATTCTA'
spacer3_boxb_start=60
spacer5_boxb_start=12
boxb_length= 19
boxb_loop_loc=[7,9,10]
primer_sequence= "AGCACAACA" #this is the beginning of the primer for amplifying, it is present in some truncated reads that need to be filtered out

annotations_file = sys.argv[1]
fastq_file=sys.argv[2]
fastq_umi=sys.argv[3]
output_file= sys.argv[4]
barcode = sys.argv[5]


# ## Load Barcode Annotations and Subset to Those We Want to Process

# In[19]:


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

for read1, read2  in zip(R1_file, R2_file):
    umi = read2[umi_start:umi_length][1]
    insert = read1[1][annotations_barcode['variable_start']:(annotations_barcode['variable_start'] + annotations_barcode['variable_length'])]

    upstream_seq=read1[1][annotations_barcode['target_start']-10:annotations_barcode['target_start']]
    downstream_seq=read1[1][(annotations_barcode['target_start']+17):(annotations_barcode['target_start']+27)]
    


    if primer_sequence not in read1 and re.fullmatch(proper_upstream, upstream_seq) and (
    (downstream_seq == '' and proper_downstream == '') or 
    re.fullmatch(proper_downstream, downstream_seq)):
        
        if annotations_barcode['oligo_name'].startswith("spacer3"):
            boxb = read1[1][spacer3_boxb_start:spacer3_boxb_start+boxb_length]
        elif annotations_barcode['oligo_name'].startswith("spacer5"):
            boxb = read1[1][spacer5_boxb_start:spacer5_boxb_start+boxb_length]

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

        boxb_loop = ''.join(boxb[loc] for loc in boxb_loop_loc)

        for num,nt in it.product(range(len(boxb_loop) + 1), ['A','G','T','C','N']):
            key = f'num_{num}_{nt}'
            if  key not in count_table[insert]:
                count_table[insert][key] = 0
        
        for nt in ['A','G','T','C','N']:
            count = boxb_loop.count(nt)
            count_table[insert][f"num_{count}_{nt}"] += 1


        for num,nt in it.product(range(len(boxb_loop)), ['A','G','T','C','N']):
            key = f'pos_{num}_{nt}'
            if  key not in count_table[insert]:
                count_table[insert][key] = 0

        for pos,char in enumerate(boxb_loop):
            count_table[insert][f'pos_{pos}_{char}'] += 1


# ## Write summary stats to output file

# In[22]:


count_table_df = pd.DataFrame.from_dict(count_table, orient='index')
count_table_df.to_csv(output_file, index_label='insert')

