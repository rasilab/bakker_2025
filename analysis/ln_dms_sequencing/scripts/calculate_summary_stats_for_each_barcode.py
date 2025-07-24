# Count number of each nucleotide within reads and the number of each nucleotide at each position for distinct inserts

# Import Libraries
import pandas as pd
import re
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import itertools as it
import sys
import os

# Define Paths to Files
annotations_file = "../annotations/sample_annotations.csv"
linked_barcode_recorder_file = f"../data/linked_barcodes/i84p1.csv"
output_file = f"../data/summary_stats/i84p1.csv"

target_edit_loc = [1,2,4,7,9,12,13,15]
target_sequence = 'ATTATGGTGTAATTCTA'
target_5_var_loc = [8,10,11,14,16]
target_3_var_loc = [0,3,5,6,8]

# Create output directory if it doesn't exist
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Read in Barcode Recorder File
edit_data = pd.read_csv(linked_barcode_recorder_file)

# Iterate through reads and parse the summary stats we want to collect
from collections import defaultdict
import numpy as np

count_table = defaultdict(lambda: defaultdict(int))

# Initialize all possible keys to avoid repeated checking
for barcode in edit_data['barcode1'].unique():
    count_table[barcode]['read_counts'] = 0
    # Pre-initialize all keys
    for num,nt in it.product(range(9), ['A','G','T','C','N']):  # 8 positions + 0
        count_table[barcode][f'num_{num}_{nt}'] = 0
    for num,nt in it.product(range(8), ['A','G','T','C','N']):  # 8 positions
        count_table[barcode][f'pos_{num}_{nt}'] = 0

print(f"Processing {len(edit_data)} rows...")

for i, (barcode, target) in enumerate(zip(edit_data['barcode1'], edit_data['recorder'])):
    if i % 100000 == 0:
        print(f"Processed {i} rows...")
    
    count_table[barcode]['read_counts'] += 1
    variable_target = ''.join(target[loc] for loc in target_edit_loc)

    # Count nucleotides more efficiently
    for nt in ['A','G','T','C','N']:
        count = variable_target.count(nt)
        count_table[barcode][f"num_{count}_{nt}"] += 1

    # Count by position more efficiently
    for pos, char in enumerate(variable_target):
        count_table[barcode][f'pos_{pos}_{char}'] += 1

# Write summary stats to output file
count_table_df = pd.DataFrame.from_dict(count_table, orient='index')
count_table_df.to_csv(output_file, index_label='barcode')

print(f"Summary statistics written to {output_file}")
print(f"Number of barcodes processed: {len(count_table)}")