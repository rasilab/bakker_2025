#!/usr/bin/env python
"""
Calculate BoxB loop editing statistics from FASTQ sequencing data.

This script processes paired-end FASTQ files to analyze BoxB editing patterns
in specific loop regions, tracking UMI counts and nucleotide compositions.
"""

import pandas as pd
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict
import itertools as it
import sys


def initialize_count_dict(sequence_length):
    """Initialize count dictionary with all possible keys."""
    count_dict = {'umi_counts': 0}
    nucleotides = ['A', 'G', 'T', 'C', 'N']
    
    # Initialize nucleotide frequency counts
    for num, nt in it.product(range(sequence_length + 1), nucleotides):
        count_dict[f'num_{num}_{nt}'] = 0
    
    # Initialize position-specific counts
    for pos, nt in it.product(range(sequence_length), nucleotides):
        count_dict[f'pos_{pos}_{nt}'] = 0
    
    return count_dict


def main():
    # Parse command line arguments
    annotations_file = sys.argv[1]
    fastq_file = sys.argv[2]
    fastq_umi = sys.argv[3]
    output_file = sys.argv[4]
    barcode = sys.argv[5]
    
    # Define constants
    UMI_START = 0
    UMI_LENGTH = 7
    SPACER3_BOXB_START = 60
    SPACER5_BOXB_START = 12
    BOXB_LENGTH = 19
    BOXB_LOOP_POSITIONS = [7, 9, 10]
    NUCLEOTIDES = ['A', 'G', 'T', 'C', 'N']
    SEQUENCE_LENGTH = len(BOXB_LOOP_POSITIONS)
    
    # Load barcode annotations
    annotations = pd.read_csv(annotations_file).set_index('reverse_complement')
    annotations_barcode = annotations.loc[barcode]
    
    # Determine BoxB extraction parameters
    is_spacer3 = annotations_barcode['oligo_name'].startswith("spacer3")
    boxb_start = SPACER3_BOXB_START if is_spacer3 else SPACER5_BOXB_START
    boxb_end = boxb_start + BOXB_LENGTH
    
    # Initialize data structures
    count_table = defaultdict(lambda: initialize_count_dict(SEQUENCE_LENGTH))
    umi_sets = defaultdict(set)
    
    # Process paired FASTQ files
    with open(fastq_file, 'r') as f1, open(fastq_umi, 'r') as f2:
        R1_file = FastqGeneralIterator(f1)
        R2_file = FastqGeneralIterator(f2)
        
        for read1, read2 in zip(R1_file, R2_file):
            # Extract UMI and BoxB sequence
            umi = read2[1][UMI_START:UMI_START + UMI_LENGTH]
            boxb = read1[1][boxb_start:boxb_end]
            
            # Extract edited sequence from specific positions
            boxb_edited_sequence = ''.join(boxb[pos] for pos in BOXB_LOOP_POSITIONS)
            
            # Skip if UMI already seen for this sequence
            if umi in umi_sets[boxb_edited_sequence]:
                continue
            
            umi_sets[boxb_edited_sequence].add(umi)
            count_dict = count_table[boxb_edited_sequence]
            count_dict['umi_counts'] += 1
            
            # Count nucleotides by frequency and position in single pass
            for nt in NUCLEOTIDES:
                count = boxb_edited_sequence.count(nt)
                count_dict[f"num_{count}_{nt}"] += 1
            
            for pos, char in enumerate(boxb_edited_sequence):
                count_dict[f'pos_{pos}_{char}'] += 1
    
    # Convert to DataFrame and save
    count_table_df = pd.DataFrame.from_dict(dict(count_table), orient='index')
    count_table_df.to_csv(output_file, index_label='boxb_edited_sequence')


if __name__ == "__main__":
    main()