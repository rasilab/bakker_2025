# useful libraries
import os
import pandas as pd
import re
import itertools as it


# configuration specific to this analysis
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment="#", dtype=object).set_index('sample_id')
print(sample_annotations)

# these rules are run locally
localrules: all

# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    joined_barcodes = expand('../data/joined_barcodes/{sample_id}.csv', 
      sample_id=sample_annotations.index),
    linked_barcode_recorder_file = expand('../data/linked_barcodes/{sample_id}.csv', 
      sample_id=sample_annotations.index),

def get_split_read_files_input(wildcards):
  """This function returns the names of R1,R2 files for combining them
  """
  filenames = [f'../data/fastq/{filename}' 
      for filename in filter(lambda x: re.search(f'{wildcards.sample_id}_', x) and re.search('_R[123]_', x) and x.endswith('.fastq'), os.listdir('../data/fastq'))]
  if len(filenames) == 0:
      raise FileNotFoundError(f"No input FASTQ files for {wildcards.sample_id}!") 
  if len(filenames) > 3:
      raise RuntimeError(f"Too many input FASTQ files for {wildcards.sample_id}: " + " | ".join(filenames)) 
  return list(sorted(filenames))

rule extract_and_tabulate_all_barcodes_recorders:
  """Extract barcodes and recorder regions and tabulate them for counting
  """
  input: 
    get_split_read_files_input
  output: '../data/joined_barcodes/{sample_id}.csv'
  params:
    barcode1_read = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'barcode_read'].tolist()[0],
    barcode1_start = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'barcode_start'].tolist()[0],
    barcode1_length = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'barcode_length'].tolist()[0],
    recorder_read = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'recorder_read'].tolist()[0],
    recorder_start = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'recorder_start'].tolist()[0],
    recorder_length = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'recorder_length'].tolist()[0]
  log: '../data/joined_barcodes/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell: 
    """
    set +o pipefail;
    # paste concatenates all input files line by line with tab separator
    paste {input} | awk -v OFS="," '
    BEGIN {{read=0;print "barcode1","recorder"}}
    NR % 4 == 2 \
    {{ 
      print substr(${params.barcode1_read}, {params.barcode1_start}, {params.barcode1_length}),\
            substr(${params.recorder_read}, {params.recorder_start}, {params.recorder_length});
      read++;
    }}' 1> {output}  2> {log}
    """

rule subset_to_linked_barcodes:
  """Subset barcode-recorder counts to only those barcodes identified in linkage sequencing
  """
  input:
    barcode_recorder_file = '../data/joined_barcodes/{sample_id}.csv',
  output:
    linked_barcode_recorder_file = '../data/linked_barcodes/{sample_id}.csv'
  params:
    linkage_ref_barcode_column = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'linkage_ref_barcode_column'].tolist()[0],
    barcode_linkage_file = lambda wildcards: sample_annotations.loc[sample_annotations.index == wildcards.sample_id, 'barcode_linkage_file'].tolist()[0],
  log: '../data/linked_barcode_counts/{sample_id}.log'
  container: 'docker://ghcr.io/rasilab/python:1.0.0'
  shell:
    """
    head -n 1 {input.barcode_recorder_file} > {output.linked_barcode_recorder_file}
    awk -F',' '
      BEGIN {{f={params.linkage_ref_barcode_column}}}
      NR == FNR {{e[$f] = 1; next}};
      e[$1]
      ' \
      {params.barcode_linkage_file} {input.barcode_recorder_file} 1>> {output.linked_barcode_recorder_file}
    """


