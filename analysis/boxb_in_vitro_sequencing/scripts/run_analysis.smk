# useful libraries
import os
import pandas as pd
import re
import itertools as it


# configuration specific to this analysis
sample_annotations = pd.read_table("../annotations/sample_annotations.csv", 
                                   sep=",", comment="#", dtype=object).set_index('sample_id')

barcodes = pd.read_table("../annotations/barcode_annotations.csv", 
                         sep=",", comment="#", dtype=object).set_index('reverse_complement')


# these rules are run locally
localrules: all


# Rules ----------------------------------------------------------------------

rule all:
  """List of all files we want at the end
  """
  input:
    split_fastq = expand("../data/split_fastq/{sample_id}/{barcode}_R{read}.fastq", 
                         sample_id=sample_annotations.index, 
                         barcode=barcodes.index, read=[1,2]),
    summary_stats = expand("../data/summary_stats/{sample_id}/{barcode}.csv", 
                         sample_id=sample_annotations.index, 
                         barcode=barcodes.index),
    summary_stats_combined = expand("../data/summary_stats_combined/{sample_id}.csv.gz", 
                         sample_id=sample_annotations.index),
    boxb_stats = expand("../data/boxb_editing_stats/{sample_id}/{barcode}.csv", 
                         sample_id=sample_annotations.index, 
                         barcode=barcodes.index),
    boxb_stats_combined = expand("../data/boxb_stats_combined/{sample_id}.csv.gz", 
                         sample_id=sample_annotations.index)
    # hairpins_free_energy = '../data/hairpins_free_energy.fold'


def get_split_read_files_input(wildcards):
  """This function returns the names of R1,R2 files for combining them
  """
  filenames = [f'../data/fastq/fastq/{filename}' 
      for filename in filter(lambda x: re.search(f'{wildcards.sample_id}_', x) and re.search('_R[123]_', x) and x.endswith('.fastq'), os.listdir('../data/fastq/fastq'))]
  if len(filenames) == 0:
      raise FileNotFoundError(f"No input FASTQ files for {wildcards.sample_id}!") 
  if len(filenames) > 3:
      raise RuntimeError(f"Too many input FASTQ files for {wildcards.sample_id}: " + " | ".join(filenames)) 
  return list(sorted(filenames))


rule split_fastq_by_barcode:
    input: get_split_read_files_input
    output: ['../data/split_fastq/{sample_id}/' + f'{barcode}_R1.fastq' for barcode in barcodes.index],
            ['../data/split_fastq/{sample_id}/' + f'{barcode}_R2.fastq' for barcode in barcodes.index]
    params:
        barcode_len = 10,
        barcode_file = "../annotations/barcode_annotations.csv",
        output_dir = "../data/split_fastq/{sample_id}"
    shell:
        """
        mkdir -p {params.output_dir}
        awk -v barcode_len={params.barcode_len} -v output_dir={params.output_dir} -v read2={input[1]} -f split_by_barcode.awk {params.barcode_file} {input[0]}
        """

rule calculate_summary_stats:
    input: 
        R1 = "../data/split_fastq/{sample_id}/{barcode}_R1.fastq",
        R2 = "../data/split_fastq/{sample_id}/{barcode}_R2.fastq",
        barcode_annotations = "../annotations/barcode_annotations.csv",
        notebook = "calculate_summary_stats_for_each_barcode.ipynb"
    output: "../data/summary_stats/{sample_id}/{barcode}.csv"
    container: "docker://ghcr.io/rasilab/python:1.0.0"
    shell:
        """
        jupyter nbconvert --to script {input.notebook}
        notebook={input.notebook}
        script="${{notebook/.ipynb/.py}}"
        python ${{script}} {input.barcode_annotations} {input.R1} {input.R2} {output} {wildcards.barcode}
        """


rule combine_summary_stats:
    input: ['../data/summary_stats/{sample_id}/' + f'{barcode}.csv' for barcode in barcodes.index],
        notebook = "combine_barcode_summary_stats.ipynb"
    output: "../data/summary_stats_combined/{sample_id}.csv.gz"
    container: "docker://ghcr.io/rasilab/r:1.0.0"
    params:
        input_folder = "../data/summary_stats/{sample_id}"
    shell:
        """
        jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
        notebook={input.notebook}
        script="${{notebook/.ipynb/.r}}"
        Rscript ${{script}} {params.input_folder} {output}
        """

rule calculate_boxb_loop_stats:
    input: 
        R1 = "../data/split_fastq/{sample_id}/{barcode}_R1.fastq",
        R2 = "../data/split_fastq/{sample_id}/{barcode}_R2.fastq",
        barcode_annotations = "../annotations/barcode_annotations.csv",
        notebook = "calculate_boxb_loop_editing.ipynb"
    output: "../data/boxb_editing_stats/{sample_id}/{barcode}.csv"
    container: "docker://ghcr.io/rasilab/python:1.0.0"
    shell:
        """
        jupyter nbconvert --to script {input.notebook}
        notebook={input.notebook}
        script="${{notebook/.ipynb/.py}}"
        python ${{script}} {input.barcode_annotations} {input.R1} {input.R2} {output} {wildcards.barcode}
        """

rule combine_boxb_stats:
    input: ['../data/boxb_editing_stats/{sample_id}/' + f'{barcode}.csv' for barcode in barcodes.index],
        notebook = "combine_barcode_summary_stats.ipynb"
    output: "../data/boxb_stats_combined/{sample_id}.csv.gz"
    container: "docker://ghcr.io/rasilab/r:1.0.0"
    params:
        input_folder = "../data/boxb_editing_stats/{sample_id}"
    shell:
        """
        jupyter nbconvert --to script --ExecutePreprocessor.kernel_name=ir {input.notebook}
        notebook={input.notebook}
        script="${{notebook/.ipynb/.r}}"
        Rscript ${{script}} {params.input_folder} {output}
        """

rule calculate_hairpin_energy:
  """Calculate delta-g for every possible hairpin in the library
  """
  input:
    hairpins="../tables/hairpin_list.fasta",
  output:
    hairpins_free_energy = '../data/hairpins_free_energy.fold'
  container: '/fh/scratch/delete90/subramaniam_a/user/rasi/singularity/viennarna_2.6.4.sif'
  shell:
    """
    RNAfold --noPS < "{input}" > "{output}"
    """