This ReadME file explains what the columns in annotations files mean.

### sample_annotations.csv

owner: initials of person who generated the libary
sample_name: descriptive name of library
sample_id: short reference name for library from lab notebooks
barcode1_read: illumina read where barcode is located
barcode1_start: which nucleotide in the read to start extracting barcode from
barcode1_length: how many nucleotides is the barcode
umi_read: illumina read where umi is located
umi_start: which nucleotide in read to start extracting umi from
umi_length: how many nucleotides is the umi
linkage_ref_barcode_column: which column in barcode_annotations.csv to match barcodes to

### barcodes_annotations.csv

barcode: sequence of barcode
reverse_complement: reverse complement of barcode sequence for matching reads
oligo_name: descriptive name of oligo
n_variants: how many variants in the oligo library with this barcode
variable_loc: which nucleotide in the read does the variable (NNN...) region start
variable_length: how many nucleotides is the variable region
target_loc: which nucleotide in the read does the a-rich target region start
target_length: how many nucleotides is the target region
target_pos_to_boxb: is the target region located 5' or 3' of the boxb hairpin
target_dist: how many nucleotides away from boxb hairpin is the a-rich target region
variable_type: what type of variable region is this- an a-rich variable region or boxb variable region
variable_subpos: for boxb variable libraries, what positions in the hairpin are varied
g_depleted: are the spacer regions the g_depleted sequence or not (there were concerns about secondary structures from g-quads)