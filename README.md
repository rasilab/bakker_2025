# Deaminase-based RNA recording enables high throughput mutational profiling of protein-RNA interactions

Rachael A. Bakker<sup>1,2</sup>, Oliver B. Nicholson<sup>1</sup>, Heungwon Park<sup>1,2</sup>, Yu-Lan Xiao<sup>3</sup>, Weixin Tang<sup>3</sup>, Arvind Rasi Subramaniam<sup>1,2,†</sup>, and Christopher P. Lapointe<sup>1,†</sup>

<sup>1</sup> Basic Sciences Division, Fred Hutchinson Cancer Center, Seattle, WA 98109, USA  
<sup>2</sup> Computational Biology Section of the Public Health Sciences Division,  Fred Hutchinson Cancer Center, Seattle, WA 98109, USA  
<sup>3</sup> Department of Chemistry and Institute for Biophysical Dynamics,  The University of Chicago, Chicago, IL 60637, USA  

<sup>†</sup> **Corresponding authors:** A.R.S: <rasi@fredhutch.org>, C.P.L: <cplapointe@fredhutch.org>



# Abstract

Protein-RNA interactions govern nearly every aspect of RNA metabolism and are frequently dysregulated in disease. 
While individual protein residues and RNA nucleotides critical for these interactions have been characterized, scalable methods that jointly map protein- and RNA-level determinants remain limited.
RNA deaminase fusions have emerged as a powerful strategy to identify transcriptome-wide targets of RNA-binding proteins by converting binding events into site-specific nucleotide edits. 
Here, we demonstrate that this 'RNA recording' approach enables high-throughput mutational scanning of protein-RNA interfaces.
Using the λN-boxB system as a model, we show that editing by a fused TadA adenosine deaminase directly correlates with binding affinity between protein and RNA variants *in vitro*.
Systematic variation of RNA sequence context reveals a strong bias for editing at UA dinucleotides by the engineered TadA8.20, mirroring wild-type TadA preferences.
We further demonstrate that stepwise recruitment of the deaminase using nanobody and protein A/G fusions maintains both sequence and binding specificity. 
Stable expression of the TadA fusion in human cells reproduces *in vitro* editing patterns across a library of RNA variants.
Finally, comprehensive single amino acid mutagenesis of λN in human cells reveals critical residues mediating RNA binding.
Together, our results establish RNA recording as a versatile and scalable tool for dissecting protein-RNA interactions at nucleotide and residue resolution, both *in vitro* and in cells.


## Instructions for running the code

- All software necessary for running the code are available as Docker images at https://github.com/orgs/rasilab/packages. These images can be downloaded to your local computer using [Docker](https://www.docker.com/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html). Most distributed clusters will already have Singularity available. You can also download and install [Singularity](https://anaconda.org/conda-forge/singularity) on your local computer using [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).
 
- We use [Snakemake](https://anaconda.org/bioconda/snakemake-minimal) for workflow management. This can be installed using Conda or might be already available in your distributed cluster.

- To reproduce the analysis on a cluster, load Singularity and Snakemake (or activate the Conda environment with these software). Ensure that all necessary folders are mounted using `--bind` in [analysis/submit_cluster.sh](./analysis/submit_cluster.sh) and [analysis/submit_local.sh](./analysis/submit_local.sh). These folder locations will be specific to your computing environment. If the correct location is not mounted, you will get `path not found` error in Snakemake workflows that use Singularity containers.

## Useful Docker containers for interactive analyses using Jupyter Notebooks

- [R and Python](https://github.com/rasilab/r_python/pkgs/container/r_python)
