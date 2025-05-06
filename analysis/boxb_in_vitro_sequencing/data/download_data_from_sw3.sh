FOLDER="../data/fastq"
mkdir -p $FOLDER
cd $FOLDER
aws s3 cp s3://fh-pi-subramaniam-a-eco/git_repos/iriboseq/analysis/rbakker/deepseq/2024-09-05_i83_tada_lambdaN_pAG_oligo_pool/ . --recursive
