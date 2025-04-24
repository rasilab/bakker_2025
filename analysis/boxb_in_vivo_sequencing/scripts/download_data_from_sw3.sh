FOLDER="../data"
mkdir -p $FOLDER
cd $FOLDER
aws s3 cp s3://fh-pi-subramaniam-a-eco/git_repos/iriboseq/analysis/rbakker/deepseq/2024-09-05_i82_in_vivo_boxb_oligo_pool/ . --recursive
