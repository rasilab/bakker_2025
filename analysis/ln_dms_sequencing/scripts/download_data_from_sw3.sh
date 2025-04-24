FOLDER="../data"
mkdir -p $FOLDER
cd $FOLDER
aws s3 cp s3://fh-pi-subramaniam-a-eco/git_repos/iriboseq/analysis/rbakker/deepseq/2025-02-13_i84_ln_dms/ . --recursive

