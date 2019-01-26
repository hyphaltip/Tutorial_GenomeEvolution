#!/usr/bin/bash
#SBATCH --ntasks 32 --mem 24G --nodes 1 --out logs/ortho_set1.log

module load OrthoFinder
module load ncbi-blast/2.8.1+
module load diamond

# process set1
orthofinder -M msa -A muscle -a 32 -t 32 -f analysis/ortho_set1

pushd analysis/ortho_set1
ln -s Results_*/Orthologs.csv .
