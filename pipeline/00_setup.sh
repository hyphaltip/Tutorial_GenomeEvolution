#!/usr/bin/bash

# setup some folders
mkdir -p analysis logs plots

pushd data
bash download.sh
popd
mkdir -p analysis/ortho_set1
pushd analysis/ortho_set1
ln -s ../../data/proteindb/*.fasta .
popd
