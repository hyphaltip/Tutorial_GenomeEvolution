#!/usr/bin/bash

mkdir proteindb
RELEASE="41"
URL="http://fungidb.org/common/downloads/release-${RELEASE}"
# get comparative proteins
while read SPECIES
do
	FILE="FungiDB-${RELEASE}_${SPECIES}_AnnotatedProteins.fasta"
	OUTFILE=proteindb/${SPECIES}.aa.fasta
	if [ ! -f $OUTFILE ]; then
		curl -o $OUTFILE $URL/${SPECIES}/fasta/data/$FILE
	fi
done < protein_species.dat

mkdir dnadb
while read SPECIES
do
	DNAFILE="FungiDB-${RELEASE}_${SPECIES}_Genome.fasta"
	OUTDNA=dnadb/${SPECIES}.dna.fasta
	GFFFILE="FungiDB-${RELEASE}_${SPECIES}.gff"
	OUTGFF=dnadb/${SPECIES}.gff

	if [ ! -f $OUTDNA ]; then
		 curl -o $OUTDNA $URL/${SPECIES}/fasta/data/$DNAFILE
	fi
	if [ ! -f $OUTGFF ]; then
		curl -o $OUTGFF $URL/${SPECIES}/gff/data/$GFFFILE
	fi
done < dna_species.dat

popd
