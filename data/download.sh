#!/usr/bin/bash

RELEASE="41"
URL="http://fungidb.org/common/downloads/release-${RELEASE}"
# get comparative proteins
PROTEINDIR=proteindb
mkdir -p $PROTEINDIR
while read SPECIES
do
	FILE="FungiDB-${RELEASE}_${SPECIES}_AnnotatedProteins.fasta"
	OUTFILE=$PROTEINDIR/${SPECIES}.aa.fasta
	if [ ! -f $OUTFILE ]; then
		curl -o $OUTFILE $URL/${SPECIES}/fasta/data/$FILE
	fi
done < protein_species.dat

DNADIR=dnadb
mkdir -p $DNADIR

IFS=,
while read SPECIES TAXONOMY
do
	DNAFILE="FungiDB-${RELEASE}_${SPECIES}_Genome.fasta"
	OUTDNA=$DNADIR/${SPECIES}.dna.fasta
	GFFFILE="FungiDB-${RELEASE}_${SPECIES}.gff"
	OUTGFF=$DNADIR/${SPECIES}.gff

	if [ ! -f $OUTDNA ]; then
		 curl -o $OUTDNA $URL/${SPECIES}/fasta/data/$DNAFILE
	fi
	if [ ! -f $OUTGFF ]; then
		curl -o $OUTGFF $URL/${SPECIES}/gff/data/$GFFFILE
	fi
done < dna_species.dat
