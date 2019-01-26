#!/usr/bin/bash

if [ ! -f FungiDB-41_AfumigatusAf293.gff ]; then
curl -O http://fungidb.org/common/downloads/Current_Release/AfumigatusAf293/gff/data/FungiDB-41_AfumigatusAf293.gff
fi
if [ ! -f FungiDB-41_AfumigatusA1163.gff ]; then
curl -O http://fungidb.org/common/downloads/Current_Release/AfumigatusA1163/gff/data/FungiDB-41_AfumigatusA1163.gff
fi
if [ ! -f FungiDB-41_AnovofumigatusIBT16806.gff ]; then
curl -O http://fungidb.org/common/downloads/Current_Release/AnovofumigatusIBT16806/gff/data/FungiDB-41_AnovofumigatusIBT16806.gff
fi

if [ ! -f FungiDB-41_AfumigatusAf293_AnnotatedProteins.fasta ]; then
	curl -O http://fungidb.org/common/downloads/Current_Release/AfumigatusAf293/fasta/data/FungiDB-41_AfumigatusAf293_AnnotatedProteins.fasta
fi
if [ ! -f FungiDB-41_AfumigatusA1163_AnnotatedProteins.fasta ]; then
curl -O http://fungidb.org/common/downloads/Current_Release/AfumigatusA1163/fasta/data/FungiDB-41_AfumigatusA1163_AnnotatedProteins.fasta
fi

if [ ! -f  FungiDB-41_AnovofumigatusIBT16806_AnnotatedProteins.fasta ]; then
	curl -O http://fungidb.org/common/downloads/Current_Release/AnovofumigatusIBT16806/fasta/data/FungiDB-41_AnovofumigatusIBT16806_AnnotatedProteins.fasta
fi

if [ ! -f FungiDB-41_AfumigatusA1163_Genome.fasta ]; then
curl -O http://fungidb.org/common/downloads/Current_Release/AfumigatusA1163/fasta/data/FungiDB-41_AfumigatusA1163_Genome.fasta
fi

if [ ! -f FungiDB-41_AfumigatusAf293_Genome.fasta ]; then
curl -O http://fungidb.org/common/downloads/Current_Release/AfumigatusAf293/fasta/data/FungiDB-41_AfumigatusAf293_Genome.fasta
fi

if [ ! -f FungiDB-41_AnovofumigatusIBT16806_Genome.fasta ]; then
curl -O http://fungidb.org/common/downloads/Current_Release/AnovofumigatusIBT16806/fasta/data/FungiDB-41_AnovofumigatusIBT16806_Genome.fasta
fi
