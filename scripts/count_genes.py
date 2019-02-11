#!/usr/bin/env python3

import csv, os, io, re, gzip

pepdir = 'data/proteindb'
outdir   = "reports"

if not os.path.exists(outdir):
    os.mkdir(outdir)

outcounts = os.path.join(outdir,"gene_counts.tab")
with open(outcounts,"w") as ofh:
    tabout = csv.writer(ofh,delimiter="\t",quoting=csv.QUOTE_MINIMAL)
    tabout.writerow(['SPECIES','GENE_COUNT'])
    for pepfile in os.listdir(pepdir):
        species = re.sub(r'\.aa\.fasta','',pepfile)
        with open(os.path.join(pepdir,pepfile)) as fa:
            unique_names = {}
            for line in fa:
                if line.startswith(">"):
                    pat = re.compile(r'^>(\S+).+')
                    m   = pat.search(line)
                    if( m ):
                        genename = m.group(1)
                        #re.sub(r'(\-t\d+_\d+\-p\d+|(\-T\d*)?\-p\d+)$',
                        #                  '',m.group(1))
                        unique_names[genename] = 1
            tabout.writerow([species,len(unique_names.keys())])
