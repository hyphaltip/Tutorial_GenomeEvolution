#!/usr/bin/env python3

import csv, os, io, re, gzip

domaindir = 'domains'
outdir   = "reports"
if not os.path.exists(outdir):
    os.mkdir(outdir)

species_counts = "reports/gene_counts.tab"

species_ct = {}
with open(species_counts,'r') as sct:
    intsv = csv.reader(sct,delimiter="\t")
    header = next(intsv)
    for row in intsv:
        species_ct[row[0]] = int(row[1])
        
for domtype in os.listdir(domaindir):
    prefixes     = {}
    domaincounts = {}
    dirpath = os.path.join(domaindir,domtype)
    for dfile in os.listdir(dirpath):
        domfile = ""
        domstream = ''
        prefix = ""
        if dfile.endswith(".domtbl"):
            domfile = os.path.join(dirpath,dfile)
            domstream = open(domfile,'r')

        elif dfile.endswith(".domtbl.gz"):
            domfile = os.path.join(dirpath,dfile)
            domstream = io.TextIOWrapper(gzip.open(domfile,'r'),newline="")

        prefix = re.sub(r'\.domtbl(\.gz)?$',"",dfile)
        print(prefix)
        if prefix not in prefixes:
            prefixes[prefix] = 1
        unique_names = {}
        for line in domstream:
            if line.startswith("#"):
                continue
            row = line.split() # not quite right, need to specify number of cols
                               # as desc has spaces that should not be parsed 
            domain = row[0]
            domlen = row[2]
            pepname= row[3]

            #            genename = re.sub(r'(\-t\d+_\d+\-p\d+|(\-T\d*)?\-p\d+)$','',pepname)
            genename = pepname
            if genename not in unique_names:
                unique_names[genename] = 1
            else:
                # skip this entry because there are multiple isoforms
                # just take first isoform for the time being (longest might be best?)
                continue
#            print(genename, ' is genename')
            peplen = row[4]
            evalue = row[12]
#            print(domain,pepname,evalue)
            if domain not in domaincounts:                
                domaincounts[domain] = {prefix: {} }
            if prefix not in domaincounts[domain]:
                domaincounts[domain][prefix] = {}
                    
            if pepname not in domaincounts[domain][prefix]:                
                domaincounts[domain][prefix][pepname] = 1
            else:
                # add one to the counts for this domain in this gene
                domaincounts[domain][prefix][pepname] = 1 + domaincounts[domain][prefix][pepname]

    with open(os.path.join(outdir,"%s.total.tab"%(domtype)),"w") as domout:
        with open(os.path.join(outdir,"%s.genes.tab"%(domtype)),"w") as genesout:
            with open(os.path.join(outdir,"%s.genes_normalized.tab" %(domtype)),"w") as normgeneout:

                outtsv = csv.writer(domout,delimiter="\t",quoting=csv.QUOTE_MINIMAL)
                geneouttsv = csv.writer(genesout,delimiter="\t",quoting=csv.QUOTE_MINIMAL)
                normgenetsv = csv.writer(normgeneout,delimiter="\t",quoting=csv.QUOTE_MINIMAL)
                
                prefixorder = sorted( prefixes.keys())
                prefixorder.insert(0,"DOMAIN")
                #        prefixorder.append("DESCRIPTION")
                outtsv.writerow(prefixorder)
                geneouttsv.writerow(prefixorder)
                normgenetsv.writerow(prefixorder)
                
                for domain in sorted( domaincounts.keys()):
                    row = [domain]
                    generow = [domain]
                    normgenerow = [domain]
                    for p in prefixorder[1:]:
                        if p in domaincounts[domain]:
                            row.append(sum(domaincounts[domain][p].values()))
                            generow.append(len (domaincounts[domain][p].keys()))
                            normgenerow.append(len (domaincounts[domain][p].keys()) / species_ct[p])
                        else:
                            row.append(0)
                            generow.append(0)
                            normgenerow.append(0)
                        
                    outtsv.writerow(row)
                    geneouttsv.writerow(generow)
                    normgenetsv.writerow(normgenerow)

