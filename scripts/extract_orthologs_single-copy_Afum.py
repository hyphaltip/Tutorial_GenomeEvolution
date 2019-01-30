#!/usr/bin/env python3

import csv, sys, os

infolder = "analysis/ortho_set1/Results"
orthofile = "Orthogroups.csv"
gene_count = "Orthogroups.GeneCount.csv"

ingroup_species = [ "AfumigatusA1163.aa",
                    "AfumigatusAf293.aa" ]

if len(sys.argv) > 1:
    # analysis folder is not the default
    infolder = sys.argv[1]

orthogroups_singlecopy = {}

with open(os.path.join(infolder,gene_count)) as genecount:
    group_ct = csv.reader(genecount,delimiter="\t",quoting=csv.QUOTE_MINIMAL)
    header = next(group_ct)
    query_columns = []

    i = 0
    for h in header:
        if h in ingroup_species:
            query_columns.append(i)
        i = i + 1
    # sanity check that header has our test species
    if len(query_columns) != len(ingroup_species):        
        print("Found number of header columbs (%d) != ingroup lits (%d)"
              %(len(query_columns), len(ingroup_species)))        
        exit()

    
    # iterate through the orthologs.GeneCount.csv file
    for orthogroup in group_ct:
        # check see if all our ingroups == 1 gene in this family
        last_sp = 0
        single_copy_this_orthogroup = 0
        for col in query_columns:
            if int(orthogroup[col]) == 1:
                # keep the streak alive - if 
                single_copy_this_orthogroup = 1
            else:
                # so value is not 1 so let's bail
                single_copy_this_orthogroup = 0
                break
        if single_copy_this_orthogroup:
            orthogroups_singlecopy[orthogroup[0]] = []


with open(os.path.join(infolder,orthofile)) as genenames:
    orthodat = csv.reader(genenames,delimiter="\t",quoting=csv.QUOTE_MINIMAL)
    header = next(orthodat)
    query_columns = []
    # same code as before as I don't want to assume columns are in same order
    # but of course they probably are
    i = 0
    for h in header:
        if h in ingroup_species:
            query_columns.append(i)
        i = i + 1
    # sanity check that header has our test species
    if len(query_columns) != len(ingroup_species):        
        print("Found number of header columbs (%d) != ingroup lits (%d)"
              %(len(query_columns), len(ingroup_species)))        
        exit()

    # each row is an orthogroup with the names of genes in it
    for orthorow in orthodat:
        # this will restrict analysis to only the rows which have
        # the single copy genes
        if orthorow[0] in orthogroups_singlecopy:
            gene_row = [orthorow[0]]
            for col in query_columns:
                # now by definition this should only have a single value
                # we should perhaps convert this to a string rather than
                # array
                gene_row.append(orthorow[col].split(", ")[0])
            print(",".join(gene_row))
