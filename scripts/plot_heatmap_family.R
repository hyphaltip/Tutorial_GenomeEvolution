library(pheatmap)

dataset = 'analysis/ortho_set1/Results/Orthogroups.GeneCount.csv'
pdf("plots/pheatmap_example.pdf")
ct = read.table(dataset,header=T,sep="\t",row.names=1)

summary(ct)
length(ct) # how many columns
length(ct$Moryzae70.15.aa) # how many rows looking at one column
length( subset(ct$Moryzae70.15.aa, ct$Moryzae70.15.aa > 0 )) # how many non-zero entries for this taxa?

colSums(ct) # number of genes partipating in these families

# get the first 40 rows
truncate_40 = head(ct,40)

# restruct to some logical 
bigger_50 = subset(ct,ct$Total > 50)

# tells you how many rows were returned
length(bigger_50$AfumigatusA1163.aa)

# now we don't print the Total column as it is a summary column, not actual data 
# so just take the first 9 columns from this dataset
bigger_50 = bigger_50[1:9]
pheatmap(bigger_50)

# top 100 
top_100 = head(ct,100)
top_100 = top_100[1:9]
pheatmap(top_100)


# the Foxysporum is blowing out the scale so we can't really see any subtly elsewhere
max(top_100$Foxysporum)
top_100 = subset(top_100,top_100$Foxysporum < 150)
pheatmap(top_100)

# lets look at top 1000
top_1000 = head(subset(ct[1:9],ct$Foxysporum < 150),1000)
pheatmap(top_1000)


pheatmap(top_1000,cluster_rows = FALSE,show_rownames=FALSE)
pheatmap(top_1000,cluster_cols = FALSE,show_rownames=FALSE)

#only Afum
Afum = top_1000[,c("AfumigatusA1163.aa","AfumigatusAf293.aa")]
pheatmap(Afum,cluster_cols = FALSE,show_rownames=FALSE)

ThreeFungi = ct[,c("AfumigatusA1163.aa","AfumigatusAf293.aa","NcrassaOR74A.aa")]
ThreeFungiBig = subset(ThreeFungi,ThreeFungi$NcrassaOR74A.aa > 5)
pheatmap(ThreeFungiBig,cluster_cols = FALSE,show_rownames=TRUE)
# only keep rows where the family is missing in Ncrassa
ThreeFungiMissingNcrassa = subset(ThreeFungi,ThreeFungi$NcrassaOR74A.aa == 0 & ThreeFungi$AfumigatusAf293.aa > 1)
pheatmap(ThreeFungiMissingNcrassa,cluster_cols = FALSE,show_rownames=FALSE)
pheatmap(ThreeFungiMissingNcrassa[1:2],cluster_cols = FALSE,show_rownames=FALSE)
