library(pheatmap)

domainset = 'reports/Pfam.genes_normalized.tab'
domainsetRaw = 'reports/Pfam.genes.tab'
pdf("plots/pheatmap_Pfam.pdf",height=15)
domainct = read.table(domainset,header=T,sep="\t",row.names=1)
domainctRaw = read.table(domainsetRaw,header=T,sep="\t",row.names=1)
colct = length(colnames(domainct))
summary(domainct)

colSums(domainct) # number of genes partipating in these families
domaintotals = rowSums(domainct)
domainct$totals = domaintotals
sorted_by_size = domainct[order(-domainct$totals),]


varMat = apply(domainct,1,var)
varMatS = varMat[order(-varMat)]

#top_100_var = head(varMatS,100)
domainct$variance = varMat
domain_order_by_variance = domainct[order(-domainct$variance),]
pheatmap(head(domain_order_by_variance[1:colct],75),main="Top 75 domains with highest Variance",
         clustering_distance_rows="correlation",scale="row")


domainctRaw$totals = rowSums(domainctRaw)
sortedraw_by_size = domainct[order(-domainctRaw$totals),]

pheatmap(head(sortedraw_by_size[1:colct],75),main="Top 75 domains by raw count",scale="row")

