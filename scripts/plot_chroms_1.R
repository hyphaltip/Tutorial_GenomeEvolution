library(ggplot2)
library(gridExtra)
library(tximport)
library(AnnotationDbi)
library(GenomicFeatures)
library(dplyr)
library(RColorBrewer)


source_version = "FungiDB-41"
species_list = read.csv("data/dna_species.dat",header=F,stringsAsFactors = FALSE)
colnames(species_list)=c("prefix","taxonomy")

drawGeneSummaryFunc <- function(pref) {
  taxonomy = sprintf("%s",species_list[species_list$prefix == pref,]$taxonomy)
  dbfile  =  sprintf("data/%s.db",pref)
  gff_file = sprintf("data/dnadb/%s.gff",pref)
  print(sprintf("prefix is %s, taxonomy is %s; gff is %s",pref,taxonomy,gff_file))
  
  if ( file.exists(dbfile) ) {
    txdb = loadDb(dbfile)
  } else {
    txdb <- makeTxDbFromGFF(gff_file,
                            dataSource= source_version,
                            organism = taxonomy)
    saveDb(txdb,dbfile)
  }
  ebg <- exonsBy(txdb, by="gene")
  gene_lengths = sum(width(reduce(ebg)))
  tbg<- transcriptsBy(txdb, by="gene")
  ibt <- intronsByTranscript(txdb)
  intronct <- elementNROWS(ibt)

  txlens <- transcriptLengths(txdb, with.cds_len=TRUE,
                                    with.utr5_len=TRUE,
                                    with.utr3_len=TRUE)

  ilens = width(ibt)
  ilens   = unlist(ilens[any(ilens>0)],use.name=FALSE)
  ict   = unlist(lapply(width(ebg),function(i) length(i)-1),use.name=FALSE)
  
  elens = unlist(width(ebg))
  ect   = unlist(lapply(width(ebg),function(i) length(i)),use.name=FALSE)
  summarystats = data.frame(
             species        = pref,
             intronlen_mean = mean(ilens),
             intronlen_sd  = sd(ilens),
             intronct_mean  = mean(ict),
             intronct_sd = sd(ict),
             exonlen_mean   = mean(elens),
             exonlen_sd    = sd(elens)
             )
  
  txdist <- ggplot(data=txlens,aes(txlens$tx_len)) + 
    geom_histogram(breaks=seq(0,10000,by = 100),fill="slategray") +
    labs(title=sprintf("%s Tx Length",pref),
         xlab="Length (bp)") + theme_minimal()
  txdist

  ChromTxs = transcripts(txdb)
  chrnames = gsub("_A_fumigatus\\S+","",seqnames(ChromTxs))
  chrnames = gsub("scf_0+","",chrnames)
  chrnames = gsub("Spom972h_","",chrnames)
  
  d = data.frame(start = start(ChromTxs), 
                 end    = end(ChromTxs),
                 chr    = chrnames,
                 strand = strand(ChromTxs),
                 txname =  ChromTxs$tx_name,
                 txid   =  ChromTxs$tx_id,
                 introncount = as.numeric(intronct),
                 length = as.numeric(width(ranges(ChromTxs)))
  )
  d <- d[order(d$chr, d$start), ]
  
  chrlens = do.call(rbind,sapply(sort(unique(d$chr)),function(i) data.frame(name=i,len=max(subset(d$end,d$chr == i))),
                                 simplify=FALSE, USE.NAMES=TRUE))
  topten = head(chrlens[order(-chrlens$len),],10)
  biggest_chroms = factor(topten$name)
  d = d[d$chr %in%  biggest_chroms, ]
  d = droplevels(d)
  
  intdist <- ggplot(data=d,aes(d$introncount)) + geom_histogram(binwidth =1, fill="maroon") +
    labs(title=sprintf("%s Intron Count",pref),xlab="Intron Count") + theme_minimal()
  print(intdist)
  
  
  d$index = rep.int(seq_along(unique(d$chr)), times = tapply(d$start,d$chr,length))

  d$pos=NA
  
  lastbase=0
  ticks = NULL
  minor = vector(,8)
  nchr = length(levels(d$chr))
  for (i in 1:nchr ) {
    if (i ==1) {
      d[d$index==i, ]$pos = d[d$index==i, ]$start
    } else {
      
      ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
      lastbase = lastbase + max(d[d$index==(i-1),"start"])
      minor[i] = lastbase
      d[d$index == i,"start"] =
        d[d$index == i,"start"]-min(d[d$index==i,"start"]) +1
      d[d$index == i,"end"] = lastbase
      d[d$index == i, "pos"] = d[d$index == i,"start"] + lastbase
    }
  }
  ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
  minorB <- tapply(d$end,d$index,max,probs=0.5)
  # reduce the number of chroms to print to the biggest 10

  p <- ggplot(d,aes(x=d$pos,y=d$length, color=chr))+geom_point(alpha=0.8,size=1,shape=16) + 
      scale_color_brewer(palette="Paired",type="seq") + scale_x_continuous(name="Chromosome", expand = c(0, 0),
                                                                         breaks = ticks,
                                                                         labels=(unique(d$chr))) +
      guides(fill = guide_legend(keywidth = 3, keyheight = 1))
    p + ggtitle(sprintf("%s Gene Length distribution",pref)) + ylab("Gene Length") + theme_minimal()
  print(p)

  p <- ggplot(d,aes(x=d$pos,y=d$introncount, color=chr))+geom_line(alpha=0.7,size=1) + 
    scale_color_brewer(palette="Paired",type="seq") + scale_x_continuous(name="Chromosome", expand = c(0, 0),
                                                                       breaks = ticks,
                                                                     labels=(unique(d$chr))) +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1)) +
    ggtitle(sprintf("%s Intron count",pref)) + ylab("Intron Count") + theme_minimal()
  print(p)
  
  return(summarystats)
}

speciesct = length(species_list)
print(species_list)
pdf(sprintf("%s/%s","plots","chrom_features.pdf"),onefile=TRUE,width=12)
sumstatslist = sapply(species_list$prefix,drawGeneSummaryFunc,simplify=FALSE, USE.NAMES=TRUE)
sumstats = do.call(rbind, sumstatslist)
pdf(sprintf("%s/%s","plots","summary_stats.pdf"),onefile=TRUE)

p <- ggplot(sumstats, aes(x=intronlen_mean,y=intronct_mean,label=species,color=species)) + geom_point() + 
  theme_minimal() + ggtitle("Intron size vs occurance in genes (mean)") + scale_color_brewer(palette="Set1")
print(p)
p <- ggplot(sumstats, aes(x=intronlen_mean,y=exonlen_mean,color=species)) + geom_point() + 
  theme_minimal() + ggtitle("Intron size vs exon size genes (mean)") + scale_color_brewer(palette="Set1")
print(p)


