library(ggplot2)
library(gridExtra)
library(tximport)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

library(GenomicFeatures)

source_version = "FungiDB-41"
species_list = read.csv("data/dna_species.dat",header=F)
colnames(species_list)=c("prefix","taxonomy")

drawSummaryFunc <- function(pref) {
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
  intronFrame = data.frame(
             prefix         = pref,
             intronlen_mean = mean(ilens),
             intronlen_var  = sd(ilens),
             intronct_mean  = mean(ict),
             intronct_var   = sd(ict))
  
  txdist <- ggplot(data=txlens,aes(txlens$tx_len)) + 
    geom_histogram(breaks=seq(0,10000,by = 100),fill="slategray") +
    labs(title=sprintf("%s Tx Length",pref),
         xlab="Length (bp)") + theme_minimal()
  txdist

  ChromTxs = transcripts(txdb)
  #gsub("_A_fumigatus\\S+","",seqnames(ChromTxs)),
  d = data.frame(start = start(ChromTxs), 
                       end    = end(ChromTxs),
                       chr    = seqnames(ChromTxs),
                       strand = strand(ChromTxs),
                       txname =  ChromTxs$tx_name,
                       txid   =  ChromTxs$tx_id,
                       introncount = intronct,
                      length = width(ranges(ChromTxs)))
  d <- d[order(d$chr, d$start), ]


  intdist <- ggplot(data=d,aes(d$introncount)) + geom_histogram(binwidth=1,fill="maroon") +
    labs(title=sprintf("%s Intron Count",pref),xlab="Intron Count") + theme_minimal()
  intdist

  d$index = rep.int(seq_along(unique(d$chr)), times = tapply(d$start,d$chr,length))

  d$pos=NA
  nchr = length(levels(d$chr))
  lastbase=0
  ticks = NULL
  minor = vector(,8)
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
  p <- ggplot(d,aes(x=d$pos,y=d$length, color=chr))+geom_point(alpha=0.8,size=1,shape=16) + 
      scale_color_brewer(palette="Set2",type="seq") + scale_x_continuous(name="Chromosome", expand = c(0, 0),
                                                                         breaks = ticks,
                                                                         labels=(unique(d$chr))) +
      guides(fill = guide_legend(keywidth = 3, keyheight = 1))
    p + labs(title=sprintf("%s Gene Length distribution",pref)) + ylab("Gene Length")
  p

  p <- ggplot(d,aes(x=d$pos,y=d$introncount, color=chr))+geom_line(alpha=0.7,size=1) + 
    scale_color_brewer(palette="Set2",type="seq") + scale_x_continuous(name="Chromosome", expand = c(0, 0),
                                                                       breaks = ticks,
                                                                     labels=(unique(d$chr))) +
    guides(fill = guide_legend(keywidth = 3, keyheight = 1)) +
    ggtitle(sprintf("%s Intron count",pref)) + ylab("Intron Count")
  p
  
  return(intronFrame)
}

pdf(sprintf("%s/%s","plots","chrom_distributions.pdf"),onefile=TRUE,width=12)
speciesct = length(species_list)
print(species_list)
intronparams <- sapply(species_list$prefix,drawSummaryFunc)
dev.off()



