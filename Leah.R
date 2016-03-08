#Load required packages
require(ggplot2)
require(reshape2)
require(gridExtra)

#Import Eldholm MiSeq data
miseq <- read.table("/Users/Leah/Desktop/miseq_filtered.txt", header = F, sep = '\t')
colnames(miseq) <- c("pops", "fst", "fet", "afd", "pos", "gene.symbol", "gene", "description", "allele", "freqs", "cats")
head(miseq)

#Import Eldholm HiSeq data
hiseq <- read.table("/Users/Leah/Desktop/hiseq.txt", header=F, sep='\t')
colnames(hiseq) <- c("pops", "fst", "afd", "fet", "pos", "gene.symbol", "gene", "description", "allele", "freqs", "gene", "refAllele", "type", "strand", "nc", "aac", "cats")
head(hiseq)

#Import Eldholm Original Dat
Eld <- read.table("/Users/Leah/Desktop/Pepperell Lab/Eld.txt", header = T, sep = '\t', na.strings = "?")
eld <- Eld[,c(1:6,17:25)]

#Manipulate p-values from Fisher's Exact Test (multiple hypothesis correction)
miseq$p <- 10** - (miseq$fet)
miseq$q <- p.adjust(miseq$p, method = "fdr")
miseq$col <- with(miseq, factor(ifelse(q > 0.01, 0, 1)))

hiseq$p <- 10** - (hiseq$fet)
hiseq$q <- p.adjust(hiseq$p, method = "fdr")
hiseq$col <- with(hiseq, factor(ifelse(q > 0.01, 0, 1)))

#Combine data sets
AFD <- rbind(data.frame(samp="miseq", position = miseq$pos, afd = miseq$afd, colour = miseq$col),
             data.frame(samp="hiseq", position = hiseq$pos, afd = hiseq$afd, colour = hiseq$col))

#Manipulate allele frequencies into new datasets
miseq.freqs <- colsplit(miseq$freqs, ",", c(0,8,12,14,28,31,34,39,42))
miseq.freqs$pos <- miseq$pos
miseq.freqs$q <- miseq$q
mi.freqs <- subset(miseq.freqs, q <= 0.01)
mi.freqs.m = melt(mi.freqs, id = c("pos", "q"))
colnames(mi.freqs.m) <- c("pos", "q", "timePoint", "alleleFreq")

hiseq.freqs <- colsplit(hiseq$freqs, ",", c(0,8,12,14,28,31,34,39,42))
hiseq.freqs$pos <- hiseq$pos
hiseq.freqs$q <- hiseq$q
hi.freqs <- subset(hiseq.freqs, q <= 0.01)
hi.freqs.m = melt(hi.freqs, id = c("pos", "q"))
colnames(hi.freqs.m) <- c("pos", "q", "timePoint", "alleleFreq")

#########################################################################
#Identify if there are differences between the two confident datasets
setdiff(hi.freqs$pos, mi.freqs$pos)
#2155214 3850630 4046322

subset(hi.freqs, pos == 2155214)
subset(miseq.freqs, pos == 2155214)
#this snp is not found in the miseq data

subset(hi.freqs, pos == 3850630)
subset(miseq.freqs, pos == 3850630)
#this snp is found in the miseq data, at lower frequency and with less confidence

subset(hi.freqs, pos == 4046322)
subset(miseq.freqs, pos == 4046322)
#this snp is not found in the miseq data


setdiff(mi.freqs$pos, hi.freqs$pos)
# 888774 3843001   56321

subset(mi.freqs, pos == 888774)
subset(hiseq.freqs, pos == 888774)
#snp is found in both datasets, but not confidently changing allele freq in hiseq

subset(mi.freqs, pos == 3843001)
subset(hiseq.freqs, pos == 3843001)
#snp not found in hiseq data

subset(mi.freqs, pos == 56321)
subset(hiseq.freqs, pos == 56321)
#not identified in hiseq data

subset(mi.freqs, pos == 6579)
subset(hiseq.freqs, pos == 6579)

#Are these SNPs reported in Eldholm? Did they verify the SNPs?

bothSeq <- union(mi.freqs$pos, hi.freqs$pos)
allSeq <- intersect(bothSeq, Eld$Ref.Pos)

tab1 <- allSeq

posMifreq <- subset(miseq, pos %in% allSeq, c(freqs,pos))
posHifreq <- subset(hiseq, pos %in% allSeq, c(freqs,pos))
posEld <- subset(Eld, Ref.Pos %in% allSeq, c(SF1.1:SF9.1, Ref.Pos))

intersect <- intersect(hiseq$pos, miseq$pos)
union <- union(hiseq$pos, miseq$pos)

fetintersect <- subset(miseq, pos %in% union, fet)
mean(fetintersect$fet)

mean(miseq[miseq$pos %in% union,]$afd)
######################################################

#Combine allele trajectories into one dataframe
AT <- rbind(data.frame(samp="miseq", position = mi.freqs.m$pos, tp = mi.freqs.m$timePoint, af = mi.freqs.m$alleleFreq),
            data.frame(samp="hiseq", position = hi.freqs.m$pos, tp = hi.freqs.m$timePoint, af = hi.freqs.m$alleleFreq))


#AFD Plot
#tiff(filename = "Plot1.tiff", width = 13.18, height = 12, units = "cm", res = 100)

pAT <- ggplot(AT, aes(x=as.numeric(as.character(tp)),y=(af))) + 
  geom_line(aes(group=position)) +
  geom_point(aes(group=position), size=3, shape=20) +
  facet_wrap(~samp, nrow = 1) +
  ylab("Minor Allele Frequency") + 
  xlab("Time point (Months)") +
  scale_x_continuous(breaks=c(seq(0,42,by=6))) + 
  scale_y_continuous(limits=c(0,1)) +
  theme_set(theme_bw(base_size = 12)) +
  labs(title = "B.") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.0001, size = rel(2)),
        axis.title.y=element_text(vjust = 1.2, size = 16),
        axis.title.x=element_text(vjust = -0.2, size = 23))
pAT
#dev.off()

#tiff("Plot2.tiff", width = 13.18, height = 6.6, units = "cm", res = 50)
pAT <- ggplot(AT, aes(x=as.numeric(as.character(tp)),y=(af))) + 
  geom_line(aes(group=position)) +
  geom_point(aes(group=position), size=3, shape=20) +
  facet_wrap(~samp, nrow = 1) +
  ylab("Minor Allele Frequency") + 
  xlab("Time point (Months)") +
  scale_x_continuous(breaks=c(seq(0,42,by=6))) + 
  scale_y_continuous(limits=c(0,1)) +
  theme_set(theme_bw(base_size = 12)) +
  labs(title = "B.") +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.position="none",
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.0001, size = rel(2)),
        axis.title.y=element_text(vjust = 1.2, size = 16),
        axis.title.x=element_text(vjust = -0.2, size = 23))
pAT
#dev.off()

#Put the plots together
require(gridExtra)

p3 <- ggplot_gtable(ggplot_build(pAFD))
p4 <- ggplot_gtable(ggplot_build(pAT))


#And the final result is a stacked plot with perfectly aligned vertical axes. Success!
#tiff(filename = "Plot3.tiff", width = 13.18, height = 12, units = "cm", res = 100)
grid.arrange(p3, p4, ncol=1)
#dev.off()