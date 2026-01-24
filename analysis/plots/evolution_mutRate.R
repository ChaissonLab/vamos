options(scipen = 999)
library(ggplot2)
library(ggExtra)
library(scales)
library(ggrepel)
library(reshape2)


######################################## process data ########################################
######################################################################
########## Tree                                                      #
######################################################################
##### all trees after filter hotspot #####
trees = read.csv('all.tree.filterHotspot.tsv', sep='\t', header=F)
colnames(trees) = c('chr','start','end','nSNP','treeLen','aveEd','aveDist')
trees$nSNP = trees$nSNP + 1
chrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
         'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
         'chr21','chr22')
trees$chr = factor(trees$chr, levels=chrs)
labels = data.frame(chr=chrs, totalSNP=0, totalTree=0, ave=0)
for (chr in chrs) {
  labels$totalSNP[labels$chr==chr] = sum(trees$nSNP[trees$chr==chr])
  labels$totalTree[labels$chr==chr] = nrow(trees[trees$chr==chr,])
  labels$ave[labels$chr==chr] = round( labels$totalSNP[labels$chr==chr] / labels$totalTree[labels$chr==chr], 0 )
}


######################################################################
########## Mutation rate                                             #
########## the input does not contain sex chromosomes                #
######################################################################
##### genomic loci #####
rate = read.csv('all.mutRate.filterHotspot.tagSegdup.tsv', sep='\t', header=T)
rate$key = paste0(rate$chr,'_',rate$start,'_',rate$end)
# create the categorical motif period variable
rate$period6 = rate$period; rate$period6[rate$period6 > 6] = '>6'
rate$period6 = factor(rate$period6, levels=c('1','2','3','4','5','6','>6'))
# filter not picked trees
rate = rate[rate$spanned != 'not-unpick',]
# change spanned tags to "full", "rescued", and "noTree"
rate$spanned[rate$TR2Tree > 10000] = 'noTree'
rate$spanned[rate$spanned == 'not-pick'] = 'rescued'

# merge catalog info
#catalog = read.csv('../vamosExpanded_v3.0_effMotifs-0.1_entropy_purity.tsv', sep='\t', header=F)
catalog = read.csv('vamosExpanded_v3.0_effMotifs-0.1_entropy_purity.tsv', sep='\t', header=F)
colnames(catalog) = c('chr','start','end','motifs','source','type','period','consensus','consGC','allGC',
                      'op1','op2','segDup','coding','exon','transcribed','HAR','entropy','purityOri','purityEff')
catalog$key = paste0(catalog$chr,'_',catalog$start,'_',catalog$end)

# merge TR catalog features into mutation rate data and rework some variables
rate = merge(rate, catalog[,c('source','consensus','consGC','allGC','op1','op2','segDup','coding','exon',
                              'transcribed','HAR','entropy','purityOri','purityEff','key')], by='key')
rate$coding[rate$coding=='partial'] = 'cds'
rate$segDup[rate$segDup=='segDup'] = 'in-segDup'
rate$segDup[rate$segDup=='partial'] = 'part-in-segDup'
rate$segDup[rate$segDup=='.'] = 'not-in-segDup'
rate$size = rate$end - rate$start
rate$snpDensity = rate$treeLen / rate$nSNP

#rate$entropy[rate$allGC < 0.5] = -rate$entropy[rate$allGC < 0.5] # change AT rich loci entropy to negative

##### bio loci #####
rateBio = read.csv('all.mutRateBio.filterHotspot.tagSegdup.tsv', sep='\t', header=T)
rateBio$key = paste0(rateBio$chr,'_',rateBio$start,'_',rateBio$end)
# create the categorical motif period variable
rateBio$period6 = rateBio$period; rateBio$period6[rateBio$period6 > 6] = '>6'
rateBio$period6 = factor(rateBio$period6, levels=c('1','2','3','4','5','6','>6'))
# filter not picked trees
rateBio = rateBio[rateBio$spanned != 'not-unpick',]
# change spanned tags to "full", "rescued", and "noTree"
rateBio$spanned[rateBio$TR2Tree > 10000] = 'noTree'
rateBio$spanned[rateBio$spanned == 'not-pick'] = 'rescued'

catalogBio = read.csv('../vamosBioTR_v3.0_effMotifs-0.1_entropy_purity.tsv', sep='\t', header=F)
colnames(catalogBio) = c('chr','start','end','motifs','source','type','period','consensus','consGC','allGC',
                         'op1','op2','segDup','coding','exon','transcribed','entropy','purityOri','purityEff')
catalogBio$key = paste0(catalogBio$chr,'_',catalogBio$start,'_',catalogBio$end)

rateBio = merge(rateBio, catalogBio[,c('source','consensus','consGC','allGC','op1','op2','segDup','coding','exon',
                                       'transcribed','entropy','purityOri','purityEff','key')], by='key')
rateBio$coding[rateBio$coding=='partial'] = 'cds'
rateBio$segDup[rateBio$segDup=='segDup'] = 'in-segDup'
rateBio$segDup[rateBio$segDup=='partial'] = 'part-in-segDup'
rateBio$segDup[rateBio$segDup=='.'] = 'not-in-segDup'
rateBio$size = rateBio$end - rateBio$start
rateBio$snpDensity = rateBio$treeLen / rateBio$nSNP

#rateBio$entropy[rateBio$allGC < 0.5] = -rateBio$entropy[rateBio$allGC < 0.5] # change AT rich loci entropy to negative

######################################################################
# Variability                                                        #
# config variability data and merge it to mutation rate              #
# note that there are locus with variability measure, but no         #
# mutation measure                                                   #
# also all samples are used for variability, but some may be removed #
# for mutation rate due to poor SNP coverage in the window           #
######################################################################
##### genomic loci #####
alleles = read.csv('../2025-03-02_assembly_stats/alleles/samples_rmDup_rmKids_GRCh38_q-0.1.annoStr.tsv', sep='\t', header=T)
alleles$key = paste0(alleles$chr,'_',alleles$start,'_',alleles$end)
rate = merge(rate, alleles[,c(6:14,60)], by='key')

#####  bio loci #####
allelesBio = read.csv('../2025-03-02_assembly_stats/alleles/samples_rmDup_rmKids_GRCh38_q-0.1.annoStr.special.tsv', sep='\t', header=T)
allelesBio$key = paste0(allelesBio$chr,'_',allelesBio$start,'_',allelesBio$end)
rateBio = merge(rateBio, allelesBio[,c(6:14,60)], by='key')

##### summary statistics #####
sum(alleles$lenAll == 1); sum(alleles$comAll == 1)

######################################################################
# config block mutation data and merge it to mutation rate           #
# only the genomic loci have block mutation data                     #
######################################################################
##### genomic loci #####
block = read.csv('all.mutRate.filterHotspot.tagBlock.closestRecomb.tsv', sep='\t', header=T)
colnames(block)[9] = 'dist2Recom'
block$key = paste0(block$chr,'_',block$start,'_',block$end)
block$aveBlockCat = 'block-no'
block$aveBlockCat[block$aveBlock > -1 & block$aveBlock < 4] = 'block-low'
block$aveBlockCat[block$aveBlock >= 4] = 'block-high'

# merge block data to rate data
rate = merge(rate, block[,4:12], by='key')
# merge catalog data to block data
block = merge(block, catalog[,6:ncol(catalog)], by='key')


######################################################################
# config parsimony data and merge everything to it                   #
# only the genomic loci have parsimony mutation data                 #
######################################################################
parsimony = read.csv('all.parsimonyBlock.merged.tsv', sep='\t', header=T)
parsimony$key = paste0(parsimony$chr,'_',parsimony$start,'_',parsimony$end)

parsimony = merge(parsimony, rate[,c(1,5:ncol(rate))], by='key')


######################################################################
# config disease data and merge it to the rateBio                    #
# only the bio loci have disease data                                #
######################################################################
#####
disease = read.csv('../disease.tsv', sep='\t', header=T)
disease$key = paste0(disease$chr,'_',disease$refined_start,'_',disease$refined_end)

rateBio = merge(rateBio, disease[,c(1,4,5,10)])

######################################################################
# summary statistics of mutation rate data                           #
######################################################################
##### basic stats #####
nrow(catalog[catalog$chr != 'chrX' & catalog$chr != 'chrY',]) # no. of autosomal loci (surveyed loci)
print(56349 + 224) # number of tree-filtered loci (large tree or low SNP density)
4050875 - nrow(rate) # number loci due to hotspot filtering
nrow(rate[rate$spanned != 'noTree',]) # number of successful loci with computed mutation rate
nrow(rate[rate$spanned == 'full',]) # number of fully spanned loci
nrow(rate[rate$spanned != 'noTree',]) / nrow(catalog[catalog$chr != 'chrX' & catalog$chr != 'chrY',])
nrow(rate[rate$spanned == 'full',]) / nrow(catalog[catalog$chr != 'chrX' & catalog$chr != 'chrY',])
nrow(catalog[catalog$chr != 'chrX' & catalog$chr != 'chrY',]) - nrow(rate[rate$spanned != 'noTree',]) # unresolved loci
# ratio of unresolved loci to the number of surveyed loci
(nrow(catalog[catalog$chr != 'chrX' & catalog$chr != 'chrY',]) - nrow(rate[rate$spanned != 'noTree',])) /
  nrow(catalog[catalog$chr != 'chrX' & catalog$chr != 'chrY',])
# unresolved loci due to poor coverage
nrow(catalog[catalog$chr != 'chrX' & catalog$chr != 'chrY',]) - nrow(rate[rate$spanned != 'noTree',]) - (56349 + 224)


##### STR vs. VNTR #####
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp1 = temp$rateNormAfter[temp$period > 6]
temp2 = temp$rateNormAfter[temp$period <= 6]
t.test(temp1, temp2); mean(temp1) / mean(temp2) # t test, how much more mutable is VNTR

##### check for stable loci with high mutation rates due to small tree branches
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp$aveDist = temp$totalDist / (temp$annotated*(temp$annotated-1) / 2)
temp$rateCat = 'rateHigh'; temp$rateCat[log10(temp$rateNormAfter) < -4] = 'rateLow'
temp$branchCat = 'average branch length \u2264 10000'; temp$branchCat[temp$aveDist > 10000] = 'average branch length > 10000'
temp$branchCat = factor(temp$branchCat, levels=c('average branch length \u2264 10000','average branch length > 10000'))
sum(temp$branchCat == 'average branch length \u2264 10000') # no. of loci with average branch length <= 10000
sum(temp$branchCat == 'average branch length > 10000') # no. of loci with average branch length > 10000
sum(temp$branchCat == 'average branch length \u2264 10000') / nrow(temp) # ratio of loci with average branch length > 10000 / all loci
# boxplot of mutation rate by number of alleles
colors = c('#fc8d59','#3288bd')
png(file='bias by small branch.png', width=3000, height=2000, res=300)
ggplot(data = temp[temp$comAll < 10,]) +
  theme_classic() +
  geom_boxplot(aes(x=branchCat, y=log10(rateNormAfter), color=branchCat)) + 
  xlab('Number of alleles by composition') + ylab('Mutation rate (log10)') +
  facet_wrap(~comAll, strip.position='bottom', scales='free_x', nrow=2) + # set multiple group labels on x axis
  scale_color_manual(name='', values=colors) + # manually set color scheme for dot color
  theme(text = element_text(size = 30),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.title = element_blank(), legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()

hist(temp[temp$comAll < 10,]$aveDist)
sum(temp$aveDist <= 10000) / nrow(temp)
# number of unique trees and number of trees with small average branch
temp$treeKey = paste0(temp$chr,'_',temp$treeS,'_',temp$treeE)
length(unique(temp$treeKey[temp$aveDist <= 10000])) / length(unique(temp$treeKey))

######################################## association ########################################
######################################################################
# density plot:  rate by GC, overlay by disease loci                 #
######################################################################
# plot rateNormAfter vs. GC, disease overlap
png(file='mutRate vs. GC content.overlay.png', width=2500, height=2500, res=300)
#ggplot(data = rbind(rate[rate$spanned != 'noTree' & rate$period != 1, c(1:22,26:45)],
ggplot(data = rbind(rate[rate$spanned != 'noTree', c(1:22,26:45)],
                    rateBio[rate$spanned != 'noTree', c(1:22,26:45)])) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=allGC), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('GC content') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data=rateBio[rateBio$spanned != 'noTree',],
             aes(x=log10(rateNormAfter), y=allGC, pch=period6, color=coding), cex=3, stroke=1) +
  scale_shape_manual(name = 'period size', values = c(1,2,3,4,5)) +
  scale_color_manual(name = 'coding', values = c('pink','red'), labels = c('non-coding','coding')) +
  guides(fill = guide_colourbar(title.position='top', direction='horizontal', barwidth=unit(5, 'cm')),
         shape = guide_legend(ncol=3, title.position='top'),
         color = guide_legend(nrow=2, title.position='top')) +
  theme(text = element_text(size = 30), legend.title = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'top', legend.direction = 'horizontal',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()

# plot rateNormAfter vs. GC, disease transparent overlay
png(file='mutRate vs. GC content.transparent.png', width=3000, height=2000, res=300)
#ggplot(data = rbind(rate[rate$spanned != 'noTree' & rate$period != 1, c(1:22,26:45)],
ggplot(data = rbind(rate[rate$spanned != 'noTree', c(1:22,26:45)],
                    rateBio[rate$spanned != 'noTree', c(1:22,26:45)])) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=allGC), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('GC content') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data=rateBio[rateBio$spanned != 'noTree',],
             aes(x=log10(rateNormAfter), y=allGC, pch=period6, color=coding), cex=2, stroke=1, alpha=0) +
  scale_shape_manual(name = 'period size', values = c(1,2,3,4,5)) +
  scale_color_manual(name = 'coding', values = c('pink','red'), labels = c('non-coding','coding')) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()

# plot not normalized rate vs. GC
png(file='mutRate vs. GC content.overlay.noNorm.png', width=3000, height=2000, res=300)
#ggplot(data = rbind(rate[rate$spanned != 'noTree' & rate$period != 1, c(1:22,26:45)],
ggplot(data = rbind(rate[rate$spanned != 'noTree', c(1:22,26:45)],
                    rateBio[rate$spanned != 'noTree', c(1:22,26:45)])) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rate), y=allGC), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('GC content') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data=rateBio[rateBio$spanned != 'noTree',],
             aes(x=log10(rate), y=allGC, pch=period6, color=coding), cex=2, stroke=1) +
  scale_shape_manual(name = 'period size', values = c(1,2,3,4,5)) +
  scale_color_manual(name = 'coding', values = c('pink','red'), labels = c('non-coding','coding')) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()

# plot rateNormBefore vs. rateNormAfter
png(file='mutRateNormBefore vs. mutRateNormAfter.png', width=3000, height=2000, res=300)
#ggplot(data = rbind(rate[rate$spanned != 'noTree' & rate$period != 1, c(1:22,26:45)],
ggplot(data = rbind(rate[rate$spanned != 'noTree', c(1:22,26:45)],
                    rateBio[rate$spanned != 'noTree', c(1:22,26:45)])) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormBefore), y=log10(rateNormAfter)), bins=50) + 
  xlab('Normalize before average (log10)') + ylab('Normalize after average (log10)') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data=rateBio[rateBio$spanned != 'noTree',],
             aes(x=log10(rateNormBefore), y=log10(rateNormAfter), pch=period6, color=coding), cex=2, stroke=1) +
  scale_shape_manual(name = 'period size', values = c(1,2,3,4,5)) +
  scale_color_manual(name = 'coding', values = c('pink','red'), labels = c('non-coding','coding')) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()


# test for GC
#temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0 & rate$period != 1,]
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp1 = temp$rateNormAfter[temp$allGC >= 0.5]
temp2 = temp$rateNormAfter[temp$allGC < 0.5]
t.test(temp1, temp2); mean(temp1) / mean(temp2) # t test
model <- lm(temp$rateNormAfter ~ temp$allGC); summary(model) # linear regression

# contrast coding bio loci vs. non-coding bio loci
mean(rateBio$rateNormAfter[rateBio$coding == '.']) / mean(rateBio$rateNormAfter[rateBio$coding != '.'])
t.test(rateBio$rateNormAfter[rateBio$coding == '.'], rateBio$rateNormAfter[rateBio$coding != '.'])
# check GC richness
hist(catalogBio$allGC[catalogBio$coding != '.'])
table(catalogBio$period[catalogBio$coding != '.'])

######################################################################
# density plot:  rate by consensus entropy                           #
######################################################################
#ggplot(data = rbind(rate[rate$spanned != 'noTree' & rate$period != 1, c(1:22,26:45)],
ggplot(data = rbind(rate[rate$spanned != 'noTree', c(1:22,26:45)],
                    rateBio[rate$spanned != 'noTree', c(1:22,26:45)])) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=entropy), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('entropy') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data=rateBio[rateBio$spanned != 'noTree',],
             aes(x=log10(rateNormAfter), y=entropy, pch=period6, color=coding), cex=2, stroke=1) +
  scale_shape_manual(name = 'period size', values = c(1,2,3,4,5)) +
  scale_color_manual(name = 'coding', values = c('pink','red'), labels = c('non-coding','coding')) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')

#temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0 & rate$period != 1,]
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp1 = temp$rateNormAfter[temp$entropy >= 0.6 & temp$entropy <= 0.7]
temp2 = temp$rateNormAfter[temp$entropy < 0.6 | temp$entropy > 0.7]
t.test(temp1, temp2); mean(temp1) / mean(temp2) # t test
model <- lm(temp$rateNormAfter ~ temp$entropy); summary(model) # linear regression

######################################################################
# density plot:  rate by consensus length                            #
######################################################################
#ggplot(data = rate[rate$spanned != 'noTree' & rate$period != 1 & rate$period < 200,]) +
ggplot(data = rate[rate$spanned != 'noTree' & rate$period < 200,]) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=period), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('motif period') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')

#temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0 & rate$period != 1,]
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp1 = temp$rateNormAfter[temp$period > 6]
temp2 = temp$rateNormAfter[temp$period <= 6]
t.test(temp1, temp2); mean(temp1) / mean(temp2) # t test
model <- lm(temp$rateNormAfter ~ temp$period); summary(model) # linear regression

######################################################################
# density plot:  rate by locus purity                                #
######################################################################
ggplot(data = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=purityEff), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('purity (effMotifs)') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')

#temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0 & rate$purityEff != 1,]
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
median(temp$purityEff)
temp1 = temp$rateNormAfter[temp$purityEff <= 0.5]
temp2 = temp$rateNormAfter[temp$purityEff > 0.5]
t.test(temp1, temp2); mean(temp1) / mean(temp2) # t test
model <- lm(temp$rateNormAfter ~ temp$purityEff); summary(model) # linear regression


######################################################################
# density plot:  rate by locus length                                #
######################################################################
temp = rate; temp$aveLen = temp$totalLen / temp$annotated
png(file='mutRate vs. length.loglog.png', width=3000, height=2000, res=300)
#ggplot(data = temp[temp$spanned != 'noTree' & temp$period != 1 & temp$aveLen < 1000,]) +
ggplot(data = temp[temp$spanned != 'noTree' & rate$rateNormAfter > 0,]) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=log10(aveLen)), bins=30) + 
  geom_hline(yintercept = log10(5), color='red', lty='dashed') +
  geom_hline(yintercept = log10(30), color='red', lty='dashed') +
  xlab('Normalized mutation rate (log10)') + ylab('Locus length by motif (log10)') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')
dev.off()

#temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0 & rate$period != 1,]
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp$aveLen = temp$totalLen / temp$annotated
temp1 = temp$rateNormAfter[temp$aveLen > 200]
temp2 = temp$rateNormAfter[temp$aveLen < 200]
t.test(temp1, temp2); mean(temp1) / mean(temp2) # t test
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp$aveLen = temp$totalLen / temp$annotated
temp = temp[temp$aveLen > 5 & temp$aveLen < 30,]
model <- lm(log10(temp$rateNormAfter) ~ log10(temp$aveLen));summary(model) # linear regression

# before normalization
temp = rate; temp$aveLen = temp$totalLen / temp$annotated
ggplot(data = temp[temp$spanned != 'noTree' & temp$period != 1 & temp$aveLen < 1000,]) +
  theme_classic() +
  #geom_bin2d(aes(x=rateNormAfter, y=size), bins=50) + 
  geom_bin2d(aes(x=log10(rate), y=aveLen), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('Locus length') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')

temp = temp[temp$spanned != 'noTree' & temp$period != 1 & temp$rateNormAfter > 0,]
temp$rateLog10 = log10(temp$rate)
model <- lm(temp$rateLog10 ~ temp$aveLen)
summary(model)







######################################################################
# check AT loci                                                      #
######################################################################
AT5Bio = rateBio[rateBio$allGC<0.05 & rateBio$period==5,
                 c('chr','start','end','size','consensus','consGC','allGC','rateNormAfter','stdLenAll','numOutlierAll',
                   'segDup','coding','exon','transcribed')]
write.table(AT5Bio, file='ATTTT.bio.tsv', quote=FALSE, row.names=FALSE, sep='\t')

AT5genome = rate[rate$allGC<0.05 & rate$period==5
                 & !is.na(rate$rateNormAfter) & rate$rateNormAfter>0.00001
                 & rate$transcribed == 'transcribed'
                 & rate$size > 40 & rate$size < 120,
                 c('chr','start','end','size','consensus','consGC','allGC','rateNormAfter','stdLenAll','numOutlierAll',
                   'segDup','coding','exon','transcribed')]
write.table(AT5genome, file='ATTTT.tsv', quote=FALSE, row.names=FALSE, sep='\t')

# plot a histogram to show the outlier
outlierExample = read.csv('selectedforhistogram.annoLen.tsv', sep='\t', header=F)
hist(as.numeric(outlierExample[1,4:ncol(outlierExample)]))
hist(as.numeric(outlierExample[2,4:ncol(outlierExample)]))

# total no. of outliers in the genome (constraining sd)
nrow(rate[rate$spanned != 'noTree',]) # number of successful loci with computed mutation ratesum(rate$stdLenAll>10)
nrow(rate[rate$stdLenAll>5 & rate$spanned != 'noTree',])
nrow(rate[rate$stdLenAll>5 & rate$numOutlierAll > 0 & rate$spanned != 'noTree',])
nrow(rate[rate$stdLenAll>5 & rate$numOutlierAll > 0 & rate$spanned != 'noTree' & rate$coding == 'cds',])

hyperexpansion = rate[rate$stdLenAll>5 & rate$numOutlierAll > 0 & rate$spanned != 'noTree' & rate$coding == 'cds',]
write.table(hyperexpansion[,2:ncol(hyperexpansion)], file='hyper-expansion.tsv', quote=FALSE, row.names=FALSE, sep='\t')

######################################################################
# rate/hotspot by coding/non-coding                                  #
######################################################################
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
fisher.test(table(temp$hotspot, temp$coding))
t.test(rate$rateNormAfter[rate$coding=='cds'],
       rate$rateNormAfter[rate$coding=='.'])
mean(rate$rateNormAfter[rate$coding=='.']) / mean(rate$rateNormAfter[rate$coding=='cds'])
######################################################################
# rate/hotspot by HAR                                                #
######################################################################
# size of HARs
(853533 / 3088269832 ) * 4000000 # so, the number of loci in HARs is proportional to the size of HARs in the genome
853533 / 3168
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
mean(temp$rateNormAfter[temp$HAR!='.']) / mean(temp$rateNormAfter[temp$HAR=='.']) # rate in HARs / rate not in HARs
ratios = c()
# permutation of HAR tags
set.seed(123)
for (i in 1:100) {
  newHAR = sample(temp$HAR)
  ratios = c( ratios, mean(rate$rateNormAfter[newHAR!='.'])/mean(rate$rateNormAfter[newHAR=='.']) )
}
ratios
temp$inHAR = temp$HAR != '.'
fisher.test(table(temp$hotspot, temp$inHAR))
t.test(rate$rateNormAfter[rate$HAR!='.'],
       rate$rateNormAfter[rate$HAR=='.'])
mean(rate$rateNormAfter[rate$HAR!='.']) / mean(rate$rateNormAfter[rate$HAR=='.'])
######################################################################
# rate/hotspot by recombination hotspot                              #
######################################################################
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp$inRecomb = temp$dist2Recom == 0
fisher.test(table(temp$hotspot, temp$inRecomb))
######################################################################
# rate/hotspot by segdup                                             #
######################################################################
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,]
temp$inSegdup = temp$segDupIdentity == -1
fisher.test(table(temp$hotspot, temp$inSegdup)) # fisher exact on hotspot vs. segdup, by TRs not windows
temp1 = temp$rateNormAfter[temp$segDupIdentity != -1]
temp2 = temp$rateNormAfter[temp$segDupIdentity == -1]
t.test(temp1, temp2); mean(temp1) / mean(temp2) # t.test on in/out segdup

# permutation test
pvalues = read.csv('p-value', sep='\t', header=F)
hist(pvalues$V1, breaks=seq(0,1,0.05)) # fisher exact p-values
hist(pvalues$V3); sum(pvalues$V3 < 87) # no. of hotspots

# check for high identity/CNV segdups
plotData = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0,
                c('key','rateNormAfter','hotspot','segDupIdentity','segDupCNV')]

plotData$identity = 'high identity\nsegdup'
plotData$identity[plotData$segDupIdentity <= 99] = 'low identity\nsegdup'
plotData$identity[plotData$segDupIdentity == -1] = 'not in\nsegdup'
plotData$CNV = 'high CNV\nsegdup'
plotData$CNV[plotData$segDupCNV <= 0.1] = 'low CNV\nsegdup'
plotData$CNV[plotData$segDupCNV == -1] = 'not in\nsegdup'
plotData$segDupAll = 'NULL'
plotData$segDupAll[plotData$segDupIdentity <= 99 & plotData$segDupCNV <= 0.1] = 'low identity\nlow CNV'
plotData$segDupAll[plotData$segDupIdentity <= 99 & plotData$segDupCNV > 0.1] = 'low identity\nhigh CNV'
plotData$segDupAll[plotData$segDupIdentity > 99  & plotData$segDupCNV <= 0.1] = 'high identity\nlow CNV'
plotData$segDupAll[plotData$segDupIdentity > 99  & plotData$segDupCNV > 0.1] = 'high identity\nhigh CNV'
plotData$segDupAll[plotData$segDupIdentity == -1] = 'not in\nsegdup'

t.test(plotData$rateNormAfter[plotData$identity == 'not in\nsegdup'],
       plotData$rateNormAfter[plotData$identity == 'low identity\nsegdup'])
# fold change of mutation rate of low identity segdup to outside segdup
mean(plotData$rateNormAfter[plotData$identity == 'low identity\nsegdup']) /
  mean(plotData$rateNormAfter[plotData$identity == 'not in\nsegdup'])
t.test(plotData$rateNormAfter[plotData$CNV == 'not in\nsegdup'],
       plotData$rateNormAfter[plotData$CNV == 'low CNV\nsegdup'])
mean(plotData$rateNormAfter[plotData$CNV == 'low CNV\nsegdup']) /
  mean(plotData$rateNormAfter[plotData$CNV == 'not in\nsegdup'])

#colors = c('#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494')

table(plotData$identity, plotData$hotspot)
prop.table(table(plotData$identity, plotData$hotspot))
png(file='mutRate_by_segDup_identity.png', width=3000, height=1000, res=300)
ggplot(plotData) +
  theme_classic() +
  geom_boxplot(aes(x=identity,y=log10(rateNormAfter),color=identity)) +
  #scale_color_manual(name='', values=colors) + # manually set color scheme for bar fill
  theme(text = element_text(size = 15), legend.position='none') +
  xlab('') + ylab('Mutation rate (log10)') +
  facet_wrap(~ hotspot, nrow=1)
dev.off()

table(plotData$CNV, plotData$hotspot)
prop.table(table(plotData$identity, plotData$hotspot))
png(file='mutRate_by_segDup_CNV.png', width=3000, height=1000, res=300)
ggplot(plotData) +
  theme_classic() +
  geom_boxplot(aes(x=CNV,y=log10(rateNormAfter),color=CNV)) +
  #scale_color_manual(name='', values=colors) + # manually set color scheme for bar fill
  theme(text = element_text(size = 15), legend.position='none') +
  xlab('') + ylab('Mutation rate (log10)') +
  facet_wrap(~ hotspot, nrow=1)
dev.off()

png(file='mutRate_by_segDup_all.png', width=4000, height=2000, res=300)
ggplot(plotData) +
  theme_classic() +
  geom_boxplot(aes(x=segDupAll,y=log10(rateNormAfter),color=segDupAll)) +
  #scale_color_manual(name='', values=colors) + # manually set color scheme for bar fill
  theme(text = element_text(size = 15), legend.position='none') +
  xlab('') + ylab('Mutation rate (log10)') +
  facet_wrap(~ hotspot, nrow=1)
dev.off()


######################################################################
# rate/hotspot by transposable elements                              #
######################################################################
temp = rate[rate$spanned != 'noTree' & rate$period != 1 ,]
temp$rateLog10 = temp$rateNormAfter; temp$rateLog10[temp$rateNormAfter > 0] = log10(temp$rateNormAfter[temp$rateNormAfter > 0])
temp$op = temp$op1 > 0

 # fisher test
t.test(temp$rateNormAfter[temp$op], temp$rateNormAfter[!temp$op])
model <- glm(op ~ rateNormAfter, family = binomial(link = "logit"), data = temp); summary(model); exp(coef(model)) # logistic model
model <- lm(temp$rateNormAfter ~ temp$op2); summary(model) # linear model







######################################## block analysis ########################################
####################
# block by recombination dist
####################
# summary statistics
# autosomal TRs with motif sizes greater than two nucleotides
nrow(catalog[catalog$period > 2 & catalog$chr != 'chrX' & catalog$chr != 'chrY',])
# autosomal TRs with motif sizes greater than two nucleotides and have at least 20 annotated genomes
nrow(block[block$period > 2 & block$chr != 'chrX' & block$chr != 'chrY',])
# block positive loci
nrow(block[block$period > 2 & block$chr != 'chrX' & block$chr != 'chrY' & block$aveBlock > -1,])
29964 / 2221932

# block yes/no vs. distance to recombination hotspot (t-test)
temp = rate[rate$period > 2 & rate$telomere == 'telo-no',]
t.test(temp$dist2Recom[temp$aveBlockCat == 'block-no'], temp$dist2Recom[temp$aveBlockCat != 'block-no'])
389428.0 / 341231.5

# test variability vs. distance to recombination hotspot
temp = rate[rate$period > 2 & rate$telomere == 'telo-no',]
temp$lenAve = temp$lenAll / temp$annotated
temp$comAve = temp$comAll / temp$annotated
t.test(temp$dist2Recom[temp$lenAve <= 0.1 & temp$lenAve <= 6], temp$dist2Recom[temp$lenAve > 0.1 & temp$lenAve <= 6])
plot(log10(temp$dist2Recom), temp$lenAve)
model <- lm(temp$lenAve ~ temp$dist2Recom)
summary(model)

temp = parsimony[parsimony$start == 155188482,]

##### test of block vs. variability
temp = rate[rate$period > 2 & rate$rateNormAfter > 0,]
temp$lenAve = temp$lenAll / temp$annotated
temp$comAve = temp$comAll / temp$annotated

# average block vs. heterozygosity by len (linear model)
plot(temp$heteLenAll[temp$aveBlockCat != 'block no'], temp$aveBlock[temp$aveBlockCat != 'block no'])
model <- lm(temp$aveBlock[temp$aveBlockCat != 'block no'] ~ temp$heteLenAll[temp$aveBlockCat != 'block no'])
summary(model)
# average block vs. heterozygosity by com (linear model)
plot(temp$heteComAll[temp$aveBlockCat != 'block no'], temp$aveBlock[temp$aveBlockCat != 'block no'])
model <- lm(temp$aveBlock[temp$aveBlockCat != 'block no'] ~ temp$heteComAll[temp$aveBlockCat != 'block no'])
summary(model)

# block yes/no vs. average allele by len (t-test)
t.test(temp$lenAve[temp$aveBlockCat != 'block-no'], temp$lenAve[temp$aveBlockCat == 'block-no'])
0.023557016 / 0.004639978
# average block vs. average allele by len (linear model)
plot(temp$aveBlock[temp$aveBlockCat != 'block no'], temp$lenAve[temp$aveBlockCat != 'block no'])
model <- lm(temp$lenAve[temp$aveBlockCat != 'block no'] ~ temp$aveBlock[temp$aveBlockCat != 'block no'])
summary(model)

# average block vs. average allele by com (linear model)
plot(temp$aveBlock[temp$aveBlockCat != 'block no'], temp$comAve[temp$aveBlockCat != 'block no'])
model <- lm(temp$comAve[temp$aveBlockCat != 'block no'] ~ temp$aveBlock[temp$aveBlockCat != 'block no'])
summary(model)


# linear model
model <- lm(temp$aveBlock[temp$dist2Recom != 0] ~ temp$dist2Recom[temp$dist2Recom != 0])
summary(model)

ggplot(temp) +
  theme_classic() +
  geom_boxplot(aes(x=aveBlockCat,y=log10(dist2Recom),color=aveBlockCat)) +
  #scale_color_manual(name='', values=colors) + # manually set color scheme for bar fill
  theme(text = element_text(size = 15), legend.position='none') +
  xlab('') + ylab('Dist to recombination hotspot (log10)')

####################
# parsimony analysis
####################
# number of block positive VNTR loci fully spanned
nrow(rate[rate$period > 6 & rate$spanned == 'full' & rate$numPattern > 0,])
# number of block positive VNTR loci fully spanned, having > 1 pattern
nrow(rate[rate$period > 6 & rate$spanned == 'full' & rate$numPattern > 1,])
3129 / 4653
# number of block positive VNTR loci with recurrent mutations
nrow(parsimony[parsimony$totalRecurrence > 0,])
3351 / 4643

parsimony$aveRecurrence = parsimony$totalRecurrence / parsimony$annotated
parsimony$normRecurrence = parsimony$aveRecurrence / (parsimony$totalLen / parsimony$annotated )
# normalized average recurrence vs. ???
plot(parsimony$normRecurrence, parsimony$aveBlock)
model <- lm(parsimony$normRecurrence ~ parsimony$aveBlock)
summary(model)
plot(parsimony$normRecurrence, parsimony$numPattern)
plot(parsimony$normRecurrence, parsimony$rateNormAfter)
model <- lm(parsimony$normRecurrence ~ parsimony$rateNormAfter)
summary(model)

temp = parsimony[parsimony$telomere == 'telo-no',]
plot(log10(temp$dist2Recom), temp$normRecurrence)
model <- lm(temp$normRecurrence ~ temp$dist2Recom)
summary(model)
t.test(temp$normRecurrence[temp$dist2Recom < 10000], temp$normRecurrence[temp$dist2Recom >= 10000])





####################
# everything on ideogram
####################
library(karyoploteR)

allChrs = c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
            'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
            'chr20','chr21','chr22')


# read mutation rate hotspots and segdup and recombination regions
hotspot = read.csv('all.sorted.hotspot.tsv', sep='\t', header=F)
colnames(hotspot) = c('chr','start','end','startSNP','endSNP','type')
# read segDup bed
segDups = read.csv('segDup.sort.merge.bed', sep='\t', header=F)
colnames(segDups)[1:3] = c('chr','start','end')
# read recomb hotspots
recomb = read.csv('decode_recomb_hotspot_both.tsv.hg38.sorted.merged.bed', sep='\t', header=F)
colnames(recomb)[1:3] = c('chr','start','end')
# read HARs
HARs = read.csv('GSE180714_HARs.bed', sep='\t', header=T)

##### plot for all chroms
png(file='rate.ideogram.png', width=6000, height=6000, res=300) # all plot
kp = plotKaryotype(genome="hg38", chromosomes=allChrs)
# recombination hotspots
kpPlotDensity(kp, toGRanges(recomb), r0=0.7, r1=1.0, window.size=200000, col="#FFAACC")
# segdup density
kpPlotDensity(kp, toGRanges(segDups), r0=0.4, r1=0.7, window.size=200000, col="#FFAACC")
# TR mutation rate hotspots
kpPlotDensity(kp, toGRanges(hotspot), r0=0.0, r1=0.2, window.size=1000, col="#ccece6")
# rate
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0, c('chr','start','end','rateNormAfter')]
colnames(temp) = c('chr','start,','end','y')
kpPoints(kp, data=toGRanges(temp), y=temp[,4], r0=0.0, r1=0.3, ymin=0.00000000001, ymax=0.001, cex=0.1)
dev.off()

##### plot for chr1 and chr2
png(file='rate.ideogram.chr1-2.png', width=3000, height=1500, res=300) # chr1-3
pp <- getDefaultPlotParams(plot.type = 1)
pp$leftmargin = 0.12; pp$rightmargin = 0.01; pp$topmargin = 20; pp$bottommargin = 20
kp = plotKaryotype(genome="hg38", chromosomes=allChrs[c(1,2)], plot.type=1, cex=2, plot.params=pp)
# recombination hotspots
kpPlotDensity(kp, toGRanges(recomb), r0=0.7, r1=1.0, window.size=200000, col="#FFAACC")
kpAddLabels(kp, labels="Recombination", r0=0.6, r1=1.0, cex=1)
# segdup density
kpPlotDensity(kp, toGRanges(segDups), r0=0.4, r1=0.7, window.size=200000, col="#FFAACC")
kpAddLabels(kp, labels="Segmental\nduplications", r0=0.2, r1=0.7, cex=1)
# TR mutation rate hotspots
kpPlotDensity(kp, toGRanges(hotspot), r0=0.0, r1=0.2, window.size=1000, col="#ccece6")
# rate
temp = rate[rate$spanned != 'noTree' & rate$rateNormAfter > 0, c('chr','start','end','rateNormAfter')]
colnames(temp) = c('chr','start,','end','y')
kpPoints(kp, data=toGRanges(temp), y=temp[,4], r0=0.0, r1=0.3, ymin=0.00000000001, ymax=0.001, cex=0.1)
kpAddLabels(kp, labels="Mutation rate", r0=-0.1, r1=0.3, cex=1)
#temp = rate[rate$period > 4, c('chr','start','end','aveBlock')]; colnames(temp) = c('chr','start,','end','y')
#kpPoints(kp, data=toGRanges(temp), y=temp[,4], r0=0.0, r1=0.2, ymin=-1, ymax=50, cex=0.1) # average block
# treeEd and alleleDist
#binsEd = toGRanges(trees[trees$aveEd < 2, c(1,2,3,6)])
#kpPoints(kp, data=binsEd, y=trees[trees$aveEd < 2,6], r0=0.21, r1=0.4, ymin=0, ymax=2, cex=0.1, col='red')
#binsDist = toGRanges(trees[trees$aveDist < 50000, c(1,2,3,7)])
#kpPoints(kp, data=binsDist, y=trees[trees$aveDist < 50000,7], r0=0.0, r1=0.19, ymin=0, ymax=50000, cex=0.1, col='blue')
dev.off()





####################
# compare to Gymrek
####################
Gymrek = read.csv('intersectGymrekbyVamos.tsv', sep='\t', header=T)
Gymrek$startDiff = abs(Gymrek$start - Gymrek$startV)
Gymrek$endDiff = abs(Gymrek$end - Gymrek$endV)
Gymrek$rateVlog10 = log10(Gymrek$rateV)
Gymrek$rateVlog10[Gymrek$rateV==0] = -9

GymrekFilter = Gymrek[Gymrek$period != 1 & Gymrek$startDiff < 2 & Gymrek$endDiff < 2, ]
GymrekFilter = na.omit(GymrekFilter)

# visual guide lines
y = seq(-7.2,-2.8,0.05)
x1 = (y-0.65) / 1.1 - 0.02
x2 = (y-0.65) / 1.1 + 0.02
x3 = (y-0.65) / 1.1
lineData = data.frame(x1=x1, y1=y, x2=x2, y2=y, x3=x3, y3=y)

png(file='V vs. G.png', width=2600, height=2600, res=300)
ggplot(data = GymrekFilter) +
  theme_classic() +
  geom_bin2d(aes(x=rateVlog10, y=rateG), bins=50) + 
  xlab('Mutation rate before normalization (log10)') + ylab('SRS-length-rate (log10)') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_y_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data=lineData, aes(x=x3,y=y3), cex=0.5, color='white') +
  annotate(x=-1.8, y=-1, label='slope=1.1', vjust=2, geom="label", cex=8) + # fit line text
  guides(fill = guide_colourbar(title.position='top', direction='horizontal', barwidth=unit(5, 'cm'))) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()

sum( GymrekFilter$rateG < -7.5 & GymrekFilter$rateVlog10 == -9 ) # bottom left corner: 76553
sum( GymrekFilter$rateG < -7.5 & GymrekFilter$rateVlog10 != -9 ) # bottom line: 23479
23479/nrow(GymrekFilter)
sum( GymrekFilter$rateVlog10 == -9 ) # all left line: 87614
sum( GymrekFilter$rateG > -7.5 & GymrekFilter$rateVlog10 == -9 ) # left line without the corner: 11061
11061/nrow(GymrekFilter)
sum( GymrekFilter$rateG > -2.5 & GymrekFilter$rateVlog10 == -9 ) # left line top part: 467
sum( GymrekFilter$rateG > -2.5 & GymrekFilter$rateVlog10 != -9 ) # top region: 1324
1324/nrow(GymrekFilter)

sum( GymrekFilter$rateG > -2 & GymrekFilter$rateVlog10 != -9 
     & GymrekFilter$rateG - GymrekFilter$rateVlog10 > 2 ) # upper region: 779
sum( GymrekFilter$rateG > -7.5 & GymrekFilter$rateG < -2
     & GymrekFilter$rateVlog10 != -9 ) # central region: 126600
sum( GymrekFilter$rateVlog10 != -9 ) # right region: 150932

model <- lm(GymrekFilter$rateG[GymrekFilter$rateG > -7.5 & GymrekFilter$rateV > 0]
            ~ GymrekFilter$rateVlog10[GymrekFilter$rateG > -7.5 & GymrekFilter$rateV > 0])
summary(model)



# check for the main divergent ones
kp = plotKaryotype(genome="hg38", chromosomes=allChrs[c(1,2,3)])
kpPlotDensity(kp, toGRanges(Gymrek[Gymrek$rateG > -2,][,1:3]), r0=0.7, r1=1.0, window.size=2000, col="#FFAACC") # recombination hotspots

# number of VNTR that are only covered in our study, not in Gymmeric
(364235 - 293206) / 364235


####################
# compare to Evan
####################
Evan = read.csv('intersectEvanbyVamos.tsv', sep='\t', header=T)
Evan$startDiff = abs(Evan$start - Evan$startV)
Evan$endDiff = abs(Evan$end - Evan$endV)
Evan$rateVlog10 = log10(Evan$rateNormAfter)
Evan$rateVlog10[Evan$rateNormAfter==0] = -9
Evan = na.omit(Evan)

temp = rate[rate$spanned != 'noTree', c(1:22,26:45)] #rate$spanned != 'noTree' & rate$period != 1, c(1:22,26:45)]
temp$rateVlog10 = log10(temp$rateNormAfter)
temp$rateVlog10[temp$rateNormAfter==0] = -9

# test of mapped loci vs. genomic
t.test(Evan$rateNormAfter, temp$rateNormAfter)
t.test(Evan[Evan$rateNormAfter > 0,]$rateNormAfter, temp[temp$rateNormAfter > 0,]$rateNormAfter)

#temp = rbind(Evan[,c('chr','start','end','rateVlog10')],temp[,c('chr','start','end','rateVlog10')])

png(file='V vs. E.png', width=1200, height=2000, res=300)
colors <- c('de novo'='#fdb462','genomic'='#8dd3c7')
ggplot() +
  theme_classic() +
  geom_histogram(data=Evan[Evan$rateVlog10 > -9,], aes(x=rateVlog10, y=after_stat(density), fill='de novo'),
                 binwidth=0.5, color='black', alpha=0.5) +
  geom_histogram(data=temp[temp$rateVlog10 > -9,], aes(x=rateVlog10, y=after_stat(density), fill='genomic'),
                 binwidth=0.5, color='black', alpha=0.5) +
  scale_fill_manual(name = '', values = colors, labels = c('de novo','genomic')) +
  xlab('Mutation rate (log10)') + ylab('Proportion') +
  #guides(fill = guide_legend(nrow = 2)) +
  theme(text = element_text(size = 25),
        legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()

# number of loci with mutation rate > -5
sum(rate$rateNormAfter> 0.00001)
sum(rate$rateNormAfter> 0.00001) / 65.3 # minimum no. of family-based samples needed to observe all rate>0.00001 loci










####################
# compare variability (heterozygosity) with mutation rate
####################
# for all loci
t.test(rate$rateNormAfter[rate$segDup != 'not in segDup'],
       rate$rateNormAfter[rate$segDup == 'not in segDup'])
# for only variable loci
t.test(rate$rateNormAfter[rate$segDup != 'not in segDup' & rate$rateNormAfter > 0],
       rate$rateNormAfter[rate$segDup == 'not in segDup' & rate$rateNormAfter > 0])


####################
# compare variability (heterozygosity) with mutation rate
####################
# data filtering and variable creation
alleles = read.csv('../2025-03-02_assembly_stats/alleles/samples_rmDup_rmKids_rmPoor_GRCh38_q-0.1.noHPRC2.annoStr.tsv', sep='\t', header=T)

# filter chrX,chrY, and homopolymers
alleles = alleles[alleles$chr != 'chrX',]; alleles = alleles[alleles$chr != 'chrY',]
#alleles = alleles[alleles$period != 1,]
alleles$key = paste0(alleles$chr,'_',alleles$start,'_',alleles$end)

overall = alleles[,c('homoLenAll','heteLenAll','homoComAll','heteComAll','lenAll','comAll','lenSDAll','key')]
colnames(overall)[1:7] = c('homoLen','heteLen','homoCom','heteCom','len','com','lenSD')

plotData = merge(rate, overall, by='key')

png(file='mutRate vs. lenSD.png', width=3000, height=2000, res=300)
ggplot(data = plotData[plotData$spanned != 'noTree' & plotData$period != 1 & plotData$lenSD < 50,]) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=lenSD), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('Length standard deviation') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 15), legend.position = 'right')
dev.off()

temp = plotData[plotData$spanned != 'noTree' & plotData$period != 1 & plotData$rateNormAfter > 0 & plotData$lenSD < 10,]
temp$rateLog10 = log10(temp$rateNormAfter)
model <- lm(temp$rateNormAfter ~ temp$lenSD)
model <- lm(temp$rateLog10 ~ temp$lenSD)
summary(model)






####################
# rate vs. average branch length
####################
rate$nPairs = rate$annotated * (rate$annotated-1) / 2
rate$avePairDist = rate$totalDist / rate$nPairs
rate$avePairEd = rate$totalEd / rate$nPairs
rate$avePairEdNorm = rate$avePairEd / (rate$totalLen / rate$annotated)

rateBio$nPairs = rateBio$annotated * (rateBio$annotated-1) / 2
rateBio$avePairDist = rateBio$totalDist / rateBio$nPairs
rateBio$avePairEd = rateBio$totalEd / rateBio$nPairs
rateBio$avePairEdNorm = rateBio$avePairEd / (rateBio$totalLen / rateBio$annotated)

# downsample for plot
set.seed(123) # for reproducibility
temp <- rate[sample(nrow(rate[rate$rateNormAfter>0,]), size = 100000), ]

plot(temp$snpDensity, log10(temp$rateNormAfter), cex=0.1)

library(ggplot2)

ggplot(data = rate[rate$spanned != 'noTree' & rate$period != 1,]) +
  theme_classic() +
  geom_bin2d(aes(y=log10(rateNormAfter), x=avePairDist), bins=50) + 
  ylab('Mutation rate (log10)') + xlab('Average pairwise branch length') +
  scale_y_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10")






# check hotspot map difference vs. max segdup identity
temp = hotspot[hotspot$idenMax == 0,]
plot(temp$idenMax, temp$meanMapDiff)
temp = hotspot[hotspot$idenMax > 0,]
plot(temp$idenMax, temp$meanMapDiff)



####################
# compare variability (heterozygosity) with mutation rate
####################
# data filtering and variable creation
alleles = read.csv('../2025-03-02_assembly_stats/alleles/samples_rmDup_rmKids_rmPoor_GRCh38_q-0.1.noHPRC2.annoStr.tsv', sep='\t', header=T)

# filter chrX,chrY, and homopolymers
alleles = alleles[alleles$chr != 'chrX',]; alleles = alleles[alleles$chr != 'chrY',]
#alleles = alleles[alleles$period != 1,]
alleles$key = paste0(alleles$chr,'_',alleles$start,'_',alleles$end)

overall = alleles[,c('homoLenAll','heteLenAll','homoComAll','heteComAll','lenAll','comAll','lenSDAll','key')]
colnames(overall)[1:7] = c('homoLen','heteLen','homoCom','heteCom','len','com','lenSD')

plotData = merge(rate, overall, by='key')

png(file='mutRate vs. lenSD.png', width=2000, height=2000, res=300)
ggplot(data = plotData[plotData$spanned != 'noTree' & plotData$period != 1 & plotData$lenSD < 50,]) +
  theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=lenSD), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('Length standard deviation') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 15), legend.position = 'top')
dev.off()



AT5genome = plotData[plotData$allGC<0.05 & plotData$period==5
                 & !is.na(plotData$rateNormAfter) & plotData$rateNormAfter>0.00001
                 & plotData$transcribed == 'transcribed'
                 & plotData$size > 40 & plotData$size < 120
                 & plotData$spanned != 'noTree',
                 c('chr','start','end','period','period6','length','spanned','rateNormAfter','consensus','allGC','lenSD')]







#### PCA playground
temp = rate[,c('chr','start','end','period','consGC','allGC','totalLen','stdLenAll',
               'op1','op2','segDup','coding','exon','transcribed','dist2Recom','rate')]
# engineer categorical variables into numerical
temp$segDupNum = 0; temp$segDupNum[temp$segDup == 'part in segDup'] = 1; temp$segDupNum[temp$segDup == 'in segDup'] = 2
temp$codingNum = 0; temp$codingNum[temp$coding == 'cds'] = 1
temp$exonNum = 0; temp$exonNum[temp$exon == 'partial'] = 1; temp$exonNum[temp$exon == 'exon'] = 2
temp$transcribedNum = 0; temp$transcribedNum[temp$transcribed == 'partial'] = 1; temp$transcribedNum[temp$transcribed == 'transcribed'] = 2
temp$chrNum = gsub("chr", "", temp$chr); temp$chrNum = as.numeric(temp$chrNum)


# prepare pca data
pca = temp[temp$rate != 0, c('chrNum','start','end','period','consGC','allGC','totalLen','stdLenAll',
                             'op1','op2','segDupNum','codingNum','exonNum','transcribedNum','dist2Recom')]
pca = temp[temp$rate != 0, c('chrNum','period','consGC','allGC','totalLen','stdLenAll',
                             'op1','op2','segDupNum','codingNum','exonNum','transcribedNum','dist2Recom')]
out = prcomp(pca, center = T, scale. = T)
print( summary(out)$importance[,1:10] )

plotData = as.data.frame(out$x)[,1:10]
plotData = cbind(temp[temp$rate != 0,],plotData)
plotData$ratelog10 = log10(plotData$rate)

set.seed(123)
ggplot(data = plotData[sample(nrow(plotData), size = 1000, replace = FALSE), ]) +
  theme_classic() +
  geom_point(aes(x=PC1, y=PC2, color=ratelog10),cex=0.8) +
  xlab('PC1') + ylab('PC2') +
  scale_color_continuous(low = "red", high = "white", name='') +
  #scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')








####################
# show covered loci vs. some feature
####################
##### plot number of fully covered, rescued, not calculated loci, by period6
plotData = data.frame(table(rate$spanned,rate$period6))
colnames(plotData) = c('type','period6','freq'); plotData$total = 0
for (i in c('1','2','3','4','5','6','>6')) plotData$total[plotData$period6==i] = sum(plotData$freq[plotData$period6==i])
plotData$props = plotData$freq / plotData$total
plotData$type = factor(plotData$type, levels=c('noTree','rescued','full'))

colors = c('#d53e4f','#ffffbf','#3288bd')

png(file='TR-tree.png', width=1000, height=2000, res=300)
ggplot(data=plotData, aes(x=period6, y=props, fill=type)) + geom_bar(stat="identity") + theme_classic() +
  theme(text = element_text(size = 22), legend.position = 'top') +
  xlab('Motif period') + ylab("Proportion") +
  geom_hline(yintercept=0.94, linetype='dashed', color='orange', linewidth=0.5) +
  geom_text(aes(0, 0.89, label='0.94', hjust=-0.1), cex=4.5, color='orange') +
  geom_hline(yintercept=0.98, linetype='dashed', color='orange', linewidth=0.5) +
  geom_text(aes(0, 1.02, label='0.98', hjust=-0.1), cex=4.5, color='orange') +
  scale_fill_manual(name='', values=colors, labels = c('no nearby tree','nearby tree','fully spanned')) +
  guides(fill = guide_legend(nrow = 3))
dev.off()

##### plot number of fully covered, rescued, not calculated loci, by STR/VNTR
plotData = data.frame(table(rate$spanned,rate$type))
colnames(plotData) = c('type','catalog','freq'); plotData$total = 0
for (i in c('STR','VNTR')) plotData$total[plotData$catalog==i] = sum(plotData$freq[plotData$catalog==i])
plotData$props = plotData$freq / plotData$total
plotData$type = factor(plotData$type, levels=c('noTree','rescued','full'))

colors = c('#d53e4f','#ffffbf','#3288bd')

png(file='TR-tree.png', width=1000, height=2000, res=300)
ggplot(data=plotData, aes(x=catalog, y=props, fill=type)) + geom_bar(stat="identity") + theme_classic() +
  theme(text = element_text(size = 22), legend.position = 'top') +
  xlab('') + ylab("Proportion") +
  geom_hline(yintercept=0.94, linetype='dashed', color='orange', linewidth=0.5) +
  geom_text(aes(0, 0.89, label='0.94', hjust=-0.1), cex=4.5, color='orange') +
  geom_hline(yintercept=0.98, linetype='dashed', color='orange', linewidth=0.5) +
  geom_text(aes(0, 1.02, label='0.98', hjust=-0.1), cex=4.5, color='orange') +
  scale_fill_manual(name='', values=colors, labels = c('no nearby tree','nearby tree','fully spanned')) +
  guides(fill = guide_legend(nrow = 3))
dev.off()


# fully spanned loci
sum(rate$spanned == 'full')
sum(rate$spanned == 'full') / nrow(rate)
sum(rate$spanned == 'rescued')
sum(rate$spanned == 'rescued') / nrow(rate)
sum(rate$spanned != 'noTree')
sum(rate$spanned != 'noTree') / nrow(rate)

sum(rate$chr == 'chrY')



####################
# histogram: fully spanned loci by segDup for each period
####################
table(rate$segDup, rate$spanned)

period6 = factor(c('1','2','3','4','5','6','>6'), levels=c('1','2','3','4','5','6','>6'))
underlayer = data.frame(period6 = rep(period6,4),
                        segDup = c(rep('not in segDup',14), rep('in segDup',14)),
                        spanned = rep(c(rep('not fully spanned',7), rep('fully spanned',7)), 2),
                        y = c(rep(310000,14),rep(15000,14)))

fisher.test(matrix(c(964425,219197,38196,22038),2,2))

png(file='fully spanned loci vs. segDup.png', width=2000, height=1500, res=300)
ggplot(data = rate) + theme_classic() +
  geom_histogram(aes(x=period6, fill=period6), stat='count', position=position_dodge()) + # histogram
  geom_point(data = underlayer, aes(x=period6, y=y), colour = "white") +
  xlab('Motif period') + ylab('Number of loci') +
  facet_wrap(segDup~spanned, strip.position='bottom', scales='free_y', ncol=2) +
  theme(strip.text = element_text(size=15, face="bold")) +
  theme(legend.position="none") +
  theme(panel.grid.major  = element_line(colour = "white", linewidth = 0.2)) +
  theme(panel.grid.minor  = element_line(colour = "white", linewidth = 0.5)) +
  theme(text = element_text(size = 15))
dev.off()

sum(pull$spanned == 'fully spanned')



####################
# density plot:  rate by GC, overlay by disease loci
####################

png(file='mutRate vs. GC content.noHPRC2.png', width=3000, height=2000, res=300)
ggplot(data = rbind(rate[!is.na(rate$rateNormAfter) & rate$period != 1,],
                    rateBio[!is.na(rateBio$rateNormAfter),])) + theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=allGC), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('GC content') +
  scale_x_continuous(breaks = seq(-10, -2, by = 1)) +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data=rateBio[!is.na(rateBio$rateNormAfter),],
             aes(x=log10(rateNormAfter), y=allGC, pch=period6, color=coding), cex=2, stroke=1) +
  scale_shape_manual(name = 'period size', values = c(1,2,3,4,5)) +
  scale_color_manual(name = 'coding', values = c('pink','red'), labels = c('non-coding','coding')) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')
dev.off()

# check suspicious disease loci of ATTTT type c(2:4,19,30)
AT5Bio = rateBio[rateBio$allGC<0.05 & rateBio$period==5, c(2:4,18,15,21:27)]
AT5Bio = rateBio[rateBio$allGC<0.5 & rateBio$period==5, c(2:4,18,15,21:27)]

# check suspicious loci of ATTTT type
AT5genome = rate[rate$allGC<0.05 & rate$period==5
                 & !is.na(rate$rateNormAfter) & rate$rateNormAfter>0.000033
                 & rate$transcribed == 'transcribed'
                 & rate$size > 40 & rate$size < 120,
                 c(2:4,17,14,20:26)]

#AT5genome = AT5genome[sample(nrow(AT5genome)),]

sum(rateBio$totalEd != 'noTree' | is.na(rateBio$totalEd))

####################
# compare ATTTT sequence: genomic vs. large mutation rate
####################
# check 4-1 or 3-2 AT case (genome-wide)
AT5genome = rate[rate$allGC<0.05 & rate$period==5
                 & !is.na(rate$rateNormAfter) #& rate$rateNormAfter>0.000033
                 & rate$transcribed == 'transcribed'
                 & rate$size > 40 & rate$size < 120,
                 c(2:4,15,18,23:27)]
AT5genome$ATratio = 0
for (i in 1:nrow(AT5genome)) {
  cons = AT5genome$consensus[i]
  nAs = length(gregexpr("A", cons, TRUE)[[1]])
  nTs = length(gregexpr("T", cons, TRUE)[[1]])
  nCs = length(gregexpr("C", cons, TRUE)[[1]])
  nGs = length(gregexpr("G", cons, TRUE)[[1]])
  if (nAs == 4 | nTs == 4) {
    AT5genome$ATratio[i] = 'AT ratio = 1:4'
  } else if (nAs == 3 | nTs == 3) {
    AT5genome$ATratio[i] = 'AT ratio = 2:3'
  } else {
    AT5genome$ATratio[i] = 0
  }
}
sum(AT5genome$ATratio=='AT ratio = 1:4')
sum(AT5genome$ATratio=='AT ratio = 2:3')
sum(AT5genome$ATratio=='AT ratio = 1:4' & AT5genome$rateNormAfter > 0.000033)
sum(AT5genome$ATratio=='AT ratio = 2:3' & AT5genome$rateNormAfter > 0.000033)

fisher.test(matrix(c(133,194,173,6),2,2))

# plot distribution of mutation rate of AT rich loci
png(file='mutRate genomic AT5.png', width=2000, height=1000, res=300)
ggplot(data = AT5genome[AT5genome$ATratio!=0,]) + theme_classic() +
  geom_histogram(aes(x=log10(rateNormAfter)), colour = 4, fill = "white", bins = 50) +
  geom_vline(xintercept=log10(0.00005),lty='dashed',color='red') +
  facet_wrap(~ATratio, strip.position='bottom', scales='free_y', ncol=1) +
  theme(strip.text = element_text(size=15, face="bold")) +
  theme(legend.position="none") +
  xlab('Mutation rate (log10)') + ylab('Number of loci') +
  theme(text = element_text(size = 15))
dev.off()

####################
# plot block affected loci: ideogram
####################
library(karyoploteR)

temp = rate[rate$block=='True',]

png(file='block.all.ideogram.png', width=4000, height=4000, res=300)
kp <- plotKaryotype(genome='hg38')
points <- toGRanges(rate[rate$block=='True' & rate$chr!='chrX' & rate$chr!='chrY',2:4])
kpPlotDensity(kp, points, window.size=100000, col="#FFAACC")
dev.off()
png(file='block.VNTR.ideogram.png', width=4000, height=4000, res=300)
kp <- plotKaryotype(genome='hg38')
points <- toGRanges(rate[rate$block=='True' & rate$chr!='chrX' & rate$chr!='chrY' & rate$type=='VNTR',2:4])
kpPlotDensity(kp, points, window.size=100000, col="#FFAACC")
dev.off()





# density: heter comp vs. heter length, overlay by block loci
png(file='mutRate vs. GC content.noHPRC2.png', width=3000, height=2000, res=300)
ggplot(data = rate[!is.na(rate$rateNormAfter),]) + theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=nAffected), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('GC content') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data = rate[!is.na(rate$rateNormAfter) & rate$block == 'True',],
             aes(x=log10(rateNormAfter), y=allGC, pch=period6, color=coding), cex=2, stroke=1) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')
dev.off()

ggplot(data = rate[rate$nAffected > 0,]) + theme_classic() +
  geom_bin2d(aes(x=nAffected, y=allGC), bins=50) + 
  xlab('Mutation rate (log10)') + ylab('GC content') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  geom_point(data = rate[!is.na(rate$rateNormAfter) & rate$block == 'True',],
             aes(x=log10(rateNormAfter), y=allGC, pch=period6, color=coding), cex=2, stroke=1) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')














# hot spot of block positive loci?
ggplot(data = temp[temp$nAffected > 0,]) + theme_classic() +
  geom_histogram() + 
  xlab('Mutation rate (log10)') + ylab('GC content') +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')

# 2d density
nrow(rate[rate$prop > 0,])
png(file='mutRate vs. block severity.png', width=3000, height=2000, res=300)
ggplot(data = rate[rate$prop > 0,]) + theme_classic() +
  geom_bin2d(aes(x=log10(rateNormAfter), y=prop), bins=50) +
  xlab('Mutation rate (log10)') + ylab('Proportion of block-affected genome') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 25),
        legend.text = element_text(size = 20), legend.position = 'right')
dev.off()

# scatter
ggplot(data = rate[rate$prop > 0,]) + theme_classic() +
  geom_point(aes(x=log10(rateNormAfter), y=prop, color=segDup), cex=0.5) + 
  xlab('Mutation rate (log10)') + ylab('Proportion of block-affected genome') +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = 'right')

fullCols = c('chr','start','end','treeS','treeE','TR2Tree','nSNP','treeLen','winS','winE','period','type','annotated',
             'totalAllele','spanned','totalEd','totalDist','totalLen','rate','rateNormBefore','rateNormAfter','hotspot',
             'segDupIdentity','segDupCNV','period6','source','consensus','consGC','allGC','op1','op2','segDup','coding',
             'exon','transcribed','HAR','entropy','purityOri','purityEff','size','snpDensity','coverAll','stdLenAll',
             'numOutlierAll','homoLenAll','heteLenAll','homoComAll','heteComAll','lenAll','comAll','numPattern','aveBlock',
             'chrR','startR','endR','dist2Recom','telomere','aveBlockCat')
pickCols = c('chr','start','end','treeS','treeE','TR2Tree','nSNP','treeLen','winS','winE','period','type','annotated',
             'totalAllele','spanned','totalEd','totalDist','totalLen','rate','rateNormBefore','rateNormAfter','hotspot',
             'segDupIdentity','segDupCNV','period6','source','consensus','consGC','allGC','op1','op2','segDup','coding',
             'exon','transcribed','entropy','purityOri','purityEff','size','snpDensity','coverAll','stdLenAll',
             'numOutlierAll','homoLenAll','heteLenAll','homoComAll','heteComAll','lenAll','comAll','numPattern','aveBlock',
             'chrR','startR','endR','dist2Recom','telomere','aveBlockCat')
write.table(rate[,pickCols], 
            file = "rate.mergeAllFeatures.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

