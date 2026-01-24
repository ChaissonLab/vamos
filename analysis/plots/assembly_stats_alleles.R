options(scipen = 999)
library(ggplot2)
#library(scales)
#library(ggrepel)
library(data.table)

####################
# config input allele data, merge TR catalog to it
####################
# data filtering and variable creation
#alleles = read.csv('samples_rmDup_rmKids_GRCh38_q-0.1.annoStr.tsv', sep='\t', header=T)
alleles = read.csv('samples_rmDup_rmKids_GRCh38_q-0.1.annoStr.3IQR.tsv', sep='\t', header=T)

# filter chrX,chrY, and homopolymers
alleles = alleles[alleles$chr != 'chrX',]; alleles = alleles[alleles$chr != 'chrY',]
#alleles = alleles[alleles$period != 1,]
# create the categorical motif period variable
alleles$period6 = alleles$period; alleles$period6[alleles$period6 > 6] = '>6'
alleles$period6 = factor(alleles$period6, levels=c('1','2','3','4','5','6','>6'))
alleles$key = paste0(alleles$chr,'_',alleles$start,'_',alleles$end)

# read TR catalog features
catalog = read.csv('../../vamosExpanded_v3.0_effMotifs-0.1.tsv', sep='\t', header=F)
colnames(catalog) = c('chr','start','end','motifs','source','type','period','consensus',
                      'consGC','allGC','op1','op2','segDup','coding','exon','transcribed')
catalog$key = paste0(catalog$chr,'_',catalog$start,'_',catalog$end)
catalog = catalog[,8:17]

# merge TR catalog features into mutation rate data
alleles = merge(alleles, catalog, by='key')


# separate datasets for each batch
overall = alleles[,c('chr','start','end','period','type','coverAll','homoLenAll','heteLenAll','homoComAll','heteComAll',
                     'lenAll','comAll','period6','consensus','consGC','allGC','op1','op2','segDup','coding','exon','transcribed')]
colnames(overall)[6:12] = c('cover','homoLen','heteLen','homoCom','heteCom','byLen','byCom')

overall$byLenRatio = overall$byLen / overall$cover
overall$byComRatio = overall$byCom / overall$cover
overall$variabilityLen = 'Constant'
#overall$variabilityLen[overall$byLenRatio > 0 & overall$byLenRatio < 0.05] = '< 0.05'
#overall$variabilityLen[overall$byLenRatio >= 0.05 & overall$byLenRatio < 0.1] = '[0.05,0.1)'
#overall$variabilityLen[overall$byLenRatio >= 0.1 & overall$byLenRatio < 0.25] = '[0.1,0.25)'
#overall$variabilityLen[overall$byLenRatio >= 0.25 & overall$byLenRatio < 0.5] = '[0.25-0.5)'
#overall$variabilityLen[overall$byLenRatio >= 0.5] = '\u2265 0.5'
#overall$variabilityLen[overall$byLen == 1] = 'Constant'
overall$variabilityLen[overall$byLenRatio > 0 & overall$byLenRatio < 0.02] = '< 0.02'
overall$variabilityLen[overall$byLenRatio >= 0.02 & overall$byLenRatio < 0.05] = '[0.02,0.05)'
overall$variabilityLen[overall$byLenRatio >= 0.05 & overall$byLenRatio < 0.1] = '[0.05,0.1)'
overall$variabilityLen[overall$byLenRatio >= 0.1 & overall$byLenRatio < 0.2] = '[0.1-0.2)'
overall$variabilityLen[overall$byLenRatio >= 0.2] = '\u2265 0.2'
overall$variabilityLen[overall$byLen == 1] = 'Constant'

overall$variabilityCom = 'Constant'
#overall$variabilityCom[overall$byComRatio > 0 & overall$byComRatio < 0.05] = '< 0.05'
#overall$variabilityCom[overall$byComRatio >= 0.05 & overall$byComRatio < 0.1] = '[0.05,0.1)'
#overall$variabilityCom[overall$byComRatio >= 0.1 & overall$byComRatio < 0.25] = '[0.1,0.25)'
#overall$variabilityCom[overall$byComRatio >= 0.25 & overall$byComRatio < 0.5] = '[0.25-0.5)'
#overall$variabilityCom[overall$byComRatio >= 0.5] = '\u2265 0.5'
#overall$variabilityCom[overall$byCom == 1] = 'Constant'
overall$variabilityCom[overall$byComRatio > 0 & overall$byComRatio < 0.02] = '< 0.02'
overall$variabilityCom[overall$byComRatio >= 0.02 & overall$byComRatio < 0.05] = '[0.02,0.05)'
overall$variabilityCom[overall$byComRatio >= 0.05 & overall$byComRatio < 0.1] = '[0.05,0.1)'
overall$variabilityCom[overall$byComRatio >= 0.1 & overall$byComRatio < 0.2] = '[0.1-0.2)'
overall$variabilityCom[overall$byComRatio >= 0.2] = '\u2265 0.2'
overall$variabilityCom[overall$byCom == 1] = 'Constant'

# overall statistics on loci variability
sum(table(overall$variabilityLen)) == nrow(overall)
table(overall$variabilityLen)
prop.table(table(overall$variabilityLen))
prop.table(table(overall$variabilityCom))

# average number of alleles plot
allelesAPR = alleles[,c('chr','start','end','period','type','cover_APR.phase1','len_APR.phase1','com_APR.phase1','period6')]
colnames(allelesAPR)[6:8] = c('cover','byLen','byCom'); allelesAPR$batch = 'APR-phase1'
allelesCPC = alleles[,c('chr','start','end','period','type','cover_CPC.phase1','len_CPC.phase1','com_CPC.phase1','period6')]
colnames(allelesCPC)[6:8] = c('cover','byLen','byCom'); allelesCPC$batch = 'CPC-phase1'
allelesHGSVC = alleles[,c('chr','start','end','period','type','cover_HGSVC.phase3','len_HGSVC.phase3','com_HGSVC.phase3','period6')]
colnames(allelesHGSVC)[6:8] = c('cover','byLen','byCom'); allelesHGSVC$batch = 'HGSVC-phase3'
allelesHPRC = alleles[,c('chr','start','end','period','type','cover_HPRC.phase2','len_HPRC.phase2','com_HPRC.phase2','period6')]
colnames(allelesHPRC)[6:8] = c('cover','byLen','byCom'); allelesHPRC$batch = 'HPRC-phase2'

# prepare plot data
plotData = rbind(allelesAPR,allelesCPC,allelesHGSVC,allelesHPRC)
plotData$byLenRatio = plotData$byLen / plotData$cover
plotData$byComRatio = plotData$byCom / plotData$cover

plotData$variabilityLen = 'Constant'
#plotData$variabilityLen[plotData$byLenRatio > 0 & plotData$byLenRatio < 0.05] = '< 0.05'
#plotData$variabilityLen[plotData$byLenRatio >= 0.05 & plotData$byLenRatio < 0.1] = '[0.05,0.1)'
#plotData$variabilityLen[plotData$byLenRatio >= 0.1 & plotData$byLenRatio < 0.25] = '[0.1,0.25)'
#plotData$variabilityLen[plotData$byLenRatio >= 0.25 & plotData$byLenRatio < 0.5] = '[0.25-0.5)'
#plotData$variabilityLen[plotData$byLenRatio >= 0.5] = '\u2265 0.5'
#plotData$variabilityLen[plotData$byLen == 1] = 'Constant'
plotData$variabilityLen[plotData$byLenRatio > 0 & plotData$byLenRatio < 0.02] = '< 0.02'
plotData$variabilityLen[plotData$byLenRatio >= 0.02 & plotData$byLenRatio < 0.05] = '[0.02,0.05)'
plotData$variabilityLen[plotData$byLenRatio >= 0.05 & plotData$byLenRatio < 0.1] = '[0.05,0.1)'
plotData$variabilityLen[plotData$byLenRatio >= 0.1 & plotData$byLenRatio < 0.2] = '[0.1-0.2)'
plotData$variabilityLen[plotData$byLenRatio >= 0.2] = '\u2265 0.2'
plotData$variabilityLen[plotData$byLen == 1] = 'Constant'

plotData$variabilityCom = 'Constant'
#plotData$variabilityCom[plotData$byComRatio > 0 & plotData$byComRatio < 0.05] = '< 0.05'
#plotData$variabilityCom[plotData$byComRatio >= 0.05 & plotData$byComRatio < 0.1] = '[0.05,0.1)'
#plotData$variabilityCom[plotData$byComRatio >= 0.1 & plotData$byComRatio < 0.25] = '[0.1,0.25)'
#plotData$variabilityCom[plotData$byComRatio >= 0.25 & plotData$byComRatio < 0.5] = '[0.25-0.5)'
#plotData$variabilityCom[plotData$byComRatio >= 0.5] = '\u2265 0.5'
#plotData$variabilityCom[plotData$byCom == 1] = 'Constant'
plotData$variabilityCom[plotData$byComRatio > 0 & plotData$byComRatio < 0.02] = '< 0.02'
plotData$variabilityCom[plotData$byComRatio >= 0.02 & plotData$byComRatio < 0.05] = '[0.02,0.05)'
plotData$variabilityCom[plotData$byComRatio >= 0.05 & plotData$byComRatio < 0.1] = '[0.05,0.1)'
plotData$variabilityCom[plotData$byComRatio >= 0.1 & plotData$byComRatio < 0.2] = '[0.1-0.2)'
plotData$variabilityCom[plotData$byComRatio >= 0.2] = '\u2265 0.2'
plotData$variabilityCom[plotData$byLen == 1] = 'Constant'

#plotData$variabilityLen = factor(plotData$variabilityLen,
#                                 levels=c('\u2265 0.5','[0.25-0.5)','[0.1,0.25)','[0.05,0.1)','< 0.05','Constant'))
#plotData$variabilityCom = factor(plotData$variabilityCom,
#                                 levels=c('\u2265 0.5','[0.25-0.5)','[0.1,0.25)','[0.05,0.1)','< 0.05','Constant'))
plotData$variabilityLen = factor(plotData$variabilityLen,
                                 levels=c('\u2265 0.2','[0.1-0.2)','[0.05,0.1)','[0.02,0.05)','< 0.02','Constant'))
plotData$variabilityCom = factor(plotData$variabilityCom,
                                 levels=c('\u2265 0.2','[0.1-0.2)','[0.05,0.1)','[0.02,0.05)','< 0.02','Constant'))

# test against coding regions (all loci)
t.test(overall$heteLen[overall$coding == '.'], overall$heteLen[overall$coding != '.'])
0.10506299 / 0.01704702
t.test(overall$heteCom[overall$coding == '.'], overall$heteCom[overall$coding != '.'])
0.11963445 / 0.02598308

# test against coding regions (only variable loci)
t.test(overall$heteLen[overall$homoLen != 1 & overall$coding == '.'], overall$heteLen[overall$homoLen != 1 & overall$coding != '.'])
0.2637881 / 0.1481231
t.test(overall$heteCom[overall$homoCom != 1 & overall$coding == '.'], overall$heteCom[overall$homoCom != 1 & overall$coding != '.'])
0.2718048 / 0.1595343

# ratio of heterozygosity by composition / heterozygosity by length
mean(overall[overall$heteLen!=0,]$heteCom / overall[overall$heteLen!=0,]$heteLen)

# filter
plotData = plotData[plotData$cover >= 20,]
for (p in c('1','2','3','4','5','6','>6')) {
  print( nrow(plotData[plotData$period6==p,]) )
}
plotData = plotData[plotData$chr != 'chrX',]
plotData = plotData[plotData$chr != 'chrY',]
#plotData = plotData[plotData$period != 1,]

# plot
colors = c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd')
png(file='variability_GRCh38_byLen.png', width=2500, height=2000, res=300)
ggplot(data=plotData) + theme_classic() +
  geom_histogram(aes(x=period6, fill=variabilityLen), stat='count') + facet_wrap(vars(batch)) +
  ggtitle('Average number of TR alleles per annotated haplotype') +
  xlab('Period size') + ylab('Frequency') +
  scale_fill_manual(name='', values=colors) +
  theme(text = element_text(size = 30), plot.title = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()

png(file='variability_GRCh38_byCom.png', width=2500, height=2000, res=300)
ggplot(data=plotData) + theme_classic() +
  geom_histogram(aes(x=period6, fill=variabilityCom), stat='count') + facet_wrap(vars(batch)) +
  ggtitle('Average number of TR alleles per annotated haplotype') +
  xlab('Period size') + ylab('Frequency') +
  scale_fill_manual(name='', values=colors) +
  theme(text = element_text(size = 30), plot.title = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5)) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
dev.off()


####################
# plot variability, by proportion for each motif size
####################
# calculate the proportion for each motif period group
plotData = data.frame(matrix(ncol = 4, nrow = 0))
colnames(plotData) = c('period6', 'variability', 'measure', 'props')
for (p in c('1','2','3','4','5','6','>6')) {
  for (v in c('\u2265 0.2','[0.1-0.2)','[0.05,0.1)','[0.02,0.05)','< 0.02','Constant')) {
    temp = data.frame(period6=c(p,p),
                      variability=c(v,v),
                      measure=c('len','com'),
                      props=c(-1,-1))
    temp$props[1] = nrow(overall[overall$period6==p & overall$variabilityLen==v,]) / nrow(overall[overall$period6==p,])
    temp$props[2] = nrow(overall[overall$period6==p & overall$variabilityCom==v,]) / nrow(overall[overall$period6==p,])
    plotData = rbind(plotData, temp)
  }
}
plotData$period6 = factor(plotData$period6, levels=c('1','2','3','4','5','6','>6'))
plotData$variability = factor(plotData$variability, levels=c('\u2265 0.2','[0.1-0.2)','[0.05,0.1)','[0.02,0.05)','< 0.02','Constant'))
plotData$measure = factor(plotData$measure, levels=c('len','com'))

colors = c('#d53e4f','#fc8d59','#fee08b','#ffffbf','#e6f598','#99d594','#3288bd')
colors = c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd')
#colors = c('#b2182b','#ef8a62','#fddbc7','#ffffff','#e0e0e0','#999999','#4d4d4d')

png(file='variability_GRCh38_props.png', width=2500, height=2000, res=300)
ggplot(data=plotData) + theme_classic() + 
  geom_bar(aes(x=measure, y=props, fill=variability), stat="identity") + 
  ggtitle('Average number of TR alleles per annotated haplotype') +
  xlab('Motif period') + ylab("Proportion") +
  scale_fill_manual('', values=colors) +
  scale_alpha_manual('', labels=c('by length','by composition'), values=c(0.6,1), ) +
  geom_text(aes(x=measure, y=1, label=measure), nudge_y=0.05, cex=6) +
  facet_wrap(~period6, strip.position='bottom', scales='free_x', nrow=1) + # set multiple group labels on x axis
  theme(text = element_text(size = 30), plot.title = element_text(size = 20, hjust = 0.5),
        legend.text = element_text(size = 20), legend.position = 'top',
        #panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        #panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
        axis.text.x = element_blank(), axis.ticks.x=element_blank()) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE), alpha='none')
dev.off()






####################
# 2d density: heter comp vs. heter length, overlay with disease loci
####################
# summary statistics increase of heterozygosity by composition from length
temp = overall[overall$heteLen > 0,]
temp$heteIncrease = temp$heteCom - temp$heteLen
temp$heteRatio = temp$heteCom / temp$heteLen
mean(temp$heteIncrease) # all loci
mean(temp$heteIncrease[temp$heteLen == 0]) # loci that are constant by length
mean(temp$heteIncrease[temp$heteLen > 0]) # loci that are variable by length
mean(temp$heteIncrease[temp$heteLen == 0 & temp$heteCom > 0]) # loci that are constant by length but variable by composition
mean(temp$heteIncrease[temp$heteLen > 0 & temp$heteCom > 0]) # loci that are variable by length and composition

mean(temp$heteRatio) # for all variable loci

# config disease data
allelesbio = read.csv('samples_rmDup_rmKids_GRCh38_q-0.1.annoStr.special.tsv', sep='\t', header=T)

allelesbio = allelesbio[allelesbio$chr != 'chrX',]
allelesbio = allelesbio[allelesbio$chr != 'chrY',]
#alleles = alleles[alleles$period != 1,]
allelesbio$period6 = allelesbio$period; allelesbio$period6[allelesbio$period6 > 6] = '>6'
allelesbio$period6 = factor(allelesbio$period6, levels=c('1','2','3','4','5','6','>6'))
allelesbio$key = paste0(allelesbio$chr,'_',allelesbio$start,'_',allelesbio$end)

catalogBio = read.csv('../../vamosBioTR_v3.0_effMotifs-0.1.tsv', sep='\t', header=F)
colnames(catalogBio) = c('chr','start','end','motifs','version','type','period','consensus',
                      'consGC','allGC','op1','op2','segDup','coding','exon','transcribed')
catalogBio$key = paste0(catalogBio$chr,'_',catalogBio$start,'_',catalogBio$end)
catalogBio = catalogBio[,8:17]

# merge TR catalog features into mutation rate data
allelesbio = merge(allelesbio, catalogBio, by='key')

overallBio = allelesbio[,c('chr','start','end','period','type','coverAll','homoLenAll','heteLenAll','homoComAll','heteComAll',
                          'lenAll','comAll','period6','consensus','consGC','allGC','op1','op2','segDup','coding','exon','transcribed')]
colnames(overallBio)[6:12] = c('cover','homoLen','heteLen','homoCom','heteCom','len','com')

# pull gene information of the biological loci
disease = read.csv('../../disease.tsv', sep='\t', header=T)

pull = overallBio[overallBio$start == disease$refined_start[1],]; pull$gene = disease$gene[1]
for (i in 2:nrow(disease)) {
  temp = overallBio[overallBio$start == disease$refined_start[i],]
  if (dim(temp)[1] != 0) {
    temp$gene = disease$gene[i]
    pull = rbind(pull,temp)
  }
}
pull$coding[pull$coding=='partial'] = 'cds'

# visual guide lines
y = seq(0,0.5,0.02)
x1 = (y+1) / 2 - 0.02
x2 = (y+1) / 2 + 0.02
x3 = (y+1) / 2
lineData = data.frame(x1=x1, y1=y, x2=x2, y2=y, x3=x3, y3=y)

# density: heter comp vs. heter length
png(file='heterozygosity com vs. heterozygosity len.png', width=3000, height=2500, res=300)
ggplot(data = overall) + theme_classic() +
  geom_bin2d(aes(x=heteCom, y=heteLen),bins=100) + 
  xlab('Heterozygosity by composition (HBC)') + ylab('Heterozygosity by length (HBL)') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  #scale_fill_continuous(type = "viridis", trans = "log10", labels = scales::comma) +
  #scale_fill_continuous(type = "viridis", trans = "log10", labels = scales::number) +
  geom_point(data=lineData, aes(x=x3,y=y3), cex=0.3, color='white') +
  geom_point(data=pull,
             aes(x=heteCom, y=heteLen, pch=period6, color=coding), cex=4, stroke=1) +
  scale_shape_manual(name = 'period size', values = c(1,2,3,4,5)) +
  scale_color_manual(name = 'coding', values = c('pink','red'), labels = c('non-coding','coding')) +
  guides(shape = guide_legend(ncol = 2)) + # This line changes the legend to two columns
  theme(text = element_text(size = 30), legend.title = element_text(size = 30),
        legend.text = element_text(size = 20), legend.position = c(0.2,0.6),
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()

# RAI1 gene in the middle
# HOXD13 gene at the bottom

# check some numbers and do a simple fisher's exact test
sum(pull$heteCom >= 0.3) / nrow(pull)
sum(overall$heteCom < 0.3) / nrow(overall)

temp = matrix(c(sum(pull$heteCom >= 0.3),sum(pull$heteCom < 0.3),sum(overall$heteCom >= 0.3),sum(overall$heteCom < 0.3)),2,2)
fisher.test(temp)


# if someone only measures length, he's gonna miss a lot of information
# (variation due to com is more dominant than length)
# but does this mean that motif mutations are more common than insertion deletions?

# if the IBD violates correlation, it is possible for reccurent mutation?
# should total length be accounted for here?








# density: # allele comp vs. allele length
ggplot(data = overall) + theme_classic() +
  geom_bin2d(aes(x=com, y=len),bins=100) + 
  xlab('Average number of alleles (by composition)') + ylab('Average number of alleles (by length)') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  #xlim(0,800) + ylim(0,800) +
  theme(legend.position = 'right')
ggplot(data = cpc.phase1) + theme_classic() +
  geom_bin2d(aes(x=com, y=len),bins=100) + 
  xlab('Average number of alleles (by composition)') + ylab('Average number of alleles (by length)') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  xlim(0,150) + ylim(0,150) +
  theme(legend.position = 'right')
ggplot(data = hgsvc.phase3) + theme_classic() +
  geom_bin2d(aes(x=com, y=len),bins=100) + 
  xlab('Average number of alleles (by composition)') + ylab('Average number of alleles (by length)') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  xlim(0,150) + ylim(0,150) +
  theme(legend.position = 'right')

# average number of alleles vs. number of annotated haplotype
ggplot(data = overall) + theme_classic() +
  geom_bin2d(aes(x=coverAll, y=lenAll),bins=200) + 
  xlab('Number of annotated haplotypes') + ylab('Average number of alleles (by length)') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(legend.position = 'right')
ggplot(data = overall) + theme_classic() +
  geom_bin2d(aes(x=coverAll, y=comAll),bins=200) + 
  xlab('Number of annotated haplotypes') + ylab('Average number of alleles (by composition)') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(legend.position = 'right')


overall[overall$start==57024495,] # WDR7
overall[overall$start==3074877,] # HTT
overall[overall$start==155188482,] # MUC1
overall[overall$start==88855423,] # ACAN



temp = overall[overall$homoLen / overall$homoCom > 1.9 & overall$homoLen / overall$homoCom < 2.1,]
temp = overall[overall$homoComAll == 0.5,]



####################
# ANOVA
####################
library(ggplot2)

# anova for all samples
scope = 'all'
anovaLen = read.csv(paste0(scope,'.anova.annoLen.tsv'), sep='\t', header=T)
anovaCom = read.csv(paste0(scope,'.anova.ed2Major.tsv'), sep='\t', header=T)
anovaLen$key = paste0(anovaLen$chr,'_',anovaLen$start,'_',anovaLen$end)
anovaCom$key = paste0(anovaCom$chr,'_',anovaCom$start,'_',anovaCom$end)

catalog = read.csv('../../vamosExpanded_v3.0_effMotifs-0.1.tsv', sep='\t', header=F)
colnames(catalog) = c('chr','start','end','motifs','source','type','period','consensus',
                      'consGC','allGC','op1','op2','segDup','coding','exon','transcribed')
catalog$key = paste0(catalog$chr,'_',catalog$start,'_',catalog$end)
catalog = catalog[,6:17]

anovaLen = merge(anovaLen, catalog, by='key')
anovaCom = merge(anovaCom, catalog, by='key')
anovaLen$log10p = log10(anovaLen$p)
anovaCom$log10p = log10(anovaCom$p)
anovaLen$sigfinicant = (anovaLen$log10p < -10) * 1
anovaCom$sigfinicant = (anovaCom$log10p < -10) * 1
anovaLen$codingYN = anovaLen$coding != '.'
anovaCom$codingYN = anovaCom$coding != '.'
anovaLen$period6 = anovaLen$period; anovaLen$period6[anovaLen$period6 > 6] = '>6'
anovaLen$period6 = factor(anovaLen$period6, levels=c('1','2','3','4','5','6','>6'))
anovaCom$period6 = anovaCom$period; anovaCom$period6[anovaCom$period6 > 6] = '>6'
anovaCom$period6 = factor(anovaCom$period6, levels=c('1','2','3','4','5','6','>6'))

# plot by STR/VNTR and len/com, excluding homopolymers
plotData = merge(anovaLen,anovaCom[,c('log10p','key')], by='key')[,c('period','key','type','log10p.x')]
plotData$test = 'by length'
temp = merge(anovaLen,anovaCom[,c('log10p','key')], by='key')[,c('period','key','type','log10p.y')]
temp$test = 'by composition'
colnames(plotData) = c('period','key','type','log10p','test')
colnames(temp) = c('period','key','type','log10p','test')
plotData = rbind(plotData, temp)
plotData$tag = paste0(plotData$type, ' ', plotData$test)
plotData$log10p[plotData$log10p<=-30] = -30
#temp = plotData[1:4,]; temp$tag = unique(plotData$tag); temp$log10p = -Inf
#plotData = rbind(plotData, temp)


# ratio: length-significant / length-variable
table(anovaLen$codingYN, anovaLen$sigfinicant); fisher.test(table(anovaLen$codingYN, anovaLen$sigfinicant))
table(anovaCom$codingYN, anovaCom$sigfinicant); fisher.test(table(anovaCom$codingYN, anovaCom$sigfinicant))
nrow(anovaLen[anovaLen$log10p < -10,]) == sum(table(anovaLen[anovaLen$log10p < -10,]$coding))
nrow(anovaLen[anovaLen$log10p < -10,]) / nrow(anovaLen)
# ratio: composition-significant / composition-variable
nrow(anovaCom[anovaCom$log10p < -10,]) == sum(table(anovaCom[anovaCom$log10p < -10,]$period6))
nrow(anovaCom[anovaCom$log10p < -10,]) / nrow(anovaCom)
# how much more by com vs. len
sum( table(anovaCom[anovaCom$log10p < -10,]$period6) ) / sum( table(anovaLen[anovaLen$log10p < -10,]$period6) )
# exclude homopolymers
sum( table(anovaLen[anovaLen$log10p < -10 & anovaLen$period > 1,]$period6) ) / nrow(anovaLen[anovaLen$period != 1,])

# n: variable by composition but length
nrow(anovaCom) - nrow(anovaLen)
# n: composition-significant, in loci that are variable by composition but length
nrow(anovaCom[anovaCom$log10p < -10,]) - 
  nrow(plotData[plotData$test == 'by composition' & plotData$log10p < -10,])

# number of significant loci on the same set of loci (the length-variable set)
# n: length-variable
nrow(plotData[plotData$test == 'by length',])
# n: length-significant
nrow(plotData[plotData$test == 'by length' & plotData$log10p < -10,])
# ratio: length-significant / length-variable
nrow(plotData[plotData$test == 'by length' & plotData$log10p < -10,]) / nrow(plotData[plotData$test == 'by length',])
# n: composition-significant
nrow(plotData[plotData$test == 'by composition' & plotData$log10p < -10,])
# ratio: composition-significant / length-significant
nrow(plotData[plotData$test == 'by composition' & plotData$log10p < -10,]) / nrow(plotData[plotData$test == 'by length' & plotData$log10p < -10,])

png(file=paste0(scope,'.annova.LenvsCom.cum.png'), width=1800, height=2700, res=300)
colors = c('#d53e4f','#fc8d59','#1d91c0','#081d58')
ggplot(data = plotData) +
  theme_classic() +
  stat_ecdf(aes(x=log10p, color=tag), linewidth=1, geom="line", position="identity", pad=F) +
  xlab('ANOVA p-value (log10)') + ylab('Cumulative density') +
  scale_color_manual(name='', values=colors) +
  theme(text = element_text(size = 30),
        legend.text = element_text(size = 23), legend.position = c(0.45,0.87),
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5)) +
  guides(color=guide_legend(nrow=4,byrow=TRUE), alpha='none') # put legend to two rows
dev.off()

# anovaLen$lengthCat = anovaLen$grandMean
# anovaLen$lengthCat[anovaLen$grandMean > 0 & anovaLen$grandMean <= 10] = '0-10'
# anovaLen$lengthCat[anovaLen$grandMean > 10 & anovaLen$grandMean <= 20] = '10-20'
# anovaLen$lengthCat[anovaLen$grandMean > 20 & anovaLen$grandMean <= 30] = '20-30'
# anovaLen$lengthCat[anovaLen$grandMean > 30 & anovaLen$grandMean <= 40] = '30-40'
# anovaLen$lengthCat[anovaLen$grandMean > 40 & anovaLen$grandMean <= 50] = '40-50'
# anovaLen$lengthCat[anovaLen$grandMean > 50] = '>50'
# anovaLen$lengthCat = factor(anovaLen$lengthCat, levels=c('0-10','10-20','20-30','30-40','40-50','>50'))
# anovaCom$lengthCat = anovaCom$grandMean
# anovaCom$lengthCat[anovaCom$grandMean > 0 & anovaCom$grandMean <= 10] = '0-10'
# anovaCom$lengthCat[anovaCom$grandMean > 10 & anovaCom$grandMean <= 20] = '10-20'
# anovaCom$lengthCat[anovaCom$grandMean > 20 & anovaCom$grandMean <= 30] = '20-30'
# anovaCom$lengthCat[anovaCom$grandMean > 30 & anovaCom$grandMean <= 40] = '30-40'
# anovaCom$lengthCat[anovaCom$grandMean > 40 & anovaCom$grandMean <= 50] = '40-50'
# anovaCom$lengthCat[anovaCom$grandMean > 50] = '>50'
# anovaCom$lengthCat = factor(anovaCom$lengthCat, levels=c('0-10','10-20','20-30','30-40','40-50','>50'))
# 
# anovaLen$lengthCat = anovaLen$grandMean
# anovaLen$lengthCat[anovaLen$grandMean > 0 & anovaLen$grandMean <= 10] = '\u2264 10'
# anovaLen$lengthCat[anovaLen$grandMean > 10] = '>10'
# anovaLen$lengthCat = factor(anovaLen$lengthCat, levels=c('\u2264 10','>10'))
# anovaCom$lengthCat = anovaCom$grandMean
# anovaCom$lengthCat[anovaCom$grandMean > 0 & anovaCom$grandMean <= 10] = '\u2264 10'
# anovaCom$lengthCat[anovaCom$grandMean > 10] = '>10'
# anovaCom$lengthCat = factor(anovaCom$lengthCat, levels=c('\u2264 10','>10'))

# png(file=paste0(scopre,'.anova.annoLen.png'), width=2500, height=2000, res=300)
# colors = c('#d53e4f','#3288bd')
# ggplot(data = anovaLen[anovaLen$log10p > -20 & anovaLen$period != 1,]) +
#   theme_classic() +
#   #stat_density(aes(x = log10p, color = lengthCat), geom = "line", position = "identity") +
#   stat_ecdf(aes(x = log10p, color = lengthCat), geom = "line", position = "identity") +
#   xlab('ANOVA p-value (log10)') +
#   #ylim(c(0,0.25)) +
#   scale_color_manual(name='Length (by motifs)', values=colors) +
#   theme(text = element_text(size = 20),
#         legend.text = element_text(size = 20), legend.position = 'top',
#         panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
#         panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
# dev.off()
# 
# png(file=paste0(scopre,'.anova.annoCom.png'), width=2500, height=2000, res=300)
# colors = c('#d53e4f','#fc8d59','#fee08b','#e6f598','#99d594','#3288bd')
# ggplot(data = anovaCom[anovaCom$log10p > -20 & anovaCom$period != 1,]) +
#   theme_classic() +
#   #stat_density(aes(x = log10p, color = lengthCat), geom = "line", position = "identity") +
#   stat_ecdf(aes(x = log10p, color = lengthCat), geom = "line", position = "identity") +
#   xlab('ANOVA p-value (log10)') +
#   #ylim(c(0,0.25)) +
#   scale_color_manual(name='Length (by motifs)', values=colors) +
#   theme(text = element_text(size = 20),
#         legend.text = element_text(size = 20), legend.position = 'top',
#         panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
#         panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
# dev.off()


# # check how many loci show significance on Com over Len
# plotData = merge(anovaLen,anovaCom[,c('log10p','key')], by='key')
# 
# # density, may not use, need to explain the bias issue
# ggplot(data = plotData) + theme_classic() +
#   geom_bin2d(aes(x=log10p.x, y=log10p.y),bins=100,) + 
#   xlab('Length p-value (log10)') + ylab('Composition p-value (log10)') +
#   xlim(c(-50,0)) + ylim(c(-50,0)) +
#   geom_vline(xintercept=-15, color="red", linetype="dashed", size = 0.5) +
#   geom_hline(yintercept=-15, color="red", linetype="dashed", size = 0.5) +
#   scale_fill_continuous(type = "viridis", trans = "log10") +
#   theme(text = element_text(size = 30),
#         legend.text = element_text(size = 20), legend.position = 'right',
#         panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
#         panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))





# check for the 50% enrichment
temp = overall[overall$heteLen > 0.48 & overall$heteLen < 0.52,]
write.table(temp, file = "selected0.5.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


hist(overall[overall$freq1+overall$freq2 > 0.95 & overall$freq1 < 0.95,]$heteLen, breaks = seq(0,1,0.01))

hist(overall[overall$freq1+overall$freq2 < 0.95 & overall$freq1 < 0.95 & overall$freq1 >0.5,]$heteLen, breaks = seq(0,1,0.01))


# random simulate freqs to see if homozygosity enrich at 50%
simu = data.frame(x1=numeric(0),x2=numeric(0),x3=numeric(0),heter=numeric(0))
for (i in seq(0,1,0.005)) {
  for (j in seq(0,1,0.005)) {
    for (k in seq(0,1,0.005)) {
      if (i+j+k == 1 & i >= j & j >= k) {
        temp = data.frame(x1=i,x2=j,x3=k,heter=1-i*i-j*j-k*k)
        simu = rbind(simu,temp)
      }
    }
  }
}
hist(simu[simu$x3 < 0.05,]$heter,breaks = seq(0,1,0.01))

ggplot(data = simu[simu$heter > 0.45 & simu$heter < 0.55,]) +
  theme_classic() +
  geom_bin2d(aes(x=x1, y=x2),bins=50,) + 
  xlim(0,1) + ylim(0,1) +
  xlab('First major allele freq') + ylab('Second major allele freq') +
  scale_fill_continuous(type = "viridis") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))


# plot the first 2 major frequencies
freq = read.csv('alleleLenFreqs.tsv', sep='\t', header=T)
freq$key = paste0(freq$chr,'_',freq$start,'_',freq$end)
freq = freq[,c('key','freq1','freq2')]

# merge TR catalog features into mutation rate data
overall$key = paste0(overall$chr,'_',overall$start,'_',overall$end)
overall = merge(overall, freq, by='key')

hist(overall[overall$byLen > 1,]$freq1)
hist(overall[overall$byLen > 1,]$freq2)

ggplot(data = overall[overall$byLen > 1 & overall$heteLen > 0.45 & overall$heteLen < 0.55,]) +
  theme_classic() +
  geom_bin2d(aes(x=freq1, y=freq2),bins=50,) + 
  xlim(0,1) + ylim(0,0.5) +
  xlab('First major allele freq') + ylab('Second major allele freq') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))

ggplot(data = overall[overall$byLen > 1,]) +
  theme_classic() +
  geom_bin2d(aes(x=freq1, y=freq2),bins=50,) + 
  xlim(0,1) +
  xlab('First major allele freq') + ylab('Second major allele freq') +
  scale_fill_continuous(type = "viridis", trans = "log10") +
  theme(text = element_text(size = 20),
        legend.text = element_text(size = 20), legend.position = 'top',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))

