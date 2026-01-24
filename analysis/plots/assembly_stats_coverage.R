options(scipen = 999)
library(ggplot2)
library(data.table)


#####
# HiFi assemblies published in year (donut plot)
# https://statdoe.com/pie-donut-chart-in-r/
#####
library(webr)
library(dplyr)
history = data.frame(batch=c('HGSVC-phase2','HPRC-phase1','CPC-phase1','APR-phase1','HGSVC-phase3','HPRC-phase2','JSA-phase1'),
                     n=c(35,47,57,50,65,232,19),
                     nRmDup=c(2,1,57,50,60,232,14),
                     nRmDupnoHPRC2=c(2,43,57,50,65,0,19),
                     year=c(2020,2023,2023,2024,2024,2024,2024))

# Number of public HiFi assemblies in major sequencing consortium
png(file='Num_ass_total.png', width=2000, height=1500, res=300)
PieDonut(history, aes(year, batch, count=n), ratioByGroup = FALSE,
         #title = '576 HiFi assemblies from 7 consortium',
         showRatioThreshold=0.001, labelpositionThreshold=0.08, labelposition=2,
         pieLabelSize=4, donutLabelSize=4, r0=0.3, r1=0.9, r2=1.2)
dev.off()

# Number of public HiFi assemblies in the current version of TRCompDB
png(file='Num_ass_TRCompDB.png', width=2000, height=1500, res=300)
PieDonut(history, aes(year, batch, count=nRmDup), ratioByGroup = FALSE,
         showRatioThreshold=0.001, labelpositionThreshold=0.1, labelposition=2,
         pieLabelSize=4, donutLabelSize=4, r0=0.3, r1=0.9, r2=1.2)
dev.off()

# Number of public HiFi assemblies in the noHPRC2 of TRCompDB
png(file='Num_ass_TRCompDB_noHPRC2.png', width=2000, height=1500, res=300)
PieDonut(history, aes(year, batch, count=nRmDupnoHPRC2), ratioByGroup = FALSE,
         showRatioThreshold=0.001, labelpositionThreshold=0.1, labelposition=2,
         pieLabelSize=4, donutLabelSize=4, r0=0.3, r1=0.9, r2=1.2)
dev.off()


#####
# bam coverage heatmap by batch and sample
#####
bamStats = read.csv('z_bamCoveredRegions_GRCh38_summary_all.tsv', sep='\t', header=T)
# remove the unpolished apr samples
bamStats = bamStats[2:nrow(bamStats),] # remove the total line
bamStats = bamStats[c(1:100,183:nrow(bamStats)),] # remove the unpolished APR-phase1 samples
bamStats$sample = paste0(bamStats$batch,'__',bamStats$sample)

# normalize each chromosome to 0-1
for (i in 28:50) {
  minimum = min(bamStats[2:nrow(bamStats),i-25])
  maximum = max(bamStats[2:nrow(bamStats),i-25])
  #bamStats[,i] = (bamStats[,i] - minimum) / (maximum - minimum)
  bamStats[2:nrow(bamStats),i] = bamStats[2:nrow(bamStats),i-25] / maximum
}

bamStatsProp = bamStats[,c(1:2,28:49)]
bamStatsProp$batch = as.factor(bamStatsProp$batch)
# change header to chr1-chr22
colnames(bamStatsProp) = colnames(bamStats)[1:24]

# reshape wide to long
longBam = reshape2::melt(bamStatsProp, id.vars = c('sample','batch'), variable.name = 'chrom', value.name = 'cover')

png(file='bamCover_by_sample_by_chr.png', width=6000, height=1000, res=300)
ggplot(data = longBam[longBam$batch != 'HPRC-phase1' & longBam$batch != 'HGSVC-phase2',]) + 
  (aes(x = chrom, y = sample, fill = cover)) +
  geom_tile() + xlab(label='chromosome') + ylab(label='') +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=0.9, hjust=1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.title = element_text(size=13), legend.text = element_text(size=10), legend.position = 'top',
        strip.text.x = element_text(size=13)) +
  scale_fill_continuous(low = "red", high = "white", name='covered proportion') +
  facet_wrap(~ batch, nrow=1, scales='free_y')
dev.off()



######################################################################
########## TR coverage heatmap by batch and sample                   #
######################################################################
trStats = read.csv('samples_1010_GRCh38_q-0.1.annoLen_trcoverage.tsv', sep='\t', header=T)
trStats$batch = 'APR-phase1'; trStats$batch[101:214] = 'CPC-phase1'
trStats$batch[215:284] = 'HGSVC-phase2'; trStats$batch[285:414] = 'HGSVC-phase3'
trStats$batch[415:508] = 'HPRC-phase1'; trStats$batch[509:972] = 'HPRC-phase2'; trStats$batch[973:1010] = 'JSA-phase1'
# reshape wide to long
longTR = reshape2::melt(trStats, id.vars = c('genome','batch'), variable.name = 'chrom', value.name = 'cover')
longTR = longTR[longTR$chrom!='chrX',]; longTR = longTR[longTR$chrom!='chrY',]

png(file='trCover_by_sample_by_chr.png', width=6000, height=1000, res=300)
ggplot(data = longTR[longTR$batch != 'HPRC-phase1' & longTR$batch != 'HGSVC-phase2',]) + 
  (aes(x = chrom, y = genome, fill = cover)) +
  geom_tile() + xlab(label='chromosome') + ylab(label='') +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=0.9, hjust=1),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.title = element_text(size=13), legend.text = element_text(size=10), legend.position = 'top',
        strip.text.x = element_text(size=13)) +
  scale_fill_continuous(low = "red", high = "white", name='covered proportion') +
  facet_wrap(~ batch, nrow=1, scales='free_y')
dev.off()

##### combine both bam and TR coverage plot #####
colnames(longBam)[1] = 'genome'
longBam$group = 'Assembly coverage'; longTR$group = 'TR coverage'
long = rbind(longBam,longTR)

png(file='cover_by_sample_by_chr.png', width=4500, height=2000, res=300)
ggplot(data = long[long$batch != 'HPRC-phase1' & long$batch != 'HGSVC-phase2',]) + 
  (aes(x = chrom, y = genome, fill = cover)) +
  geom_tile() + xlab(label='chromosome') + ylab(label='') +
  theme(axis.text.x = element_text(size=7, angle=45, vjust=0.9, hjust=1),
        axis.title.x = element_text(size=14, face='bold'),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.title = element_text(size=13), legend.text = element_text(size=10), legend.position = 'top',
        strip.text.x = element_text(size=15)) +
  scale_fill_continuous(low = "red", high = "white", name='proportion') +
  #facet_grid(batch ~ group, scales="free_y", switch = 'y')
  facet_wrap(group ~ batch, nrow=2, scales='free_y')
dev.off()


######################################################################
########## boxplot                                                   #
######################################################################
trStats$average = apply(trStats[,2:23], 1, mean)

##### boxplot for just TRs by chrom and batch #####
colors = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
batch_order = c('CPC-phase1','JSA-phase1','HGSVC-phase2','HGSVC-phase3','APR-phase1','HPRC-phase1','HPRC-phase2')
longTR$batch = factor(longTR$batch, levels=batch_order)
png(file='trCover_by_sample_by_chr.boxplot.png', width=6000, height=3000, res=300)
ggplot(data = longTR[longTR$batch != 'HPRC-phase1' & longTR$batch != 'HGSVC-phase2',]) + 
  theme_classic() +
  geom_boxplot(aes(x=batch, y=cover, fill=batch, color=batch)) +
  xlab(label='') + ylab(label='Coverage') +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=15), axis.ticks.y = element_blank(),
        #panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        #panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
        legend.title = element_text(size=15), legend.text = element_text(size=10), legend.position = 'top',
        strip.text.x = element_text(size=15)) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  facet_wrap(~ chrom, nrow=2, scales='free_y')
dev.off()

##### boxplot for just TRs by batch (average chromosome) #####
colors = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
batch_order = c('CPC-phase1','JSA-phase1','HGSVC-phase2','HGSVC-phase3','APR-phase1','HPRC-phase1','HPRC-phase2')
count = data.frame(batch = batch_order,
                   count = c(114,38,70,130,100,94,464),
                   pos=rep(1.008,7))
png(file='trCover_by_sample.boxplot.png', width=2500, height=1500, res=300)
ggplot(data = trStats[trStats$batch != 'HPRC-phase1' & trStats$batch != 'HGSVC-phase2',]) + 
  theme_classic() +
  geom_boxplot(aes(x=factor(batch, level=batch_order), y=average, fill=batch, color=batch)) +
  xlab(label='') + ylab(label='Average TR coverage on chr1-chr22') +
  ylim(0.85,1.015) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  geom_text(data=count[count$batch != 'HPRC-phase1' & count$batch != 'HGSVC-phase2',],
            aes(x=factor(batch, level=batch_order), y=pos, label=paste0('n = ',count)),
            cex=6) + # numbers for counts
  coord_flip() +
  theme(axis.text.x = element_text(size=20), axis.title.x = element_text(size=20),
        axis.text.y = element_text(size=20,angle=45,hjust=0.7), axis.title.y = element_text(size=20),
        legend.position = 'none',
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5))
dev.off()







