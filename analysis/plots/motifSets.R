library(ggplot2)
library(scales)
library(ggrepel)

##################################################
# read and process data
##################################################
# read and process GRCh38 data
data = read.csv('vamosExpanded_v3.0_remainingInfoGRCh38.tsv', sep='\t', header=T); data$ref = 'GRCh38'
# read and process T2T data
temp = read.csv('vamosExpanded_v3.0_remainingInfoT2T-CHM13v2.0.tsv', sep='\t', header=T); temp$ref = 'T2T-CHM13'
data = rbind(data,temp); temp = 0

# create extra columns
data$info = round(data$remain / data$total, 2)
data$infoU = round(data$remainU / data$totalU, 2)
data$size = data$end - data$start + 1
data$tag = as.character(expression('\u2264 10 unique motifs')); data$tag[data$totalU>10] = '> 10 unique motifs'
data$SV = 'STR'; data$SV[data$period>6] = 'VNTR'
data$period6 = data$period; data$period6[data$period>6] = '>6'
data$tag <- factor(data$tag, levels=c(as.character(expression('\u2264 10 unique motifs')),
                                      as.character(expression('> 10 unique motifs'))))
data$period6 <- factor(data$period6, levels=c('1','2','3','4','5','6','>6'))

#data$ref[data$ref == 'GRCh38'] = 'GRCh38 (vamos v3.0 catalog expanded by TRExplore-v1.0)'
#data$ref[data$ref == 'T2T-CHM13'] = 'T2T-CHM13 (vamos v3.0 catalog)'


nrow(data[data$totalU>5,])
hist(data[data$totalU>5,]$info); mean(data[data$totalU>5,]$info); median(data[data$totalU>5,]$info)
hist(data[data$totalU>5,]$infoU); mean(data[data$totalU>5,]$infoU); median(data[data$totalU>5,]$infoU)


##################################################
# plot of TR size (bp) on the genome
##################################################
sizeData = data.frame(period=rep(c('1','2','3','4','5','6','>6'),2), ref=c(rep('GRCh38',7),rep('T2T-CHM13',7)))
sizeData$period <- factor(sizeData$period, levels=c('1','2','3','4','5','6','>6'))

counts = c(); sizes = c()
# obtain total loci counts and sizes (bp) by reference and period size
for (ref in c('GRCh38','T2T-CHM13')) {
  for (p in 1:6) {
    counts = c( counts, length(data$size[data$period==p & data$ref == ref]) )
    sizes = c( sizes, sum(data$size[data$period==p & data$ref == ref]) )
  }
  counts = c( counts, length(data$size[data$period>6 & data$ref == ref]) )
  sizes = c( sizes, sum(data$size[data$period>6 & data$ref == ref]) )
}
sizeData$count = counts; sizeData$size = sizes; sizeData$prop = 0

sizeData$prop[sizeData$ref=='GRCh38'] = round(sizeData$size[sizeData$ref=='GRCh38'] / 30882698.32, 2)
sizeData$prop[sizeData$ref=='T2T-CHM13'] = round(sizeData$size[sizeData$ref=='T2T-CHM13'] / 31172755.01, 2)
sizeData$percentage = paste0(sizeData$prop,'%')

# summary statistics
nrow(data[data$ref == 'GRCh38',]); sum(data[data$ref == 'GRCh38',]$size)
nrow(data[data$ref == 'T2T-CHM13',]); sum(data[data$ref == 'T2T-CHM13',]$size)
sum(sizeData$prop[sizeData$ref=='GRCh38']) # GRCh38 total size in percentage
sum(sizeData$count[sizeData$ref=='GRCh38']) # GRCh38 total count

sum(sizeData$prop[sizeData$ref=='T2T-CHM13']) # T2T total size in percentage
sum(sizeData$count[sizeData$ref=='T2T-CHM13']) # T2T total count

# plot
colors = c('#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')

png(file='set_sizes.png', width=2500, height=2000, res=300)
ggplot(data=data) + theme_classic() +
  geom_histogram(aes(x=period6, fill=period6), stat='count', position=position_dodge()) + # histogram
  geom_text(data=sizeData, aes(x=period, y=count, label=count), cex=6, vjust=-1) + # numbers for counts
  geom_point(data=sizeData, aes(x=period, y=size, color=period), cex=3) + # points for size
  #geom_text_repel(data=sizeData, aes(x=period, y=size, label=percentage), cex=6) + # numbers for proportion
  geom_text(data=sizeData, aes(x=period, y=size, label=percentage), nudge_y=0.3, cex=6) +
  theme(text = element_text(size = 18), legend.position='none') + 
  theme(panel.grid.major = element_line(colour = "white", linewidth = 0.2)) +
  theme(panel.grid.minor = element_line(colour = "white", linewidth = 0.5)) +
  facet_wrap(~ref, strip.position='bottom', scales='free_x', nrow=2) + # set multiple group labels on x axis
  scale_fill_manual(name='', values=colors) + # manually set color scheme for bar fill
  scale_color_manual(name='', values=colors) + # manually set color scheme for dot color
  coord_cartesian(ylim=c(10000,100000000), expand=TRUE, clip='on') +
  scale_y_log10('Frequency (bars)',
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                sec.axis = sec_axis(~ . * 1, name = 'Total size in bp (dots)',
                                    breaks = trans_breaks("log10", function(x) 10^x),
                                    labels = trans_format("log10", math_format(10^.x)))) +
  xlab('Motif period')
  #annotate(geom='text', x=1.5, y=110000000,
  #         label='Numbers indicate counts or percentage of genome', hjust=0, size=7)
  #annotate(geom='text', x=1.5, y=55000000,
  #         label='Total: 922,780 loci, 2.13% of genome', hjust=0, size=7)
dev.off()




# plot only T2T
png(file='set_sizes_T2T.png', width=2000, height=1200, res=300)
ggplot(data=data[data$ref=='T2T-CHM13',]) + theme_classic() +
  geom_histogram(aes(x=period6, fill=period6), stat='count', position=position_dodge()) + # histogram
  geom_text(data=sizeData[sizeData$ref=='T2T-CHM13',], aes(x=period, y=count, label=count), cex=6, vjust=-1) + # numbers for counts
  geom_point(data=sizeData[sizeData$ref=='T2T-CHM13',], aes(x=period, y=size, color=period), cex=3) + # points for size
  #geom_text_repel(data=sizeData, aes(x=period, y=size, label=percentage), cex=6) + # numbers for proportion
  geom_text(data=sizeData[sizeData$ref=='T2T-CHM13',], aes(x=period, y=size, label=percentage), nudge_y=0.3, cex=6) +
  theme(text = element_text(size = 18), legend.position='none') + 
  #facet_wrap(~ref, strip.position='bottom', scales='free_x') + # set multiple group labels on x axis
  scale_fill_manual(name='', values=colors) + # manually set color scheme for bar fill
  scale_color_manual(name='', values=colors) + # manually set color scheme for dot color
  coord_cartesian(ylim=c(10000,100000000), expand=TRUE, clip='on') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  xlab('Motif period') + ylab('Frequency/Total size in bp')
#annotate(geom='text', x=1.5, y=110000000,
#         label='Numbers indicate counts or percentage of genome', hjust=0, size=7)
#annotate(geom='text', x=1.5, y=55000000,
#         label='Total: 922,780 loci, 2.13% of genome', hjust=0, size=7)
dev.off()


##################################################
# plot of remaining information after efficient selection
##################################################
plotData = data[,c('period','info','SV','tag','ref')]; colnames(plotData)[2] = 'prop'; plotData$type = 'All motifs'
temp = data[,c('period','infoU','SV','tag','ref')]; colnames(temp)[2] = 'prop'; temp$type = 'Unique motifs'
plotData = rbind(plotData,temp)

n_fun <- function(x){
  return(data.frame(y = 1.05, label = paste0("n = ",length(x)/2)))
}
nrow(plotData[plotData$ref=='T2T-CHM13',]) / 2

png(file='information_remained.png', width=1200, height=960, res=300)
ggplot(data=plotData[plotData$ref=='T2T-CHM13',]) + theme_classic() +
  geom_boxplot(aes(x=SV, y=prop, fill=type), cex=0.2, outlier.size = 0.1) +
  theme(text = element_text(size = 13), legend.position='top') +
  xlab('') + ylab('Proportion remained') +
  facet_wrap(~tag, strip.position='bottom', scales='free_x') +
  stat_summary(fun.data = n_fun, geom = 'text', aes(x=SV, y=prop), cex=3) +
  scale_fill_manual(name='', values=c('#7fc97f', '#386cb0'))
dev.off()

png(file='information_remained_exclude1.png', width=2500, height=1500, res=300)
ggplot(data=plotData[plotData$period!=1,]) + geom_boxplot(aes(x=xlab, y=prop, fill=type)) +
  theme_classic() + theme(text = element_text(size = 15), legend.position='top') +
  xlab('') + ylab('Proportion remained') +
  stat_summary(fun.data = n_fun, geom = 'text', aes(x=xlab, y=prop)) +
  scale_fill_manual(name='', values=c('#7fc97f', '#386cb0'))
dev.off()




##################################################
# genomic features of the expanded v3.0 catalog
##################################################
# read catalog info
catalog = read.csv('../vamosExpanded_v3.0_effMotifs-0.1_entropy.tsv', sep='\t', header=F)
colnames(catalog) = c('chr','start','end','motifs','source','type','period','consensus',
                      'consGC','allGC','op1','op2','segDup','coding','exon','transcribed','entropy')
catalog$key = paste0(catalog$chr,'_',catalog$start,'_',catalog$end)
catalog$TE = catalog$op1 > 0 | catalog$op2 > 0
catalog$segDup2 = catalog$segDup != '.'
catalog$coding2 = catalog$coding != '.'

library(VennDiagram)
setAll = catalog$key
setTE = catalog$key[catalog$TE]
setSegdup = catalog$key[catalog$segDup2]
setCoding = catalog$key[catalog$coding2]

set1 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set2 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set3 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")
set4 <- paste(rep("word_" , 200) , sample(c(1:1000) , 200 , replace=F) , sep="")

colors = c('#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')
venn.diagram(
  x = list(setTE, setSegdup, setCoding),
  category.names = c('Transposable\nelement', 'Segmental\nduplication', 'Coding'),
  #x = list(setAll, setTE, setSegdup, setCoding),
  #category.names = c('all', 'TE', 'segDup', 'coding'),
  filename = 'venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png", height = 2000, width = 2000, resolution = 300, compression = "lzw",
  
  # Circles
  lwd = 2, lty = 'blank', fill = colors[c(1,2,4)],
  
  # Numbers
  cex = 1.8, fontfamily = "sans", #fontface = "bold",
  
  # Set names
  cat.cex = 1.8, cat.fontface = "bold", cat.default.pos = "outer",
  cat.pos = c(-20, 20, 135),
  cat.dist = c(0.085, 0.085, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)





################## old
##########
# plot of motif set purity comparing old GRCh38 and new GRCh38
##########

purity = read.csv('purityOld2New.tsv', sep='\t', header=T)

purity$period6 = purity$period; purity$period6[purity$period6 > 6] = '>6'
purity$ratioOri = purity$oldOri / purity$newOri; purity$ratioEff = purity$oldEff / purity$newEff
purity$ratioOri[purity$ratioOri > 5] = 5 # cap large values at 5
purity$ratioEff[purity$ratioEff > 5] = 5 # cap large values at 5
purity$ratioOri[is.na(purity$ratioOri)] = 1 # replace 0/0 case by 1
purity$ratioEff[is.na(purity$ratioEff)] = 1 # replace 0/0 case by 1

purity$period6 = factor(purity$period6, levels=c('1','2','3','4','5','6','>6'))

plotData = purity[,c('startO','period6','ratioOri')]; colnames(plotData)[3] = 'ratio'
plotData$group = 'Original set'
temp = purity[,c('startO','period6','ratioEff')]; colnames(temp)[3] = 'ratio'
temp$group = 'Efficient set (q=0.1)'
plotData = rbind(plotData,temp)

plotData$group = factor(plotData$group, levels=c('Original set', 'Efficient set (q=0.1)'))

png(file='purity_ratio.png', width=2500, height=1500, res=300)
ggplot(data=plotData[plotData$ratio != 0,]) + theme_classic() +
  geom_boxplot(aes(x=period6, y=ratio, fill=group)) +
  geom_hline(yintercept=1, linetype='dashed', color = "red") +
  xlab('Motif period in new set') + ylab('Ratio of average motif distance (old/new)') +
  theme(text = element_text(size = 15), legend.position = 'top') +
  scale_fill_manual(name='', values=c('#7fc97f', '#386cb0'))
dev.off()

groups = c('Old set higher purity', 'Balanced', 'New set higher purity')

plotWhichBetter = function(groups) {
  period6 = c('1','2','3','4','5','6','>6')
  
  plotData = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(plotData) = c('period6', 'props', 'group')
  for (p in period6) {
    periodData = purity[purity$period6 == p,]; total = nrow(periodData)
    temp = data.frame( period6=rep(p,3), props=1:3, group=groups )
    
    temp$props[1] = nrow(periodData[periodData$ratioOri < 0.5,]) / total
    temp$props[2] = nrow(periodData[periodData$ratioOri >= 0.5 & periodData$ratioOri <= 2,]) / total
    temp$props[3] = nrow(periodData[periodData$ratioOri > 2,]) / total
    
    plotData = rbind(plotData,temp)
  }
  
  plotData$period6 = factor(plotData$period6, levels=period6)
  plotData$group = factor(plotData$group, levels=groups)
  
  colors = c('#d53e4f','#ffffbf','#3288bd')
  plot1 = ggplot(data=plotData, aes(x=period6, y=props, fill=group)) + theme_classic() +
    geom_bar(stat="identity") + xlab('Motif period in new set') + ylab("Proportion") +
    scale_fill_manual('', values=colors) +
    theme(text = element_text(size = 30),
          legend.text = element_text(size=20), legend.position = 'top')
  
  return(plot1)
}

plot1 = plotWhichBetter(c('Ratio (old/new) < 0.5', 'Balanced', 'Ratio (old/new) > 2'))
plot1








