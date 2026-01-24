library(ggplot2)

pcaFunc = function(data, info, batches) {
  # drop rows with any NA
  data = na.omit(data)
  # remove constant loci
  delIndex = which( data[,4] == rowMeans(data[,4:ncol(data)]) )
  keepIndex = setdiff(1:nrow(data), delIndex)
  data = data[keepIndex,]
  
  # prepare the pca data
  pcaData = data[,4:ncol(data)]
  # organize data for pca: features as cols, subjects as rows
  pcaData = t(pcaData)
  
  out = prcomp(pcaData, center = T, scale. = T)
  print( summary(out)$importance[,1:10] )
  
  plotData = as.data.frame(out$x)[,1:10]
  plotData$population = c(info[4:nrow(info),2])
  plotData$batch = c(info[4:nrow(info),1])
  plotData$batch = factor(plotData$batch, levels=batches)
  
  return(plotData)
}
color = c('AMR'='#e41a1c','EUR'='#377eb8','EAS'='#4daf4a','SAS'='#984ea3','AFR'='#ff7f00','MDE'='#4d4d4d')
color = c('AMR'='#e41a1c','EUR'='#b3cde3','EAS'='#4daf4a','SAS'='#984ea3','AFR'='#fecc5c','MDE'='#4d4d4d')
pch = c('CPC'=21, 'HGSVC2'=22, 'HGSVC3'=23, 'HPRC'=24)


##########
# cluster - GRCh38 full data, PC1/PC2
##########
prefix = 'samples_rmDup_GRCh38_q-0.1.topCount.rmHomo.rmNAConst'
info = t(read.csv(paste0(prefix,'.tsv'), sep='\t', header=T, nrows=2))
data = read.csv(paste0(prefix,'.tsv'), sep='\t', header=F, skip=3)
colnames(data) = rownames(info)

plotData = pcaFunc(data, info,
                   c('APR-phase1', 'CPC-phase1', 'HGSVC-phase2', 'HGSVC-phase3',
                     'HPRC-phase1', 'HPRC-phase2', 'JSA-phase1'))

plotData$PC1 = -plotData$PC1
plotData$PC2 = -plotData$PC2


png(file=paste0(prefix,'.png'), width=2500, height=2000, res=300)
ggplot(data=plotData[plotData$batch != 'HPRC-phase1' & plotData$batch != 'HGSVC-phase2',]) + theme_classic() +
  geom_point(aes(x=PC1, y=PC2, color=population, pch=batch), size=3) +
  scale_color_manual(values=color) + #scale_shape_manual(values=pch) +
  guides(color = guide_legend(ncol = 2)) + # This line changes the legend to two columns
  theme(text=element_text(size=30), plot.title=element_text(hjust=0.5),
        panel.grid.minor = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.25),
        panel.grid.major = element_line(color = "lightgray", linetype = "dotted", linewidth = 0.5),
        axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        legend.position=c(0.55,0.75), legend.title=element_blank()) #+ xlim(-50, 150) + ylim(-50, 25)
#geom_hline(yintercept = 0, colour="red") +
#geom_vline(xintercept = 0, colour="blue") +
#coord_fixed(ratio = 1)
dev.off()









##########
# cluster - T2T full data, PC1/PC2
##########
prefix = 'q-0.1_full.topCount.rmNAConst.T2T'
info = t(read.csv(paste0(prefix,'.tsv'), sep='\t', header=T, nrows=4))
data = read.csv(paste0(prefix,'.tsv'), sep='\t', header=F, skip=5)
colnames(data) = rownames(info)
info[,2][info[,2] == 'HPRC2'] = 'HPRC'

plotData = pcaFunc(data, info, c('CPC','HGSVC3','HPRC','HGSVC2','APRC'))

plotData$PC1 = -plotData$PC1
plotData$PC2 = -plotData$PC2

png(file=paste0(prefix,'.png'), width=2500, height=2000, res=300)
ggplot(data=plotData) + theme_classic() +
  geom_point(aes(x=PC1, y=PC2, color=population, pch=batch), size=3) +
  scale_color_manual(values=color) + #scale_shape_manual(values=pch) +
  #+ xlim(-50, 150) + ylim(-50, 25) +
  theme(text=element_text(size=20), plot.title=element_text(hjust=0.5),
        legend.position.inside=c(0.25,0.25), legend.title=element_blank())
#geom_hline(yintercept = 0, colour="red") +
#geom_vline(xintercept = 0, colour="blue") +
#coord_fixed(ratio = 1)
dev.off()






# organize data by AFR/non-AFR
temp = data[,1:3]
samples = c('chr','start','end')
for (i in 4:nrow(info)) {
  sample = rownames(info)[i]
  if (info[i,1] == 'AFR') {
    samples = c(samples, sample)
    temp = cbind(temp, data[,sample])
  }
}
for (i in 4:nrow(info)) {
  sample = rownames(info)[i]
  if (info[i,1] != 'AFR') {
    samples = c(samples, sample)
    temp = cbind(temp, data[,sample])
  }
}
colnames(temp) = samples
temp[is.na(temp)] <- 0 # fill NA as 0
data = temp

# t-test
data$p = 0
for (i in 1:nrow(data)) {
  x = data[i,4:63]
  y = data[i,64:ncol(data)-1]
  data$p[i] = t.test(x,y)$p.value
}

select = data[data$p < 0.00000005,]
select = select[order(select$p, decreasing = TRUE),]

# normalize data
temp = t(apply(data[,4:(ncol(data)-1)], 1, function(x)(x-min(x))/(max(x)-min(x))))
#temp = t(apply(data[,4:(ncol(data)-1)], 1, rescale, to=c(-1,1)))

data = cbind(data[,c(1:3,ncol(data))],temp)

# prepare the heatmap data for chr1
chr = 'chr21'
coors = data[data$chr==chr,2]
samples = colnames(data)[c(5:ncol(data))]
#samples = colnames(data)[c(5:217)]
plotData = expand.grid(Coors=as.factor(coors), Samples=samples)
temp = c()
for (s in samples) temp = c(temp,data[data$chr==chr,s])
plotData$length = round(temp,2)

# Heatmap 
ggplot(plotData, aes(Samples, Coors, fill=length)) + geom_tile() +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  #theme(legend.text = element_text(20), legend.title = element_text(size=20)) +
  #scale_fill_gradientn(colors = hcl.colors(20, "RdYlGn"))
  #scale_fill_gradient(low = "#ffff99", high = "#f0027f")
  scale_fill_gradient(low = "white", high = "red")



