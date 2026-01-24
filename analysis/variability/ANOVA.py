import enum
import sys
import datetime
import logging
import statistics
import numpy as np

from scipy import stats

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)

vamosFeaturesFile = sys.argv[1]
outFile = sys.argv[2]

out = open(outFile, 'w')
# read in the vamos feature file
with open(vamosFeaturesFile) as f:
    for i,line in enumerate(f):
        fields = line.strip().split('\t')
        if i == 0:
            samples = fields[3:]
            continue
        if i == 1:
            batches = fields[3:]
            continue
        if i == 2:
            populations = fields[3:]
            #populationSet = list(set(populations))
            populationSet = ['EUR','AMR','AFR','SAS','EAS','MDE']
            out.write('\t'.join(map(str,['chr','start','end','grandMean']+populationSet+\
                ['withinGroupVariance','betweenGroupVariance','withinGroupVariance_mean','betweenGroupVariance_mean','df','f','p','effectSize'])) +'\n')
            continue

        chr,start,end = fields[:3]
        features = fields[3:]
        featuresSelected = []

        if i % 100000 == 0: logging.info(f'{i} lines processed')

        # group data by population group
        dataDict = { population:[] for population in populationSet }
        for j,batch in enumerate(batches):
            #if batch != 'HPRC-phase2': continue
            sample,population,feature = samples[j],populations[j],features[j]
            if feature == 'NA':
                #dataDict[population].append(np.nan)
                pass
            else:
                dataDict[population].append(int(feature))
                featuresSelected.append(int(feature))

        # skip constant loci
        if len(set([ a for a in featuresSelected ])) == 1: continue

        # calculate the within-group variations
        groupMeans, Ns, withinGroupVariances = [], [], []
        for population,data in dataDict.items():
            skip,withinGroupVariance = '',0

            # only work on loci having >=5 annotated genomes in each population
            if len(data) < 5:
                skip = 'skip'
                break

            groupMean = np.mean(np.array(data))
            Ns.append(len(data))
            groupMeans.append(groupMean)
            for j,d in enumerate(data):
                withinGroupVariance += (d-groupMean) * (d-groupMean)
            withinGroupVariances.append(withinGroupVariance)

        # only work on loci having >5 annotated genomes in each population
        if skip == 'skip': continue
        withinGroupVariance = sum(withinGroupVariances)

        # calculate the between-group variations
        dataAll, anovaData = [], []
        for population,data in dataDict.items():
            dataAll += data
            anovaData.append(data)
        grandMean = np.mean(np.array(dataAll))
        betweenGroupVariance = 0
        for j,groupMean in enumerate(groupMeans):
            betweenGroupVariance += Ns[j]*(groupMean-grandMean)*(groupMean-grandMean)
        # calculate adjusted variance by the degree of freedom
        withinGroupVariance_mean = withinGroupVariance / (len(dataAll)-len(populationSet))
        betweenGroupVariance_mean = betweenGroupVariance / (len(populationSet)-1)

        # one-way anova
        f, p = stats.f_oneway(anovaData[0], anovaData[1], anovaData[2], anovaData[3], anovaData[4], anovaData[5])
        df = f'{len(dataDict)-1},{len(dataAll)-len(dataDict)}'
        # effect size
        effectSize = max(groupMeans) / min(groupMeans)
        groupMeans = [round(d,2) for d in groupMeans]
        withinGroupVariance = round(withinGroupVariance, 2)
        betweenGroupVariance = round(betweenGroupVariance, 2)
        withinGroupVariance_mean = round(withinGroupVariance_mean, 2)
        betweenGroupVariance_mean = round(betweenGroupVariance_mean, 2)

        f = round(f, 2)
        effectSize = round(effectSize, 2)

        if effectSize == 1: continue
        out.write('\t'.join(map(str,[chr,start,end,grandMean]+groupMeans+\
            [withinGroupVariance,betweenGroupVariance,withinGroupVariance_mean,betweenGroupVariance_mean,df,f,p,effectSize])) +'\n')

out.close()

logging.info('End of Program\n')
