import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import colorcet as cc
import lib.general as general
import edlib

class TR:
    """class for one TR locus.

    Attributes:
        chr (str): chromosome of TR.
        start (str): start coordinate of TR.
        end (str): end coordinate of TR.
        motifsFull (dict[str:int]): indexing table of the full list of motifs
        motifsUsed (dict[str:int]): indexing table of the used list of motifs
        annosByUsed (dict[str:list[str]]): annotations of all samples (by index
            of the used motif list, unannotated allele as ['.'], full deletion
            allele as [])
        annosByFull (dict[str:list[str]]): annotations of all samples (by index
            of the full motif list, unannotated allele as ['.'], full deletion
            allele as [])
        constant (bool): if the locus is constant across all samples
        maxLen (int): maxmimum annotation length of all alleles
        space (bool): if annotation of any sample is extended to the max length
        edMat (dict[tuple[(str,str)]:int]): edit distance (of annotation string)
            matrix for all samples

    """

    def __init__(self):
        """the __init__ method
        """

        # parsed in method "parseDiploidVCFOneLine"
        self.chr = ''
        self.start = ''
        self.end = ''
        self.motifsUsed = {}
        self.annosByUsed = {}
        self.annosByFull = {}
        self.constant = False
        self.maxLen = 0
        # parsed in method "readMotifsFull"
        self.motifsFull = {}
        # changed in method "appendLength"
        self.space = False
        # assigned in method "getEditDistance"
        self.edMat = {}


    def readMotifsFull(self, motifsFull:list[str]):
        """Encode the full motif list with numbering

        Args:
            motifs (list[str]): input list of full motifs
        """
        for i,m in enumerate(motifsFull): self.motifsFull[m] = i


    def parseDiploidVCFOneLine(self, samples:list[str], line:str):
        """parse one entry line from a diploid vamos vcf

        Args:
            samples (list[str]): ordered list of all sample names for gt parsing
            line (str): one line from the combined vamos vcf
        """

        fields = line.strip().split()
        self.chr,self.start,_,_,_,_,_,info,_ = fields[:9]
        end,ru,_,annos = info.split(';')
        self.end = end.split('=')[1]
        ru = ru.split('=')[1].split(',')
        # encode the list of motifs used for annotation with numbering
        for i,m in enumerate(ru): self.motifsUsed[m] = i

        gts = [ gt.split('/') for gt in fields[9:] ]
        alleles = annos.split('=')[1].split(',')

        # check if the locus is constant over all samples
        if alleles.count(alleles[0]) == len(alleles): self.constant = True

        # parse the h1/h2 allele
        for i,sample in enumerate(samples):

            if gts[i][0] == '.':
                temp = ['.']
            elif gts[i][0] == 'DEL':
                temp = []
            else:
                temp = alleles[int(gts[i][0])-1].split('-')

            self.annosByUsed[sample+'_h1'] = temp
            if self.annosByFull:
                if temp != ['.']:
                    temp = [ self.motifsFull[ru[int(t)]] for t in temp ]
                self.annosByFull[sample+'_h1'] = temp

            if gts[i][1] == '.':
                temp = ['.']
            elif gts[i][1] == 'DEL':
                temp = []
            else:
                temp = alleles[int(gts[i][1])-1].split('-')

            self.annosByUsed[sample+'_h2'] = temp
            if self.annosByFull:
                if temp != ['.']:
                    temp = [ self.motifsFull[ru[int(t)]] for t in temp ]
                self.annosByFull[sample+'_h2'] = temp

        # calculate the max annotation length among all alleles
        self.maxLen = max([ len(a) for s,a in self.annosByUsed.items() ])


    def getEditDistance(self):
        """calculate pairwise edit distance of all annotation strings
        apply ascii encoding if <=256 unique motifs and edlib for speeding up
        """

        checkDict = {}
        for i,s in enumerate(self.annosByUsed.keys()):
            for j,t in enumerate(self.annosByUsed.keys()):
                if i >= j: continue
                anno1, anno2 = self.annosByUsed[s], self.annosByUsed[t]
                temp1, temp2 = '-'.join(anno1), '-'.join(anno2)
                if anno1 in [['.'],[]] or anno2 in [['.'],[]]:
                    ed = 'NULL'
                else:
                    if (temp1, temp2) in checkDict:
                        ed = checkDict[(temp1, temp2)]
                    else:
                        if len(self.motifsUsed) <= 256:
                            anno1 = [chr(int(a)) for a in anno1]
                            anno2 = [chr(int(a)) for a in anno2]
                            ed = edlib.align(anno1, anno2)['editDistance']
                        else:
                            ed = general.alignGlobal(anno1,anno2,True,0,1,1)
                        checkDict[(temp1, temp2)] = ed
                        checkDict[(temp2, temp1)] = ed
                self.edMat[(s,t)] = ed
                self.edMat[(t,s)] = ed

    def appendLength(self, sort=True):
        """Extend all annotations to match length of the longest annotation
        """

        # sort by length
        if sort:
            temp = sorted(self.annosByUsed.items(), \
                    key=lambda item: len(item[1]), reverse=True)
            self.annosByUsed = {k: v for k, v in temp}

        nMotifs = len(self.motifsUsed)
        whiteMotif = str(nMotifs+1)
        maxLen = self.maxLen + 1

        # extend alleles to maxLen
        for sample in self.annosByUsed.keys():
            if len(self.annosByUsed[sample]) < maxLen:
                self.annosByUsed[sample].extend( [whiteMotif]*(maxLen - \
                                            len(self.annosByUsed[sample])) )
                self.space = True


    def heat(self, outPlot:str):
        """generate the allele plot

        Args:
            outPlot (str): path for the output plot
        """

        samples = list(self.annosByUsed.keys())
        # remove uncovered samples
        samples = [ s for s in samples if '.' not in self.annosByUsed[s] ]
        annos = np.array([ [int(m) for m in self.annosByUsed[s]] \
                                   for s in samples ])

        nMotifs = len(self.motifsUsed)
        whiteMotif = str(nMotifs+1)
        maxLen = self.maxLen + 1

        # set up and output color code
        cmap = sns.color_palette(cc.glasbey_hv, nMotifs)
        cmap.append((1,1,1))
        out = open(outPlot.replace('.png','.tsv'), 'w')
        for i,m in enumerate(self.motifsUsed):
            out.write('\t'.join([cmap.as_hex()[i],m])+'\n')
        out.close()

        numTicks = len(samples)
        yticks = np.linspace(0, len(samples) - 1, numTicks, dtype=int)
        yticklabels = [ samples[i] for i in yticks ]

        plt.figure(figsize=(15, 12))
        sns.set(font_scale=2)
        ax = sns.heatmap(annos, cmap=cmap, cbar=False)
        #ax = sns.heatmap(np_annp, cmap=cmap, yticklabels=yticklabels)
        #ax.set_yticklabels(yticklabels, rotation=0, fontsize="3")
        #print(yticklabels)
        plt.savefig(outPlot, bbox_inches='tight', dpi=300, format='png')
        plt.close()

