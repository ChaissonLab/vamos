import sys
import datetime
import logging

import lib.general as general

logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.DEBUG)


def parseSingleVcfOneLine(line, includeRS):

    chrom,pos,id,ref,alt,qual,filter,info,form,gt = line.strip().split('\t')
    if 'ALTANNO_H2' in info:
        end,ru,_,altH1,_,altH2,_ = info.split(';')[:-1]
    else:
        end,ru,_,altH1,_ = info.split(';')[:-1]
        altH2 = 'ALTANNO_H2=.'
    altH1 = altH1.split('=')[1]
    altH2 = altH2.split('=')[1]
    form = form.split(':')
    if includeRS and 'RS' not in form:
        logging.error(f'No "RS" field found in the input vcf(s)!')
        sys.exit(1)
    gt1,gt2 = gt.split(':')[form.index('GT')].split('/')
    value1 = [[altH1, altH2][int(gt1)-1]]
    value2 = [[altH1, altH2][int(gt2)-1]]
    if includeRS:
        rs1,rs2 = gt.split(':')[form.index('RS')].split(',')
        value1.append(rs1)
        value2.append(rs2)

    constant = [chrom, pos, id, ref, alt, qual, filter, end, ru]

    return(chrom, pos, constant, value1, value2)


# read header lines, ordered chrom list, and sample id for one vamos vcf
def advanceToChromForOneVcf(fileObject):
    meta, contigInfo = [], {}
    for line in fileObject:
        if line.startswith('##'):
            if line.startswith('##contig='):
                contigInfo[line.split(',')[0].split('ID=')[1]] = line
            if line.startswith('##fileformat') or line.startswith('##source'):
                meta.append(line)
        elif line.startswith('#CHROM'): # reach the "#Chrom" line
            sampleID = line.strip().split('\t')[-1]
            return(meta, contigInfo, sampleID)
        else:
            sys.exit('Error on input vcf header!')

# obtain the overall header, ordered chrom list, and sample id of all vamos vcf
def getHeader(headerInfo, pairingInfo, chromOrders, includeRS):
    firstChrom, ordering, contigInfoOverall = [], {}, {}
    sampleIDs = [ h[2] for h in headerInfo ]
    for meta,contigInfo,sampleID in headerInfo:
        chroms = list(contigInfo.keys())
        contigInfoOverall.update(contigInfo)
        if len(chroms) > 0: firstChrom.append(chroms[0]) # store the "first" chrom of all samples
        for i in range(len(chroms)-1): # count the chrom orders in all samples
            if chroms[i] not in ordering: ordering[chroms[i]] = {}
            if chroms[i+1] not in ordering[chroms[i]]:
                ordering[chroms[i]][chroms[i+1]] = 0
            ordering[chroms[i]][chroms[i+1]] += 1

    '''
    # find the most frequent "first" chrom in all samples
    firstChrom = max(set(firstChrom), key=firstChrom.count)

    # find the best chrom order in all samples
    for chrom,chromDict in ordering.items():
        ordering[chrom] = max(chromDict, key=chromDict.get)

    # config the final chrom list
    chroms = [firstChrom]
    currentChrom, nextChrom = firstChrom, ordering[firstChrom]
    while True:
        chroms.append(ordering[currentChrom])
        currentChrom = nextChrom
        if currentChrom not in ordering: break # currentChrom is the last chrom
        nextChrom = ordering[currentChrom]
    '''
    chroms = chromOrders.split(',')
    for chrom in chroms:
        if chrom not in contigInfoOverall:
            sys.exit(f'No input vcf headers contain information of {chrom}. \
                    \nPlease double check input vcfs or remove this chromosome from the "-c/--chromOrders" option.')

    logging.info(f"configured chromosome orders: {','.join(chroms)}")

    # config the output sampleIDs
    sampleIDsOut, dipSampleIDs = [], []
    for i,pair in enumerate(pairingInfo):
        if pair == 'Dip':
            sampleIDsOut.append(sampleIDs[i])
            dipSampleIDs.append(sampleIDs[i])
        elif pair == 'Hap1':
            sampleIDsOut.append(f'{sampleIDs[i]}/{sampleIDs[i+1]}')
            dipSampleIDs.append(f'{sampleIDs[i]}/{sampleIDs[i+1]}')
            dipSampleIDs.append(f'{sampleIDs[i]}/{sampleIDs[i+1]}')
        else:
            pass

    logging.info(f"configured diploid sample IDs: {','.join(sampleIDsOut)}")

    # config the final header
    header = headerInfo[0][0] + [contigInfoOverall[c] for c in chroms]
    header.append('##ALT=<ID=TR,Description="Allele comprised of TR repeat units">\n')
    header.append('##FILTER=<ID=PASS,Description="All filters passed">\n')
    header.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    header.append('##INFO=<ID=RU,Number=1,Type=String,Description="Comma separated motif sequences list in the reference orientation">\n')
    header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    if includeRS:
        header.append('##INFO=<ID=ALTNT,Number=A,Type=String,Description="Raw nucleotide sequence for all alleles, separated by ,">\n')
    header.append('##INFO=<ID=ALTANNO,Number=A,Type=String,Description="Motif representation for all alleles, separated by ,">\n')
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype by motif annotations, index as orders in the ALTANNO field">\n')
    if includeRS:
        header.append('##FORMAT=<ID=RS,Number=1,Type=String,Description="Genotype by raw nucleotide sequence, index as orders in the ALTNT field">\n')
    header.append('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT'] + sampleIDsOut) +'\n')

    return(header, chroms, sampleIDs, dipSampleIDs)


# check if need to read the next line for one vcf
def advanceLineForOneVcf(curChrom, curStart, minChrom, minStart, parsedLine, vcfFileObject, includeRS):
    if curChrom is None or (curChrom == minChrom and curStart == minStart):
        nextLine = vcfFileObject.readline()
        if nextLine is None or nextLine == '':
            return None
        else:
            return parseSingleVcfOneLine(nextLine, includeRS)
    else:
        return parsedLine


def combineVcfs(inVcfs, outVcf, chromOrders, includeRS):

    # config input vcf files and the hap/dip pairing information
    inVcfsInfo = [ l.strip().split(',') for l in open(inVcfs) ]
    vcfs, pairingInfo = [], []
    for vcf in inVcfsInfo:
        if len(vcf) == 1:
            logging.info(f'input single diploid vcf: {vcf[0]}')
            pairingInfo += ['Dip']
        elif len(vcf) == 2:
            logging.info(f'input pair of haploid vcfs: {vcf[0]} & {vcf[1]}')
            pairingInfo += ['Hap1', 'Hap2']
        else:
            sys.exit('Error on input! Please check the input instructions.')
        vcfs += vcf
    vcfFileObjects = [ open(f) for f in vcfs ]

    # obtain the overall header, ordered chrom list, and sample IDs
    headerInfo = [ advanceToChromForOneVcf(fo) for fo in vcfFileObjects ]
    header, chroms, sampleIDs, dipSampleIDs = getHeader(headerInfo, pairingInfo, chromOrders, includeRS)

    out = open(outVcf, 'w')
    # output the header
    out.write(''.join(header))

    nVcfs, minChrom, minStart = len(vcfs), chroms[0], None
    curChroms = [ None for i in range(0,nVcfs) ]
    curStarts = [ None for i in range(0,nVcfs) ]
    parsedLines = [ None for i in range(0,nVcfs) ]
    logging.info(f'handling chromosome: {minChrom}')

    while True:
        # update the parsed lines of each vcf
        parsedLines = [advanceLineForOneVcf(curChroms[i], curStarts[i], minChrom, minStart, parsedLines[i], \
                                            vcfFileObjects[i], includeRS) for i in range(0,nVcfs)]
        # break if no more remaining lines in any vcf
        if sum([ l is None for l in parsedLines ]) == nVcfs: break

        # update the current chrom and starts of each vcf
        curChroms = [ parsedLines[i][0] if parsedLines[i] else None for i in range(nVcfs) ]
        curStarts = [ parsedLines[i][1] if parsedLines[i] else None for i in range(nVcfs) ]
        # update the minChrom and minStart
        temp = minChrom
        minChrom = min([c for c in curChroms if c], key=chroms.index)
        minStart = str(min([int(s) for i,s in enumerate(curStarts) if curChroms[i]==minChrom]))
        # switch to a new chrom detected
        if minChrom != temp: logging.info(f'handling chromosome: {minChrom}')

        # determine which vcf to read (minChrom, minStart is the locus to read)
        isMatching = [ curChroms[i] == minChrom and curStarts[i] == minStart for i in range(0,nVcfs) ]

        # read all vcf alleles and get the genotype index of each allele
        allelesGT, allelesRS, values = [], [], []
        for i,sampleID in enumerate(sampleIDs):
            # the vcf has matched locus with "minChrom, minStart"
            if isMatching[i]:
                chrom, start, constant, value1, value2 = parsedLines[i]
                alleleGT1, alleleGT2 = value1[0], value2[0]
                # for one diploid vcf
                if pairingInfo[i] == 'Dip':
                    if alleleGT1 not in allelesGT: allelesGT.append(alleleGT1)
                    if alleleGT2 not in allelesGT: allelesGT.append(alleleGT2)
                    if includeRS: # if RS value should be included
                        alleleRS1, alleleRS2 = value1[1], value2[1]
                        if alleleRS1 not in allelesRS: allelesRS.append(alleleRS1)
                        if alleleRS2 not in allelesRS: allelesRS.append(alleleRS2)
                        rs1, rs2 = allelesGT.index(alleleRS1)+1, allelesGT.index(alleleRS2)+1
                # for hap1 of a pair of hap1/hap2 vcfs
                elif pairingInfo[i] == 'Hap1':
                    if alleleGT1 not in allelesGT: allelesGT.append(alleleGT1)
                    gt1 = allelesGT.index(alleleGT1)+1
                    if includeRS: # if RS value should be included
                        alleleRS1 = value1[1]
                        if alleleRS1 not in allelesRS: allelesRS.append(alleleRS1)
                        rs1 = allelesRS.index(alleleRS1)+1
                # for hap2 of a pair of hap1/hap2 vcfs
                else:
                    if alleleGT2 not in allelesGT: allelesGT.append(alleleGT2)
                    gt2 = allelesGT.index(alleleGT2)+1
                    if includeRS: # if RS value should be included
                        alleleRS2 = value2[1]
                        if alleleRS2 not in allelesRS: allelesRS.append(alleleRS2)
                        rs2 = allelesRS.index(alleleRS2)+1
            # the vcf has no matched locus with "minChrom, minStart"
            else:
                # for one diploid vcf
                if pairingInfo[i] == 'Dip':
                    gt1, gt2 = '.', '.'
                    if includeRS: # if RS value should be included
                        rs1, rs2 = '.', '.'
                # for hap1 of a pair of hap1/hap2 vcfs
                elif pairingInfo[i] == 'Hap1':
                    gt1 = '.'
                    if includeRS: # if RS value should be included
                        rs1 = '.'
                # for hap2 of a pair of hap1/hap2 vcfs
                else:
                    gt2 = '.'
                    if includeRS: # if RS value should be included
                        rs2 = '.'
            # record the diploid genotypes
            if pairingInfo[i] in ['Dip','Hap2']:
                if includeRS:
                    values.append(f'{gt1}/{gt2}:{rs1}/{rs2}')
                else:
                    values.append(f'{gt1}/{gt2}')
        allelesGT = ','.join([ a.replace(',','-') for i,a in enumerate(allelesGT) ])
        info = f"{constant[7]};{constant[8]};SVTYPE=VNTR;ALTANNO={allelesGT}"
        if includeRS:
            allelesRS = ','.join(allelesRS)
            info = f"{constant[7]};{constant[8]};SVTYPE=VNTR;ALTNT={allelesRS};ALTANNO={allelesGT}"
            out.write('\t'.join(constant[:7] + [info,'GT:RS'] + values) + '\n')
        else:
            out.write('\t'.join(constant[:7] + [info,'GT'] + values) + '\n')

    out.close()
    for fo in vcfFileObjects: fo.close()

