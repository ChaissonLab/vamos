import os
import logging

from lib.tr import TR
import lib.general as general


logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)


def plot(inVCF, outDir, useLoci, width, height, ylabel, sort):

    if useLoci:
        candidateLoci = general.parseBed(useLoci)
    else:
        candidateLoci = {}

    if not os.path.exists(outDir): os.makedirs(outDir)

    with open(inVCF) as f:
        for line in f:
            if line.startswith('##'): continue
            if line.startswith('#'):
                samples = line.strip().split()[9:]
                continue

            # skip chromosome(s) without candidate plotting loci
            chr = line.strip().split()[0]
            if chr not in candidateLoci: continue
            tr = TR()
            tr.parseDiploidVCFOneLine(samples, line)

            # skip loci not in candidate plotting loci list
            if candidateLoci:
                if (tr.start,tr.end) not in candidateLoci[tr.chr]:
                    continue

            logging.info(f'plotting locus: {tr.chr}:{tr.start}-{tr.end}')

            tr.appendLength(sort)
            tr.heat(f'{outDir}/{tr.chr}_{tr.start}-{tr.end}.png', width, height, ylabel)

