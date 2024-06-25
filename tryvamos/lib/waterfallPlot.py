# -*- coding: UTF-8 -*-
import os
import logging

from lib.tr import TR
import lib.general as general


logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', \
        datefmt='%Y-%m-%d %H:%M:%S', \
        level=logging.INFO)


def plot(inVCF, outDir, useLoci, sort):

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

            tr = TR()
            tr.parseDiploidVCFOneLine(samples, line)

            # skip unwanted loci
            if candidateLoci:
                if tr.chr not in candidateLoci:
                    continue
                else:
                    if (tr.start,tr.end) not in candidateLoci[tr.chr]:
                        continue

            tr.appendLength(sort)
            tr.heat(f'{outDir}/{tr.chr}_{tr.start}-{tr.end}.png')

