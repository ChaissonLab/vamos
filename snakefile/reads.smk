# -*- coding: utf-8 -*-

import os, sys, re

batchDict = {}
if not os.path.exists('analysis'): os.makedirs('analysis')
if not os.path.exists('log'): os.makedirs('log')
if not os.path.exists('status'): os.makedirs('status')

manifest = config['input']['manifest']
# config input manifest
with open(manifest) as f:
    for line in f:

        sample, seqFile = line.strip().split(',')
        batchDict[sample] = seqFile


rule all:
    input:
        ['analysis/split/split.done'],
        [f'analysis/{sample}/phase/{sample}.phase.done' for sample in batchDict.keys()],
        [f'analysis/{sample}/anno/{sample}.anno.done' for sample in batchDict.keys()]


# split genomic regions for phasing
rule split:
    input:
        inVNTRs = config['databases']['evntr'],
        inMotifs = config['databases']['emotif']
    output:
        done = 'analysis/split/split.done'
    params:
        python = config['softwares']['python'],
        split = config['tools']['split'],
        winSize = config['parameters']['winSize'],
        outDir = 'analysis/split',
        rule = 'split',
        status = 'status/all__run.log',
        cluster = config['cluster']['split'],
        stdout = 'log/split.o.%j',
        stderr = 'log/split.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {params.python} \
        {params.split} \
        --inVNTRs {input.inVNTRs} \
        --inMotifs {input.inMotifs} \
        --winSize {params.winSize} \
        --outDir {params.outDir}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

# longshot phasing
rule phase:
    input:
        bam = lambda wildcards: batchDict[wildcards.sample],
        splitDone = 'analysis/split/split.done'
    output:
        done = 'analysis/{sample}/phase/{sample}.phase.done'
    params:
        python = config['softwares']['python'],
        phase = config['tools']['phase'],
        winDir = 'analysis/split',
        longshot = config['softwares']['longshot'],
        ref = config['databases']['GRCh38'],
        samtools = config['softwares']['samtools'],
        outDir = 'analysis/{sample}/phase',
        rule = 'phase__{sample}',
        status = 'status/all__run.log',
        cluster = config['cluster']['phase'],
        stdout = 'log/phase_{sample}.o.%j',
        stderr = 'log/phase_{sample}.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {params.python} \
        {params.phase} \
        --inBam {input.bam} \
        --winDir {params.winDir} \
        --longshot {params.longshot} \
        --ref {params.ref} \
        --samtools {params.samtools} \
        --thread 16 \
        --outDir {params.outDir}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

# vamos annotation
rule anno:
    input:
        phaseDone = 'analysis/{sample}/phase/{sample}.phase.done'
    output:
        done = 'analysis/{sample}/anno/{sample}.anno.done'
    params:
        python = config['softwares']['python'],
        anno = config['tools']['anno'],
        sample = lambda wildcards: wildcards.sample,
        vamos = config['softwares']['vamos'],
        filter = config['parameters']['filter'],
        rule = 'anno__{sample}',
        status = 'status/all__run.log',
        cluster = config['cluster']['anno'],
        stdout = 'log/anno_{sample}.o.%j',
        stderr = 'log/anno_{sample}.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {params.python} \
        {params.anno} \
        --sample {params.sample} \
        --vamos {params.vamos} \
        --filter {params.filter}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''



