# -*- coding: utf-8 -*-

import os, sys, re

batchDict = {}
if not os.path.exists('analysis'): os.makedirs('analysis')
if not os.path.exists('log'): os.makedirs('log')
if not os.path.exists('status'): os.makedirs('status')


# check of input parameters
mode_of_analysis = config['parameter']['mode_of_analysis']
type_of_input = config['parameter']['type_of_input']

if mode_of_analysis not in ['raw', 'assembly']:
    sys.exit('Error: "mode_of_analysis" must be set as "raw" or "assembly"')

if type_of_input not in ['fasta', 'fastq', 'bam']:
    sys.exit('Error: "type_of_input" must be set as "fasta" or "fastq" or "bam"')

if type_of_input in ['fasta','fastq']:
    type_of_aligner = config['parameter']['type_of_aligner']
    if type_of_aligner not in ['mm2', 'lra']:
        sys.exit('Error: "type_of_aligner" must be set as "lra" or "mm2"')


# read input manifest
with open(config['input']['manifest']) as f:
    for line in f:

        sample, seq = line.strip().split(',')
        bam = 'analysis/%s/aln/%s.aln.bam' %(sample, sample)
        if not os.path.exists('analysis/%s/aln' %(sample)): os.makedirs('analysis/%s/aln' %(sample))
        if type_of_input == 'bam' and not os.path.exists(bam):
            os.symlink(seq, bam)
            os.symlink(seq+'.bai', bam+'.bai')
        batchDict[sample] = [seq, bam]


# set flow of analysis
if mode_of_analysis == 'raw':

    # analysis of raw sequencing reads, with alignment step
    if type_of_input in ['fasta', 'fastq']:
        rule all:
            input:
                ['analysis/split/split.done'],
                [f'analysis/{sample}/aln/{sample}.aln.done' for sample in batchDict.keys()],
                [f'analysis/{sample}/phase/{sample}.phase_with_aln.done' for sample in batchDict.keys()],
                [f'analysis/{sample}/anno/{sample}.anno.raw_with_aln.done' for sample in batchDict.keys()]

    # analysis of raw sequencing reads, without alignment step
    else:
        rule all:
            input:
                ['analysis/split/split.done'],
                [f'analysis/{sample}/phase/{sample}.phase_without_aln.done' for sample in batchDict.keys()],
                [f'analysis/{sample}/anno/{sample}.anno.raw_without_aln.done' for sample in batchDict.keys()]

else:

    # analysis of assembly contigs, with alignment step
    if type_of_input in ['fasta', 'fastq']:
        rule all:
            input:
                [f'analysis/{sample}/aln/{sample}.aln.done' for sample in batchDict.keys()],
                [f'analysis/{sample}/anno/{sample}.anno.ass_with_aln.done' for sample in batchDict.keys()]

    # analysis of assembly contigs, without alignment step
    else:
        rule all:
            input:
                [f'analysis/{sample}/anno/{sample}.anno.ass_without_aln.done' for sample in batchDict.keys()]


# alignment
rule aln:
    input:
        fa = lambda wildcards: batchDict[wildcards.sample][0]
    output:
        bam = 'analysis/{sample}/aln/{sample}.aln.bam',
        done = 'analysis/{sample}/aln/{sample}.aln.done'
    params:
        aligner = config['parameter']['type_of_aligner'],
        ref = config['database']['reference'],
        rule = 'aln__{sample}',
        status = 'status/all__run.log',
        cluster = config['cluster']['aln'],
        stdout = 'log/aln_{sample}.o.%j',
        stderr = 'log/aln_{sample}.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        if [ "{params.aligner}" = "lra" ]; then
            lra align \
            -CONTIG {params.ref} \
            {input.fa} \
            -t 16 -p s -H | samtools sort -@4 > {output.bam}
        else
            minimap2 \
            -a {params.ref} \
            {input.fa} \
            -t 16 | samtools sort -@4 > {output.bam}
        fi

        samtools index -@16 {output.bam}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

# split genomic regions for phasing
rule split:
    input:
        vntr = config['database']['vntr'],
    output:
        done = 'analysis/split/split.done'
    params:
        repo = config['parameter']['vamos_repo'],
        winSize = config['parameter']['window_size'],
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

        python \
        {params.repo}/snakefile/pyscript/split_for_phasing.py \
        --inVNTRs {input.vntr} \
        --winSize {params.winSize} \
        --outDir {params.outDir}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

# longshot phasing (with alignment)
rule phase_with_aln:
    input:
        bam = lambda wildcards: batchDict[wildcards.sample][1],
        alnDone = 'analysis/{sample}/aln/{sample}.aln.done',
        splitDone = 'analysis/split/split.done'
    output:
        done = 'analysis/{sample}/phase/{sample}.phase_with_aln.done'
    params:
        repo = config['parameter']['vamos_repo'],
        winDir = 'analysis/split',
        ref = config['database']['reference'],
        thread = 16,
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

        python \
        {params.repo}/snakefile/pyscript/phase.py \
        --inBam {input.bam} \
        --winDir {params.winDir} \
        --longshot longshot \
        --ref {params.ref} \
        --samtools samtools \
        --thread {params.thread} \
        --outDir {params.outDir}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

# longshot phasing (without alignment)
rule phase_without_aln:
    input:
        bam = lambda wildcards: batchDict[wildcards.sample][1],
        splitDone = 'analysis/split/split.done'
    output:
        done = 'analysis/{sample}/phase/{sample}.phase_without_aln.done'
    params:
        repo = config['parameter']['vamos_repo'],
        winDir = 'analysis/split',
        ref = config['database']['reference'],
        thread = 16,
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

        python \
        {params.repo}/snakefile/pyscript/phase.py \
        --inBam {input.bam} \
        --winDir {params.winDir} \
        --longshot longshot \
        --ref {params.ref} \
        --samtools samtools \
        --thread {params.thread} \
        --outDir {params.outDir}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''


# readwise annotation by efficient motifs (for raw sequencing reads, with alignment)
rule anno_raw_with_aln:
    input:
        phaseDone = 'analysis/{sample}/phase/{sample}.phase_with_aln.done'
    output:
        done = 'analysis/{sample}/anno/{sample}.anno.raw_with_aln.done'
    params:
        repo = config['parameter']['vamos_repo'],
        sample = lambda wildcards: wildcards.sample,
        thread = 16,
        filter = config['parameter']['min_depth'],
        rule = 'anno__{sample}',
        status = 'status/all__run.log',
        cluster = config['cluster']['anno_raw'],
        stdout = 'log/anno_{sample}.o.%j',
        stderr = 'log/anno_{sample}.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        python \
        {params.repo}/snakefile/pyscript/readwise_anno.py \
        --sample {params.sample} \
        --vamos {params.repo}/src/vamos \
        --thread {params.thread} \
        --filter {params.filter}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''


# readwise annotation by efficient motifs (for raw sequencing reads, without alignment)
rule anno_raw_without_aln:
    input:
        phaseDone = 'analysis/{sample}/phase/{sample}.phase_without_aln.done'
    output:
        done = 'analysis/{sample}/anno/{sample}.anno.raw_without_aln.done'
    params:
        repo = config['parameter']['vamos_repo'],
        sample = lambda wildcards: wildcards.sample,
        thread = 16,
        filter = config['parameter']['min_depth'],
        rule = 'anno__{sample}',
        status = 'status/all__run.log',
        cluster = config['cluster']['anno_raw'],
        stdout = 'log/anno_{sample}.o.%j',
        stderr = 'log/anno_{sample}.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        python \
        {params.repo}/snakefile/pyscript/readwise_anno.py \
        --sample {params.sample} \
        --vamos {params.repo}/src/vamos \
        --thread {params.thread} \
        --filter {params.filter}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''


# locuswise annotation (for assembly contigs, with alignment)
rule anno_ass_with_aln:
    input:
        vntr = config['database']['vntr'],
        alnDone = 'analysis/{sample}/aln/{sample}.aln.done',
        bam = lambda wildcards: batchDict[wildcards.sample][1]
    output:
        vcf = 'analysis/{sample}/anno/{sample}.anno.vcf',
        done = 'analysis/{sample}/anno/{sample}.anno.ass_with_aln.done'
    params:
        repo = config['parameter']['vamos_repo'],
        sample = '{sample}',
        rule = 'anno__{sample}',
        status = 'status/all__run.log',
        cluster = config['cluster']['anno_ass'],
        stdout = 'log/anno_{sample}.o.%j',
        stderr = 'log/anno_{sample}.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {params.repo}/src/vamos --locuswise \
        -b {input.bam} \
        -r {input.vntr} \
        -s {params.sample} \
        -o {output.vcf} \
        -t 16

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

# locuswise annotation (for assembly contigs, without alignment)
rule anno_ass_without_aln:
    input:
        vntr = config['database']['vntr'],
        bam = lambda wildcards: batchDict[wildcards.sample][1]
    output:
        vcf = 'analysis/{sample}/anno/{sample}.anno.vcf',
        done = 'analysis/{sample}/anno/{sample}.anno.ass_without_aln.done'
    params:
        repo = config['parameter']['vamos_repo'],
        sample = '{sample}',
        rule = 'anno__{sample}',
        status = 'status/all__run.log',
        cluster = config['cluster']['anno_ass'],
        stdout = 'log/anno_{sample}.o.%j',
        stderr = 'log/anno_{sample}.e.%j'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {params.repo}/src/vamos --locuswise \
        -b {input.bam} \
        -r {input.vntr} \
        -s {params.sample} \
        -o {output.vcf} \
        -t 16

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''


