# -*- coding: utf-8 -*-

import os

"""

"""

if not os.path.exists('analysis'): os.makedirs('analysis')
if not os.path.exists('analysis/boundary'): os.makedirs('analysis/boundary')
if not os.path.exists('analysis/oriMotif'): os.makedirs('analysis/oriMotif')
if not os.path.exists('analysis/effMotif'): os.makedirs('analysis/effMotif')
if not os.path.exists('analysis/effMotif/split'): os.makedirs('analysis/effMotif/split')
if not os.path.exists('analysis/effMotif/ILP'): os.makedirs('analysis/effMotif/ILP')
if not os.path.exists('analysis'): os.makedirs('analysis')
if not os.path.exists('log'): os.makedirs('log')
if not os.path.exists('status'): os.makedirs('status')

splits = int(config['parameter']['splits'])
splitsIDX = range(1,splits+1)

chrs = config['parameter']['chrs'].split(',')
manifestSample, manifestChr, samples = {}, {}, []
with open(config['input']['manifest']) as f:
    for line in f:
        sample,inFA,inBAM,inTRF,inRM = line.strip().split(',')
        samples.append(sample)
        manifestSample[(sample)] = [inFA, inBAM, inTRF, inRM]

for chr in chrs: manifestChr[(chr)] = [chr]

rule all:
    input:
        [f'analysis/boundary/{sample}/{sample}.boundaryDecompose.done' for sample in samples],
        [f'analysis/boundary/{sample}/{sample}.boundaryDecomposeLift.done' for sample in samples],
        [f'analysis/boundary/{sample}/{sample}.boundaryDecomposeLiftRefine.done' for sample in samples],
        [f'analysis/boundary/{chr}/{chr}.boundaryDecomposeLiftRefineAggregate.done' for chr in chrs],
        f'analysis/boundary/boundaryFinal.done',
        [f'analysis/oriMotif/{sample}/{sample}.oriMotifLift.done' for sample in samples],
        [f'analysis/oriMotif/{sample}/{sample}.oriMotifLiftDecompose.done' for sample in samples],
        f'analysis/oriMotif/oriMotifFinal.done',
        f'analysis/effMotif/effMotifSplit.done',
        [f'analysis/effMotif/ILP/effMotif.{splitsIDX}.tsv' for splitsIDX in splitsIDX],
        f'analysis/effMotif/effMotifFinal.done'


############################################################
# major step 1: build TR boundaries from TRF and/or RM calls
############################################################
'''
This rule combines and decomposes the Tandem Repeat Finder and RepeatMasker TR
calls. The python program also supports only TRF or RM input (need to modify the
rule syntax a bit).
'''
rule boundaryDecompose:
    input:
        fasta = lambda wildcards: manifestSample[(wildcards.sample)][0],
        TRF = lambda wildcards: manifestSample[(wildcards.sample)][2],
        RM = lambda wildcards: manifestSample[(wildcards.sample)][3]
    output:
        boundaryDecompose = 'analysis/boundary/{sample}/{sample}.boundaryDecompose.tsv',
        done = 'analysis/boundary/{sample}/{sample}.boundaryDecompose.done'
    params:
        rule = 'boundaryDecompose__{sample}',
        status = 'status/run.log',
        cluster = config['cluster']['boundaryDecompose'],
        stdout = 'log/boundaryDecompose__{sample}.%j.out',
        stderr = 'log/boundaryDecompose__{sample}.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {config[software][python3]} \
        {config[lib]}/{config[tool][boundaryDecompose]} \
        {input.fasta} \
        {output.boundaryDecompose} \
        --inTRF {input.TRF} \
        --inRM {input.RM}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
This rule lifts the boundaryDecompose regions to give the reference regions.
Reference is dependent on the reference used by the input bam files.
'''
rule boundaryDecomposeLift:
    input:
        done = 'analysis/boundary/{sample}/{sample}.boundaryDecompose.done',
        fasta = lambda wildcards: manifestSample[(wildcards.sample)][0],
        bam = lambda wildcards: manifestSample[(wildcards.sample)][1],
        boundaryDecompose = 'analysis/boundary/{sample}/{sample}.boundaryDecompose.tsv'
    output:
        boundaryDecomposeLift = 'analysis/boundary/{sample}/{sample}.boundaryDecomposeLift.tsv',
        done = 'analysis/boundary/{sample}/{sample}.boundaryDecomposeLift.done'
    params:
        rule = 'boundaryDecomposeLift__{sample}',
        status = 'status/run.log',
        cluster = config['cluster']['boundaryDecomposeLift'],
        stdout = 'log/boundaryDecomposeLift__{sample}.%j.out',
        stderr = 'log/boundaryDecomposeLift__{sample}.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {config[software][python3]} \
        {config[lib]}/{config[tool][boundaryDecomposeLift]} \
        {input.boundaryDecompose} \
        {input.fasta} \
        {input.bam} \
        {output.boundaryDecomposeLift}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
The boundaryDecomposeLift regions may have overlaps again on the reference. This
rule refines the boundaryDecomposeLift regions by picking for each overlapping
group the entry whose contig is most frequently aligned to the same chromosome.
Output is split by chromosomes to speed up the rule boundaryDecomposeRefineFinal.
'''
rule boundaryDecomposeLiftRefine:
    input:
        done = 'analysis/boundary/{sample}/{sample}.boundaryDecomposeLift.done',
        boundaryDecomposeLift = 'analysis/boundary/{sample}/{sample}.boundaryDecomposeLift.tsv'
    output:
        done = 'analysis/boundary/{sample}/{sample}.boundaryDecomposeLiftRefine.done'
    params:
        outPre = 'analysis/boundary/{sample}/{sample}.boundaryDecomposeLiftRefine',
        rule = 'boundaryDecomposeLiftRefine__{sample}',
        status = 'status/run.log',
        cluster = config['cluster']['boundaryDecomposeLiftRefine'],
        stdout = 'log/boundaryDecomposeLiftRefine__{sample}.%j.out',
        stderr = 'log/boundaryDecomposeLiftRefine__{sample}.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {config[software][python3]} \
        {config[lib]}/{config[tool][boundaryDecomposeLiftRefine]} \
        {input.boundaryDecomposeLift} \
        {params.outPre} \
        --byChr

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
TR calls from different samples may not align perfectly by the coordinates
and/or consensuses. The program aggregates TR calls from different samples,
grouping calls with any overlaps and clustering by coordinates and consensus
length. Clusters are regrouped by overlaps to give the final TR boundary.
Note the python program assumes all input calls are on the same chromosome.
'''
rule boundaryDecomposeLiftRefineAggregate:
    input:
        done = expand('analysis/boundary/{sample}/{sample}.boundaryDecomposeLiftRefine.done', sample=samples)
    output:
        boundaryDecomposeLiftRefineAggregate = 'analysis/boundary/{chr}/boundaryDecomposeLiftRefineAggregate.{chr}.tsv',
        cluster = 'analysis/boundary/{chr}/boundaryDecomposeLiftRefineAggregate.{chr}.cluster.txt',
        done = 'analysis/boundary/{chr}/{chr}.boundaryDecomposeLiftRefineAggregate.done'
    params:
        chr = lambda wildcards: manifestChr[(wildcards.chr)][0],
        csv = 'analysis/boundary/{chr}/boundaryDecomposeLiftRefineAggregate.{chr}.csv',
        rule = 'boundaryDecomposeLiftRefineAggregate__{chr}',
        status = 'status/run.log',
        cluster = config['cluster']['boundaryDecomposeLiftRefineAggregate'],
        stdout = 'log/boundaryDecomposeLiftRefineAggregate__{chr}.%j.out',
        stderr = 'log/boundaryDecomposeLiftRefineAggregate__{chr}.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        cut -d',' {config[input][manifest]} -f1 > temp.{params.chr}
        cat temp.{params.chr} | while read sample; do
            touch analysis/boundary/${{sample}}/${{sample}}.boundaryDecomposeLiftRefine.{params.chr}.tsv
        done
        ls -d analysis/boundary/*/*boundaryDecomposeLiftRefine.{params.chr}.tsv > temp1.{params.chr}
        paste -d',' temp.{params.chr} temp1.{params.chr} > {params.csv}
        rm temp.{params.chr} temp1.{params.chr}

        {config[software][python3]} \
        {config[lib]}/{config[tool][boundaryDecomposeLiftRefineAggregate]} \
        {params.csv} \
        {output.boundaryDecomposeLiftRefineAggregate} \
        --outCluster {output.cluster}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
This rule combines results of different chromosomes from the rule
boundaryDecomposeLiftRefineAggregate
'''
rule boundaryFinal:
    input:
        done = expand('analysis/boundary/{chr}/{chr}.boundaryDecomposeLiftRefineAggregate.done', chr=chrs),
        boundaryDecomposeLiftRefineAggregate = expand('analysis/boundary/{chr}/boundaryDecomposeLiftRefineAggregate.{chr}.tsv', chr=chrs)
    output:
        boundaryFinal = 'analysis/boundary/boundaryFinal.tsv',
        done = 'analysis/boundary/boundaryFinal.done'
    params:
        rule = 'boundaryFinal',
        status = 'status/run.log',
        cluster = config['cluster']['boundaryFinal'],
        stdout = 'log/boundaryFinal.%j.out',
        stderr = 'log/boundaryFinal.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        cat {input.boundaryDecomposeLiftRefineAggregate} > {output.boundaryFinal}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

############################################################
# major step 2: build oriMotifs from TR boundaries
############################################################
'''
This rule lifts the TR sequences from each sample assembly input using the
TR boundaries built from previous steps.
'''
rule oriMotifLift:
    input:
        done = 'analysis/boundary/boundaryFinal.done',
        boundaryFinal = 'analysis/boundary/boundaryFinal.tsv',
        fasta = lambda wildcards: manifestSample[(wildcards.sample)][0],
        bam = lambda wildcards: manifestSample[(wildcards.sample)][1]
    output:
        oriMotifLift = 'analysis/oriMotif/{sample}/{sample}.oriMotifLift.fasta',
        done = 'analysis/oriMotif/{sample}/{sample}.oriMotifLift.done'
    params:
        rule = 'oriMotifLift__{sample}',
        status = 'status/run.log',
        cluster = config['cluster']['oriMotifLift'],
        stdout = 'log/oriMotifLift__{sample}.%j.out',
        stderr = 'log/oriMotifLift__{sample}.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {config[software][python3]} \
        {config[lib]}/{config[tool][oriMotifLift]} \
        {input.boundaryFinal} \
        {input.fasta} \
        {input.bam} \
        {output.oriMotifLift}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
This rule decomposes the TR sequences lifted from assemblies by the
corresponding TR consensus.
'''
rule oriMotifLiftDecompose:
    input:
        done = 'analysis/oriMotif/{sample}/{sample}.oriMotifLift.done',
        boundaryFinal = 'analysis/boundary/boundaryFinal.tsv',
        oriMotifLift = 'analysis/oriMotif/{sample}/{sample}.oriMotifLift.fasta'
    output:
        oriMotifLiftDecompose = 'analysis/oriMotif/{sample}/{sample}.oriMotifLiftDecompose.tsv',
        done = 'analysis/oriMotif/{sample}/{sample}.oriMotifLiftDecompose.done'
    params:
        rule = 'oriMotifLiftDecompose__{sample}',
        status = 'status/run.log',
        cluster = config['cluster']['oriMotifLiftDecompose'],
        stdout = 'log/oriMotifLiftDecompose__{sample}.%j.out',
        stderr = 'log/oriMotifLiftDecompose__{sample}.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {config[software][python3]} \
        {config[lib]}/{config[tool][oriMotifLiftDecompose]} \
        {input.boundaryFinal} \
        {input.oriMotifLift} \
        {output.oriMotifLiftDecompose}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
This rule aggregates decomposed TRs from multiple samples as the oriMotif set.
'''
rule oriMotifFinal:
    input:
        done = expand('analysis/oriMotif/{sample}/{sample}.oriMotifLiftDecompose.done', sample=samples)
    output:
        done = 'analysis/oriMotif/oriMotifFinal.done'
    params:
        files = 'analysis/oriMotif/oriMotifFinal.list',
        outPre = 'analysis/oriMotif/oriMotifFinal',
        rule = 'oriMotifFinal',
        status = 'status/run.log',
        cluster = config['cluster']['oriMotifFinal'],
        stdout = 'log/oriMotifFinal.%j.out',
        stderr = 'log/oriMotifFinal.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        ls -d analysis/oriMotif/*/*.oriMotifLiftDecompose.tsv > {params.files}

        {config[software][python3]} \
        {config[lib]}/{config[tool][oriMotifFinal]} \
        {params.files} \
        {config[parameter][chrs]} \
        {params.outPre} \
        --inTE {config[lib]}/{config[database][transposable]} \
        --inCentro {config[lib]}/{config[database][centromere]}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

############################################################
# major step 3: build effMotifs from oriMotifs
############################################################
'''
This rule evenly splits VNTR loci into {splits} files.
'''
rule effMotifSplit:
    input:
        done = 'analysis/oriMotif/oriMotifFinal.done'
    output:
        split = expand('analysis/effMotif/split/oriMotif.{splitsIDX}.tsv', splitsIDX=splitsIDX),
        done = 'analysis/effMotif/effMotifSplit.done'
    params:
        oriMotifFinal = 'analysis/oriMotif/oriMotifFinal.all.tsv',
        splits = config['parameter']['splits'],
        prefix = 'analysis/effMotif/split/oriMotif',
        suffix = 'tsv',
        rule = 'effMotifSplit',
        status = 'status/run.log',
        cluster = config['cluster']['effMotifSplit'],
        stdout = 'log/effMotifSplit.%j.out',
        stderr = 'log/effMotifSplit.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {config[software][python3]} \
        {config[lib]}/{config[tool][effMotifSplit]} \
        {params.oriMotifFinal} \
        {params.splits} \
        {params.prefix} \
        {params.suffix} \
        --skip 1 --random

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
This rule applies ILP solver to obtain effMotif for each oriMotif split.
'''
rule effMotifILP:
    input:
        splitDone = 'analysis/effMotif/effMotifSplit.done',
        oriMotif = 'analysis/effMotif/split/oriMotif.{splitsIDX}.tsv',
    output:
        effMotif = 'analysis/effMotif/ILP/effMotif.{splitsIDX}.tsv'
    params:
        threads = config['parameter']['threadsILP'],
        timeLimit = config['parameter']['timeLimitILP'],
        q = config['parameter']['q'],
        rule = 'effMotifILP__{splitsIDX}',
        status = 'status/effMotifILP.log',
        cluster = config['cluster']['effMotifILP'],
        stdout = 'log/effMotifILP__{splitsIDX}.%j.out',
        stderr = 'log/effMotifILP__{splitsIDX}.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        {config[software][python3]} \
        {config[lib]}/{config[tool][effMotifILP]} \
        {input.oriMotif} \
        {output.effMotif} \
        -t {params.threads} \
        -l {params.timeLimit} \
        -q {params.q}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''

'''
This function aggregates splited effMotif files into one.
'''
rule effMotifFinal:
    input:
        files = expand("analysis/effMotif/ILP/effMotif.{splitsIDX}.tsv", splitsIDX=splitsIDX),
        header = config['lib']+'/'+config['database']['header']
    output:
        effMotifFinal = 'analysis/effMotif/effMotifFinal.all.tsv',
        done = 'analysis/effMotif/effMotifFinal.done'
    params:
        rule = 'effMotifFinal',
        status = 'status/run.log',
        cluster = config['cluster']['effMotifFinal'],
        stdout = 'log/effMotifFinal.%j.out',
        stderr = 'log/effMotifFinal.%j.err'
    shell:
        '''
        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tstart" >> {params.status}
        export LANG=en_US.UTF-8

        cat {input.header} {input.files} > {output.effMotifFinal}

        touch {output.done}

        echo -e "`date +"%Y-%m-%d %H:%M:%S"`\t{params.rule}\tfinish" >> {params.status}
        '''
