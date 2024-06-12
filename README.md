
# Vamos: VNTR Annotation tool using efficient MOtif Sets

Vamos is a tool to perform motif annotation of VNTR sequences from long-read assemblies
or mapped long-read BAMs. Variant calls are produced from BAMs by prephasing sequences
if no phase information is available, a max-cut heuristic is applied to phase reads. Annotation
is then called on consensus sequences constructed from partial-order alignment of reads
in each haplotype (or all reads if in an autozygous region).


The motifs used to annotate TR sequences are selected from TR annotations of long-read
assemblies.An optimization routine is used to select a subset of motifs from all observed
motifs at each locus so that the sequences used for annotation likely reflect true
motif variation and not sequencing error. 

Vamos guarantees that the encoding sequence is winthin a bounded edit distance of
the original sequence for the genomes used to compile the motif database.

For example, a VNTR sequence ACGGT|ACTGT|ACGT may be encoded to a more compact representation: ACGGT|AC**G**GT|ACGT using efficient motif set [ACGGT, ACGT].
The edit distance between the original VNTR sequence and encoding sequence is 1. 


## Updates
We have released a version 2.1 of the motif set. This has roughly 1.2m loci and uses
an additional motif-harmonization step before efficient motif selection. We additionally
supply loci for CHM13.

## Getting Started

To install vamos, g++ (>= 8.3.0), htslib, abpoa, edlib and alglib are required.
htslib and abpoa can be installed through bioconda.
Static libraries libalglib.a and libedlib.a are distributed along with vamos.

Install required libraries through conda
```
conda create --name vamos python=3.10
conda activate vamos
conda install -c bioconda --file requirements.txt
```

Or download the latest code from github
```
git clone https://github.com/ChaissonLab/vamos.git 
cd vamos*/src/ && make
```

Next you should download a locus list. These are in BED format, with the coordinates of the tandem repeat as the BED coordinates, and the list of observed motifs from the Human Pangenome Reference Consortium given as the extra field.

https://zenodo.org/records/11625069

The latest motif sets as of v2.1 can be downloaded via:

For annotation on GRCh38:
 curl "https://zenodo.org/records/11625069/files/vamos.motif.hg38.v2.1.e0.1.tsv.gz?download=1" > vamos.motif.hg38.v2.1.e0.1.tsv.gz; gunzip vamos.motif.hg38.v2.1.e0.1.tsv.gz
For annotation on CHM13
 curl "https://zenodo.org/records/11625069/files/vamos.motif.CHM13.v2.1.orig.tsv.gz?download=1" > vamos.motif.CHM13.v2.1.orig.tsv.gz; gunzip vamos.motif.CHM13.v2.1.orig.tsv.gz



For running vamos on a haplotype-resolved assembly:
```
vamos --contig -b assembly.hap1.mapped_to_grch38.bam -r vamos.motif.hg38.v2.1.e0.1.tsv -s sample_name -o assembly.hap1.vcf -t 8
vamos --contig -b assembly.hap2.mapped_to_grch38.bam -r vamos.motif.hg38.v2.1.e0.1.tsv -s sample_name -o assembly.hap2.vcf -t 8
```
For running vamos on aligned reads, you can use the following command:
```
vamos --read -b ../example/demo.aln.bam -r vamos.motif.hg38.v2.1.e0.1.tsv -s NA24385_CCS_h1 -o reads.vcf -t 8
```
If the reads are pre-phased using HapCut or WhatsHap, and contain the HA SAM tag, this phasing
will be used to call variants from each haplotype. If the reads are unphased, a max-cut heuristic will
be used to prephase reads before calling variants.



## Table of Contents

- [Introduction](#introduction)
  - [Efficient motif set](#emotif)
- [Installation](#install)
  - [Building vamos from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [General usage](#usage)
  - [To generate annotation for aligned reads at each VNTR locus](#read)
  - [To generate annotation for haplotype-resolved assembly at each VNTR locus](#assembly)
- [Commands and options](#cmd)
- [Output](#output)
  - [VCF](#VCF)
- [Combine VCFs](#combine)

## <a name="introduction"></a>Introduction

### <a name="emotif"></a>Efficient motif set
We defined VNTR loci and motifs using a collection of 32 haplotype-resolved LRS genomes constructed by Human Genome Structural Variation Consortium.
692,882 loci of simple repeating sequences on the GRCh38 assembly were obtained from the table browser tool of the UCSC Genome Browser.
For each assembly, VNTR sequences were lifted-over and decomposed into motifs by Tandem Repeats Finder (TRF). Post-filtering step leaves 467104 well-resolved VNTR loci. 
We propose efficient motif set as a smallest set of motifs, such that the string decompositions of the assembly alleles are bounded by a given edit distance. 
[snakefile/configs/vntr_region_motifs.e.bed.gz](https://github.com/ChaissonLab/vamos/blob/master/snakefile/configs/vntr_region_motifs.e.bed.gz) provides efficent motifs for 467104 VNTR loci.
[snakefile/configs/vntr_region_motifs.o.bed.gz](https://github.com/ChaissonLab/vamos/blob/master/snakefile/configs/vntr_region_motifs.o.bed.gz) provides original motifs for 467104 VNTR loci.

## <a name="install"></a>Installation
### <a name="build"></a>Building vamos from source files
You can build vamos from source files. 
Make sure to install the required libraries: g++ (>= 8.3.0), htslib, abpoa, edlib and alglib before compiling. 

Install required libraries
```
conda create --name vamos python=3.10
conda activate vamos
conda install -c bioconda --file requirements.txt
```

Download the [latest release](https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.2.0.tar.gz)
```
wget https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.2.0.tar.gz
tar -zxvf vamos-v1.2.0.tar.gz
cd vamos-v1.2.0/src; make
```
Or, you can use `git clone` command to download the source code.
This gives you the latest version of vamos, which might be still under development.
```
git clone https://github.com/ChaissonLab/vamos.git 
cd vamos/src; make
```

### <a name="binary"></a>Pre-built binary executable file for Linux/Unix 
If you meet any compiling issue, please try the pre-built binary file:
```
wget https://github.com/ChaissonLab/vamos/releases/download/vamos-v1.2.0/vamos-v1.2.0_x64-linux
tar -zxvf vamos-v1.2.0_x64-linux
```

## <a name="usage"></a>General usage
### <a name="read"></a>To generate annotation for aligned reads at each VNTR locus
```
vamos --read -b ../example/demo.aln.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o reads.vcf -t 8
```
### <a name="assembly"></a>To generate annotation for haplotype-resolved assembly
```
vamos --contig -b assembly.hap1.mapped_to_grch38.bam -r emotifs.d10.64h.bed -s sample_name -o assembly.hap1.vcf -t 8
vamos --contig -b assembly.hap2.mapped_to_grch38.bam -r emotifs.d10.64h.bed -s sample_name -o assembly.hap2.vcf -t 8
```


## <a name="cmd"></a>Commands and options
```
Usage: vamos [subcommand] [options] [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads]
Version: v1.2.0
subcommand:
vamos --contig [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads]
vamos --read [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] [-p phase_flank]
vamos -m [verison of efficient motif set]

   Input:
       -b   FILE         Input indexed bam file.
       -r   FILE         File containing region coordinate and motifs of each VNTR locus.
                         The file format: columns `chrom,start,end,motifs` are tab-delimited.
                         Column `motifs` is a comma-separated (no spaces) list of motifs for this VNTR.
       -s   CHAR         Sample name.
   Output:
       -o   FILE         Output vcf file.
   Dynamic Programming:
       -d   DOUBLE       Penalty of indel in dynamic programming (double) DEFAULT: 1.0.
       -c   DOUBLE       Penalty of mismatch in dynamic programming (double) DEFAULT: 1.0.
       -a   DOUBLE       Global accuracy of the reads. DEFAULT: 0.98.
       --naive           Specify the naive version of code to do the annotation, DEFAULT: faster implementation.
   Phase reads:
       -p   INT          Range of flanking sequences which is used in the phasing step. DEFAULT: 3000 bps.
   Downloading motifs:
       -m  MOTIF         Prints a command to download a particular motif set. Current supported motif set is: q20.
                         This motif set is selected at a level of Delta=20 from 64 haplotype-resolvd assemblies (Ebert et al., 2021)
                         This may be copied and pasted in the command line, or executed as: vamos -m q20
   Others:
       -t   INT          Number of threads, DEFAULT: 1.
       --debug           Print out debug information.
       -h                Print out help message.
```

### <a name="VCF"></a>VCF for locuswise and single_seq modes
Vamos generates annotation for each VNTR locus in VCF file under `--contig` and `--read` modes.


The following shows an example of the VCF file
```
##fileformat=VCFv4.1
##source=vamos_v1.1.0
##contig=<ID=chr1,length=248956422>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=RU,Number=1,Type=String,Description="Comma separated motif sequences list in the reference orientation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=ALTANNO_H1,Number=1,Type=String,Description=""Motif representation for the h1 alternate allele>
##INFO=<ID=ALTANNO_H2,Number=1,Type=String,Description=""Motif representation for the h2 alternate allele>
##INFO=<ID=LEN_H1,Number=1,Type=Integer,Description=""Length of the motif annotation for the h1 alternate allele>
##INFO=<ID=LEN_H2,Number=1,Type=Integer,Description=""Length of the motif annotation for the h2 alternate allele>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=VNTR,Description="Allele comprised of VNTR repeat units">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA24385_CCS_h1
chr1	191351	.	N	<VNTR>	.	PASS	END=191386;RU=CACCACAGAAAACAGAG,CACCACAGAAAACAGAGC;SVTYPE=VNTR;ALTANNO_H1=1,0;LEN_H1=2;	GT	1/1
```

## <a name="combine"></a>Generate multi-sample vcf
A *python* script [snakefile/pyscript/combine_vcf.py](https://github.com/ChaissonLab/vamos/blob/master/snakefile/pyscript/combine_vcf.py) is developed that combines vcfs output by *vamos* into multi-sample vcf.
```
python combine_vcf.py -i vcf.list -o combine.vcf
```
The input `vcf.list` is simply a list of all vcfs to be combined, one line for each vcf. The combined vcf records all unique alleles of the input samples in the `info` field and records the genotype of each sample accordingly.
