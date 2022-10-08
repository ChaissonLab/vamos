
# Vamos: VNTR Annotation tool using efficient MOtifs Sets

## Updates 

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

Download the [latest release](https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.1.0.tar.gz)
```
wget https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.1.0.tar.gz
tar -zxvf vamos-v1.1.0.tar.gz && cd vamos-v1.1.0
```

Or download the latest code from github
```
git clone https://github.com/ChaissonLab/vamos.git 
```

Make from source and run with test data:
```
cd vamos*/src/ && make
```
For running vamos on a haplotype-resolved assembly:
```
vamos --contig -b assembly.hap1.mapped_to_grch38.bam -r emotifs.d10.64h.bed -s sample_name -o assembly.hap1.vcf -t 8
vamos --contig -b assembly.hap2.mapped_to_grch38.bam -r emotifs.d10.64h.bed -s sample_name -o assembly.hap2.vcf -t 8
```
For running vamos on aligned reads (phased or unphased):
```
vamos --read -b ../example/demo.aln.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o reads.vcf -t 8
```
Note, if the reads are pre-phased (e.g. by HapCut or WhatsHap) and have the HA SAM tag, the phasing heuristic will not be applied. 


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
Vamos is a tool to perform run-length encoding of VNTR sequences using a set of selected motifs from all motifs observed at that locus.
Vamos guarantees that the encoding sequence is winthin a bounded edit distance of the original sequence. 
For example, a VNTR sequence ACGGT|ACTGT|ACGT may be encoded to a more compact representation: ACGGT|AC**G**GT|ACGT using efficient motif set [ACGGT, ACGT].
The edit distance between the original VNTR sequence and encoding sequence is 1. 

Vamos can generate annotation for haplotype-resolved assembly at each VNTR locus, given a set of motifs at that VNTR locus. 
Vamos can generate annotation for aligned reads (phased or unphased) at each VNTR locus. 

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

Download the [latest release](https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.1.0.tar.gz)
```
wget https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.1.0.tar.gz
tar -zxvf vamos-v1.1.0.tar.gz
cd vamos-v1.1.0/src; make
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
wget https://github.com/ChaissonLab/vamos/releases/download/vamos-v1.1.0/vamos-v1.1.0_x64-linux.tar.gz
tar -zxvf vamos-v1.1.0_x64-linux.tar.gz
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
Version: v1.1.0
subcommand:
vamos --contig [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] 
vamos --read [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] [-p phase_flank] 
   Input: 
       -b   FILE         input indexed bam file (when using --readwise and --locuswise) or fasta file (when using --single_seq). 
       -r   FILE         file containing region coordinate and motifs of each VNTR locus. 
                         The file format: columns `chrom,start,end,motifs` are tab-delimited. 
                         Column `motifs` is a comma-separated (no spaces) list of motifs for this VNTR. 
       -s   CHAR         sample name. 
   Output: 
       -o   FILE         output bed (when using --readwise) or vcf (when using --locuswise and --single_seq) file. 
   Dynamic Programming: 
       -d   DOUBLE       penalty of indel in dynamic programming (double) DEFAULT: 1.0. 
       -c   DOUBLE       penalty of mismatch in dynamic programming (double) DEFAULT: 1.0. 
       -a   DOUBLE       Global accuracy of the reads. DEFAULT: 0.98. 
       --naive           specify the naive version of code to do the annotation, DEFAULT: faster implementation. 
   Aggregate Annotation: 
       -f   DOUBLE       filter out noisy read annotations, DEFAULT: 0.0 (no filter). 
   Phase reads: 
       -p   INT          the range of flanking sequences which is used in the phasing step. DEFAULT: 3000 bps. 
   Downloading motifs:
       -m  MOTIF.        Prints a command to download a particular motif set. Current supported motif sets are: d10e32
                         Delta=10 generated from 32 haplotype-resolvd assemblies (Ebert et al., 2021)
                         This may be copied and pasted in the command line, or executed as: vamos -m d10e32
   Others: 
       -t   INT          number of threads, DEFAULT: 1. 
       --debug           print out debug information. 
       -h                print out help message. 
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
