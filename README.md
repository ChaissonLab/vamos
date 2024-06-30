

# *vamos*: <ins>V</ins>NTR <ins>a</ins>nnotation tool using efficient <ins>mo</ins>tif <ins>s</ins>ets

*vamos* is a tool to perform motif annotation of VNTR sequences from long-read assemblies or mapped long-read BAMs. Variant calls are produced from BAMs by prephasing sequences if no phase information is available, a max-cut heuristic is applied to phase reads. Annotation is then called on consensus sequences constructed from partial-order alignment of reads in each haplotype (or all reads if in an autozygous region).

The motifs used to annotate TR sequences are selected from TR annotations of long-read assemblies. An optimization routine is used to select a subset of motifs from all observed motifs at each locus so that the sequences used for annotation likely reflect true motif variation and not sequencing error. 

Vamos guarantees that the encoding sequence is winthin a bounded edit distance of the original sequence for the genomes used to compile the motif database.

For example, a VNTR sequence ACGGT|ACTGT|ACGT may be encoded to a more compact representation: ACGGT|AC**G**GT|ACGT using efficient motif set [ACGGT, ACGT]. The edit distance between the original VNTR sequence and encoding sequence is 1. 


## Latest Updates
We have released a version 2.1 of the motif set. This has roughly 1.2m loci and uses an additional motif-harmonization step before efficient motif selection. We additionally supply loci for CHM13.

## Getting Started
To install *vamos*, g++ (>= 8.3.0), htslib, abpoa, edlib and alglib are required.
htslib and abpoa can be installed through bioconda.
Static libraries libalglib.a and libedlib.a are distributed along with *vamos*.

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

The latest motif sets as of v2.1 can be downloaded via: https://zenodo.org/records/11625069

For annotation on GRCh38 using the *vamos* **efficient** motifs:
```
 curl "https://zenodo.org/records/11625069/files/vamos.motif.hg38.v2.1.e0.1.tsv.gz?download=1" > vamos.motif.hg38.v2.1.e0.1.tsv.gz; gunzip vamos.motif.hg38.v2.1.e0.1.tsv.gz
 ```
For annotation on CHM13 using the *vamos* **efficient** motifs
```
 curl "https://zenodo.org/records/11625069/files/vamos.motif.CHM13.v2.1.e0.1.tsv.gz?download=1" > vamos.motif.CHM13.v2.1.e0.1.tsv.gz; gunzip vamos.motif.CHM13.v2.1.e0.1.tsv.gz
```

## Running *vamos*
For running *vamos* on a haplotype-resolved assembly:
```
vamos --contig -b assembly.hap1.mapped_to_grch38.bam -r vamos.motif.hg38.v2.1.e0.1.tsv -s sample_name -o assembly.hap1.vcf -t 8
vamos --contig -b assembly.hap2.mapped_to_grch38.bam -r vamos.motif.hg38.v2.1.e0.1.tsv -s sample_name -o assembly.hap2.vcf -t 8
```
For running *vamos* on aligned reads:
```
vamos --read -b ../example/demo.aln.bam -r vamos.motif.hg38.v2.1.e0.1.tsv -s NA24385_CCS_h1 -o reads.vcf -t 8
```
If the reads are pre-phased using HapCut or WhatsHap, and contain the HA SAM tag, this phasing will be used to call variants from each haplotype. If the reads are unphased, a max-cut heuristic will be used to prephase reads before calling variants.

## Output
*vamos* generates single-sample diploid vcf (by the ```--read``` mode) or haploid vcf (by the ```--contig``` mode).

Explanation of the "INFO" field of an output vcf:
| Field | Explanation |
|:-------|:---|
| END     | Ending position of the locus |
| RU      | All repeating units (motifs) of the locus ("," separated) |
| SVTYPE  | VNTR |
| ALTANNO_H1 | Motif annotations (allels) of haplotyp1, motifs are indexed from "0" as ordered in the "RU" field |
| LEN_H1 | Total count of motifs of haplotyp1 allele |
| ALTANNO_H2 | Motif annotations (allels) of haplotyp2, motifs are indexed from "0" as ordered in the "RU" field |
| LEN_H2 | Total count of motifs of haplotyp2 allele |


The following shows example entries of a diploid single-sample vcf file
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	11226	.	N	<VNTR>	.	PASS	END=11462;RU=AGAGTGGTGGCCAGCGCCCCCTGCTGGCGCCGGGGCACTGCAGGGCCCTCTTGCTT,ACTGTATAGTGGTGGCACGCCGCCTGCTGGCAGCTAGGGACATTGCAGGGTCCTCTTGCTC,AGAGTGGTGGCCACCGCCCCCTGCTGGCGCCGGGGCACTGCAGGGTCCTCTTGCTT,ACTGTATAGTGGTGGCACGCCGCCTGCTGGCAGCTACGGACATTGCAGGGTCCTCTTGCTC,ACTGTATAGTGGTGGCACGCCGCCTGCTGGCAGCTAGGGACATTGCAGGGTCCTCTTGCTCA;SVTYPE=VNTR;ALTANNO_H1=0,0,4,0;LEN_H1=4;	GT	1/1
chr1	15796	.	N	<VNTR>	.	PASS	END=15849;RU=CTT,CTC,CTG,CAG,CATG;SVTYPE=VNTR;ALTANNO_H1=0,1,2,1,0,2,0,0,0,1,3,0,2,1,0,4,2;LEN_H1=17;ALTANNO_H2=0,1,2,1,0,2,2,2,0,1,3,0,2,1,0,4,2;LEN_H2=17;	GT	1/2
```

## Analysis using *vamos* output
Please refer to [*tryvamos*](https://github.com/ChaissonLab/vamos/blob/master/tryvamos) for more details.
