# *tryvamos*: <ins>t</ins>andem <ins>r</ins>epeat anal<ins>y</ins>sis using *<ins>vamos</ins>* outputs

*tryvamos* is a toolkit designed to perform quick visualization, result transformation, and statistical analysis for tandem repeats using annotations generated from *vamos*.


## Preparation
```
conda create --name tryvamos python=3.10
conda activate tryvamos
conda install --file requirements.txt

python tryvamos.py -h
```

## Usage

### Generate diploid vcf for a panel of samples from haploid/diploid vcfs of individual samples
Combination of vcfs of individual samples is usually the first step after annotation. Most functionality of *tryvamos* uses multi-sample diploid vcf(s) as input. To generate a multi-sample diploid vcf, paths of single-sample vcfs need to be written in a "csv" file, with each line having the two haploid vcfs of one sample (separated by ",") or a single diploid vcf of one sample.
```
python tryvamos.py combineVCF example/example.list example/example.vcf
```
Explanation of "INFO" field of a multi-sample vcf

| Field | Explanation |
|:-------|:---|
| END     | Ending position of the locus |
| RU      | All repeating units (motifs) of the locus ("," separated) |
| SVTYPE  | VNTR |
| ALTANNO | All motif annotations (allels) of the locus ("," separated), motifs are indexed from "0" as ordered in the "RU" field |

Diploid genotype of a sample is recorded as indices of the alleles of each haplotype as ordered in the "ALTANNO" field (start from 1). Unannotated haplotype is denoted as ".".
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HPRC_HG002_h1/HPRC_HG002_h2	HPRC_HG005_h1/HPRC_HG005_h2
chr1	9240	.	N	<VNTR>	.	PASS	END=9259;RU=GATGAGCAGC,ATGAGCAGC,GATGAGCAG;SVTYPE=VNTR;ALTANNO=0-2	GT	1/1	1/1
chr1	52377	.	N	<VNTR>	.	PASS	END=52421;RU=GT;SVTYPE=VNTR;ALTANNO=0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0,0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0	GT	1/1	2/1
```
### Generate feature matrix from combined vcf
Parsing combined vcf of large number of samples for all TR annotations is time consuming. To facilitate repetitive analysis using the combined vcf, we developed the quickFeature module in *tryvamos* to generate sample (row) by feature (column) matrix for convenient high-dimentional data analysis. Currently the following features are supported: "TR length by number of motifs" (annoLen), "TR allele by motif annotations" (annoStr), "counts of the most frequent motif" (topCount), and "reconstructed nucleotide sequence using motif annotations" (nt).
```
python tryvamos.py quickFeature example/example.vcf example/example.topCount.tsv --feature topCount
```
### Visualize landscape of TR alleles by waterfall plot
Waterfall plot is a plot where TR alleles are visualized as concatenation of variable motifs. It is a good tool to visualize the motif and allele landscape of a TR locus in a population. The waterfallPlot plot module in *tryvamos* provide a quick way to generate waterfall plots for all or selected TR loci.
```
python tryvamos.py waterfallPlot example/example.vcf example/plot --useLoci example/loci_for_waterfall.bed
```
### Genomew-wide tests of significance of TRs on two panel of samples

