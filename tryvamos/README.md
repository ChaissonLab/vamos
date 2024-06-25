# *tryvamos*: <ins>t</ins>andem <ins>r</ins>epeat anal<ins>y</ins>sis using *<ins>vamos</ins>* outputs

*tryvamos* is a tool package designed to perform quick visualization, result transformation, and statistical analysis for tandem repeats using annotations generated from *vamos*.


## Preparation
```
conda create --name tryvamos python=3.10
conda activate tryvamos
conda install --file requirements.txt

python tryvamos.py -h
```

## Usage

### Generate diploid vcf for a panel of samples from haploid/diploid vcfs of indivudual samples
Combination of vcfs of individual samples is usually the first step after annotation. Most functionality of *tryvamos* uses a multi-sample diploid vcf as input. To generate a multi-sample diploid vcf, paths of sample-wise vcfs need to be written in a "csv" file, with each line having the two haploid vcfs of one sample (separated by ",") or a single diploid vcf of one sample.
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
## Generate feature matrix from combined vcf
