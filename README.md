# Vamos: VNTR Annotation tool using MOtifs Selection
## Getting Started
```sh
# Install vamos (g++, htslib, abpoa, edlib and alglib required)
# htslib and abpoa are installed by bioconda 
# libalglib.a and libedlib.a are distributed along with vamos)
git clone https://github.com/ChaissonLab/vamos.git
conda create --name vamos
conda activate vamos
conda install -c bioconda htslib  
conda install -c bioconda abpoa 
cd vamos/src && make

# help page
vamos -h

# run
vamos -i in.bam -v vntrs.bed -m motifs.csv -o out.vcf -s sampleName -f
```

## Output VCF files
Vamos generates a separate VCF record for each VNTR with information about VNTR's location and genotype. The records for alternate VNTR alleles are demarcated by <VNTR> symbolically.

The header of the VCF file contains a detailed description of each record.

### Example

The following VCF entry describes the state of a VNTR locus in a sample with
name sfs

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sfs
1	18000	.	N	<VNTR>	.	PASS	END=19000;RU=AAATTTTTTGGGCCC,ATTTTGGGCCCCC,AAAAAACCCCCCT;SVTYPE=VNTR;ALTANNO_H1=MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_0,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1;ALTANNO_H2=MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_0,MOTIF_0,MOTIF_0,MOTIF_1,MOTIF_0,MOTIF_1,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_2,MOTIF_0,MOTIF_0,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_0,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_0,MOTIF_1,MOTIF_1,MOTIF_0,MOTIF_0,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_2,MOTIF_2,MOTIF_2,MOTIF_1,MOTIF_1,MOTIF_1,MOTIF_12,MOTIF_12,MOTIF_12,MOTIF_12,MOTIF_12,MOTIF_12,MOTIF_12,MOTIF_12,MOTIF_12;	PASS	GT	1/2
```

This line tells us that 
## Special Notes


## Known issues


## Preliminary results
