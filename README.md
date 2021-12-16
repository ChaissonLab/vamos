## Getting Started
```sh
# Install vamos (g++ and htslib seqan edlib alglib required)
# htslib and seqan are installed by bioconda
# libalglib.a and libedlib.a are distributed along with vamos
git clone https://github.com/ChaissonLab/vamos.git
conda create --name vamos
conda activate vamos
conda install -c bioconda htslib  
conda install -c bioconda seqan 
cd vamos/src && make

# help page
vamos -h

# run
vamos -i in.bam -v vntrs.bed -m motifs.csv -o out.vcf -s sampleName
```

## Output files


## Special Notes


## Known issues


## Preliminary results
