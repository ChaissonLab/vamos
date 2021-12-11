#include <stdio.h>     
#include <stdlib.h>   
#include <getopt.h>
#include <string.h>
#include "io.cpp"
#include "vntr.h"
#include "vcf.h"
using namespace std;

/* vamos -in read.bam -vntr vntrs.bed -motif motifs.csv -o out.vcf */

void printUsage() 
{
	printf("Usage: vamos [-h] [-i in.bam] [-v vntrs.bed] [-m motifs.csv] [-o output.vcf] [-s sample_name]\n");
	printf("Options:\n");
	printf("       -i  FILE      input alignment file (bam format)\n");
	printf("       -v  FILE      the Tab-delimited coordinate of each VNTR locus - chrom\tstart\tend, each row represents a VNTR locus\n");
	printf("       -m  FILE      the Comma-delimited motif sequences list for each VNTR locus, each row represents a VNTR locus\n");
	printf("       -o  FILE      output vcf file\n");
	printf("       -s  CHAR      the sample name\n");
	printf("       -h            print out help message\n");
} 


int main (int argc, char **argv)
{
	int c;
	char * input_bam;
	char * vntr_bed;
	char * motif_csv;
	char * out_vcf;
	char * sampleName;

	while (1)
	{
		static struct option long_options[] =
		{
			/* These options don’t set a flag. We distinguish them by their indices. */
			{"input",         no_argument,       0, 'i'},
			{"vntr",          no_argument,       0, 'v'},
			{"motif",         required_argument, 0, 'm'},
			{"output",        required_argument, 0, 'o'},
			{"sampleName",    required_argument, 0, 's'},
			{"help",          no_argument,       0, 'h'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		int sz;

		c = getopt_long (argc, argv, "i:b:m:o:s:h", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1)
			break;

		switch (c)
		{
		case 0:
			/* If this option set a flag, do nothing else now. */
			if (long_options[option_index].flag != 0)
				break;
			printf ("option %s", long_options[option_index].name);
			if (optarg)
				printf (" with arg %s", optarg);
			printf ("\n");
			break;

		case 'i':
			printf ("option -input with `%s'\n", optarg);
			sz = strlen(optarg);
			input_bam = (char *) malloc(sz);
			strcpy(input_bam, optarg);
			break;

		case 'v':
			printf ("option -vntr with `%s'\n", optarg);
			sz = strlen(optarg);
			vntr_bed = (char *) malloc(sz);
			strcpy(vntr_bed, optarg);
			break;

		case 'm':
			printf ("option -motif with `%s'\n", optarg);
			sz = strlen(optarg);
			motif_csv = (char *) malloc(sz);
			strcpy(motif_csv, optarg);
			break;

		case 'o':
			printf ("option -output with `%s'\n", optarg);
			sz = strlen(optarg);
			out_vcf = (char *) malloc(sz);
			strcpy(out_vcf, optarg);
			break;

		case 's':
			printf ("option -sampleName with `%s'\n", optarg);
			sz = strlen(optarg);
			sampleName = (char *) malloc(sz);
			strcpy(sampleName, optarg);
			break;

		case 'h':
			printUsage();
			exit(EXIT_SUCCESS);

		case '?':
			printf("Unknown option: %c\n", optopt);
			exit(EXIT_FAILURE);

		case ':':
			cerr << "[ERROR] missing option argument" << endl;
			exit(EXIT_FAILURE);

		default:
			abort ();
		}
	}
	if (optind < argc) {
		printf("non-option ARGV-elements: ");
		while (optind < argc)
			printf("%s ", argv[optind++]);
		printf("\n");
	}

	char version[] = "V1.0.0";
	vector<VNTR> vntrs;
	VcfWriter vcfWriter(input_bam, version, sampleName);

	/* read VNTR bed file */
	readVNTRFromBed(vntr_bed, vntrs);

	/* read motif csv file */
	readMotifsFromCsv(motif_csv, vntrs);

	/* process each VNTR */
	for (auto &it: vntrs) 
	{
		readSeq(input_bam, it);
		it.motifRepresentationForOneVNTR(); 
		it.concensusMotifRepForOneVNTR();
	}

	ofstream out(out_vcf);
    if (out.fail()) 
    {
        cerr << "Unable to open file " << out_vcf << endl;
        exit(EXIT_FAILURE);
    }
	VcfWriteHeader(out, vcfWriter);
	VCFWriteBody(vntrs, vcfWriter, out);
	out.close();

	exit(EXIT_SUCCESS);
}

