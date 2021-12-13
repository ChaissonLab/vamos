#include <stdio.h>     
#include <stdlib.h>   
#include <getopt.h>
#include <string>
#include <unistd.h>
#include "io.h"
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
	IO io;

	while (1)
	{
		static struct option long_options[] =
		{
			/* These options don’t set a flag. We distinguish them by their indices. */
			{"input",         required_argument,       0, 'i'},
			{"vntr",          required_argument,       0, 'v'},
			{"motif",         required_argument,       0, 'm'},
			{"output",        required_argument,       0, 'o'},
			{"sampleName",    required_argument,       0, 's'},
			{"help",          no_argument,             0, 'h'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		c = getopt_long (argc, argv, "i:v:m:o:s:h", long_options, &option_index);

		/* Detect the end of the options. */
		if (c == -1) {
			printUsage();
			exit(EXIT_SUCCESS);
		}

		switch (c)
		{
		case 'i':
			printf ("option -input with `%s'\n", optarg);
			io.input_bam = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.input_bam, optarg);
			break;

		case 'v':
			printf ("option -vntr with `%s'\n", optarg);
			io.vntr_bed = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.vntr_bed, optarg);
			break;

		case 'm':
			printf ("option -motif with `%s'\n", optarg);
			io.motif_csv = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.motif_csv, optarg);
			break;

		case 'o':
			printf ("option -output with `%s'\n", optarg);
			io.out_vcf = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.out_vcf, optarg);
			break;

		case 's':
			printf ("option -sampleName with `%s'\n", optarg);
			io.sampleName = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.sampleName, optarg);
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
			printUsage();
			exit(EXIT_SUCCESS);
			// abort ();
		}
	}
	if (optind < argc) {
		printf("non-option ARGV-elements: ");
		while (optind < argc)
			printf("%s ", argv[optind++]);
		printf("\n");
	}

	vector<VNTR> vntrs;

	/* read VNTR bed file */
	io.readVNTRFromBed(vntrs);

	/* read motif csv file */
	io.readMotifsFromCsv(vntrs);

	/* process each VNTR */
	for (auto &it: vntrs) 
	{
		io.readSeq(it);
		it.motifRepresentationForOneVNTR(); 
		it.concensusMotifRepForOneVNTR();
	}

	io.outputVCF(vntrs);
	exit(EXIT_SUCCESS);
}

