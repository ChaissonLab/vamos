#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <getopt.h>
#include "input.h"
#include "vntr.h"
using namespace std;

/* vamos -in read.bam -vntr vntrs.bed -motif motifs.csv -o out.vcf */

void printUsage() 
{
    printf("Usage: vamos [-h] [-i in.bam] [-v vntrs.bed] [-m motifs.csv] [-o output.vcf] [-s sample]\n");
    printf("Options:\n");
    printf("       -i  FILE      input alignment file (bam format)\n");
    printf("       -v  FILE      the Tab-delimited coordinate of each VNTR locus - chrom\tstart\tend, each row represents a VNTR locus\n");
    printf("       -m  FILE      the Comma-delimited motif sequences list for each VNTR locus, each row represents a VNTR locus\n");
    printf("       -o  FILE      output vcf file\n");
    printf("       -s  CHAR      the sample name\n")
    printf("       -h            print out help message\n");
} 


int main (int argc, char **argv)
{
  int c;
  string input_bam, vntr_bed, motif_csv, out_vcf;

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
      input_bam = optarg;
      break;

    case 'v':
      printf ("option -vntr with `%s'\n", optarg);
      vntr_bed = optarg;
      break;

    case 'm':
      printf ("option -motif with `%s'\n", optarg);
      motif_csv = optarg;
      break;

    case 'o':
      printf ("option -output with `%s'\n", optarg);
      out_vcf = optarg;
      break;

    case 's':
      printf ("option -sampleName with `%s'\n", optarg);
      break;

    case 'h':
      printUsage() 
      break;

    case '?':
      printf("Unknown option: %c\n", optopt);
      return 1;

    case ':':
      printf(stderr, "[ERROR] missing option argument\n");
      return 1;

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

  vector<VNTR> vntrs;

  /* read VNTR bed file */
  readVNTRFromBed (vntr_bed, vntrs);

  /* read motif csv file */
  readMotifsFromCsv (motif_csv, vntrs);

  /* process each VNTR */
  for (auto &it: vntrs) 
  {
    it.readSeqFromBam (input_bam);
    it.motifRepresentationForOneSeq(); 
    it.concensusMotifRepresentation();
  }

  exit(EXIT_SUCCESS);
}

