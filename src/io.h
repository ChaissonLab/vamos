#ifndef IO_H_
#define IO_H_
#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <string>
#include <tuple> 
#include <vector>
#include "read.h"
#include "vntr.h"
#include "vcf.h"

using namespace std;

class IO
{
public:
    char * region_and_motifs;
    char * input_fasta;
	char * input_bam;
	char * vntr_bed;
	char * motif_csv;
	char * out_vcf;
	char * sampleName;
	char * version;
	OutWriter outWriter;
  // Not the best place to put this, but since the IO is batched and we
  // do not want to store the flanking sequences at each locus for
  // the duration of the run of the program, it's passed through
  // to here for now. Eventually if sites are processed immediately
  // after overlapping reads are read, this can move.
  int phaseFlank;
	IO () 
	{
		version = (char *) malloc(7);
		strcpy(version, "1.2.0");
		region_and_motifs = NULL;
		input_bam = NULL;
		vntr_bed = NULL;
		motif_csv = NULL;
		out_vcf = NULL;
		sampleName = NULL;
		phaseFlank = 0;
	};

	~IO() 
	{
        free(region_and_motifs);
        free(version);
        free(input_bam);
		free(vntr_bed);
		free(motif_csv);
		free(out_vcf);
		free(sampleName);
	};

    int readRegionAndMotifs (vector<VNTR*> &vntrs);
  
	int readMotifsFromCsv (vector<VNTR *> &vntrs);

	int read_tsv(vector<vector<string>> &items);

	void readVNTRFromBed (vector<VNTR *> &vntrs);

	/* get the sequences from input_bam_file that overlapping with chr:start-end */
	void readSeqFromBam (vector<VNTR *> &vntrs, int nproc, int cur_thread, int sz);

	// void readSeqFromBam (vector<READ*> &reads, string &chr, const uint32_t &ref_VNTR_start, 
 //                       const uint32_t &ref_VNTR_end, const uint32_t &VNTR_len, string &region);
	
	// int outputVCF (vector<VNTR *> &vntrs);

	int writeVCFHeader_locuswise(ofstream& out);

	int writeVCFBody_locuswise(ofstream& out, vector<VNTR *> &vntrs, int tid, int nproc);

	int writeBEDHeader_readwise(ofstream& out);

	int writeBEDBody_readwise(ofstream& out, vector<VNTR *> &vntrs, int tid, int nproc);

	void writeFa(ofstream& out, vector<VNTR *> &vntrs);

	void readSeqFromFasta(vector<VNTR *> &vntrs);
};

#endif
