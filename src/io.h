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
	char * input_bam;
	char * vntr_bed;
	char * motif_csv;
	char * out_vcf;
	char * sampleName;
	char * version;
	VcfWriter vcfWriter;

	IO () 
	{
		version = (char *) malloc(7);
		strcpy(version, "V1.0.0");
		input_bam = NULL;
		vntr_bed = NULL;
		motif_csv = NULL;
		out_vcf = NULL;
		sampleName = NULL;
	};

	~IO() 
	{
		free(version);
		free(input_bam);
		free(vntr_bed);
		free(motif_csv);
		free(out_vcf);
		free(sampleName);
	};

	int readMotifsFromCsv (vector<VNTR *> &vntrs);

	int read_tsv(vector<vector<string>> &items);

	void readVNTRFromBed (vector<VNTR *> &vntrs);

	/* get the sequences from input_bam_file that overlapping with chr:start-end */
	void readSeqFromBam (vector<VNTR *> &vntrs, int nproc, int cur_thread, int sz);

	// void readSeqFromBam (vector<READ*> &reads, string &chr, const uint32_t &ref_VNTR_start, 
 //                       const uint32_t &ref_VNTR_end, const uint32_t &VNTR_len, string &region);
	
	// int outputVCF (vector<VNTR *> &vntrs);

	int writeVCFHeader(ofstream& out);

	int writeVCFBody(ofstream& out, vector<VNTR *> &vntrs, int tid, int nproc);
};

#endif
