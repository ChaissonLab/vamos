#ifndef VCF_H_
#define VCF_H_
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
#include <ostream>
#include <stdlib.h> 
#include "vntr.h"
#include "htslib/sam.h"

using namespace std;

class VcfWriter
{
public:
	char * sampleName;
	char * version;
	vector<string> target_names; // reference names
	vector<uint32_t> contigLengths;
	int32_t ncontigs;
	bool set;

	VcfWriter () {
		sampleName = NULL;
		version = NULL;
		set = 0;
	};

	VcfWriter (char * input_bam_file, char * Version, char * SampleName) : version(Version), sampleName(SampleName)
	{
	    samFile * fp_in = hts_open(input_bam_file, "r"); //open bam file
	    bam_hdr_t * bamHdr = sam_hdr_read(fp_in); //read header

	    int32_t ncontigs = bamHdr->n_targets;

	    for (int32_t i = 0; i < ncontigs; ++i)
	    {
	    	contigLengths.push_back(bamHdr->target_len[i]);
	    	target_names.push_back(string(bamHdr->target_name[i]));
	    }

	    bam_hdr_destroy(bamHdr);
	    sam_close(fp_in); 
	    set = 1;
	};

	~VcfWriter () {};

    void writeHeader(ostream &out);

    void writeBody(vector<VNTR *> &vntrs, ostream &out);
};

#endif