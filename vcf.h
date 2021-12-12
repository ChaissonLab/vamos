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
	char ** target_names;
	uint32_t * contigLengths;
	int32_t ncontigs;
	bool set;

	VcfWriter () {set = 0;};
	VcfWriter (char * input_bam_file, char * Version, char * SampleName) : version(Version), sampleName(SampleName)
	{
	    samFile * fp_in = hts_open(input_bam_file, "r"); //open bam file
	    bam_hdr_t * bamHdr = sam_hdr_read(fp_in); //read header

	    int32_t ncontigs = bamHdr->n_targets;
	    target_names = (char **) malloc(ncontigs);
	    contigLengths = (uint32_t *) malloc(ncontigs);

	    for (int32_t i = 0; i < ncontigs; ++i)
	    {
	    	contigLengths[i] = bamHdr->target_len[i];
	    	target_names[i] = bamHdr->target_name[i];
	    }

	    bam_hdr_destroy(bamHdr);
	    sam_close(fp_in); 
	    set = 1;
	};

	~VcfWriter () 
	{
		free(target_names);
		free(contigLengths);
	};

    void writeHeader(ostream &out);
    void writeBody(vector<VNTR> &vntrs, ostream &out);
};

#endif