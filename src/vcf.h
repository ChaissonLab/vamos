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

class OutWriter
{
public:
	char * sampleName;
	char * version;
	vector<string> target_names; // reference names
	vector<uint32_t> contigLengths;
	int32_t ncontigs;
	bool set;

	OutWriter () {
		sampleName = NULL;
		version = NULL;
		set = 0;
		ncontigs = 0;
	};

	~OutWriter () 
	{
		// free(version);
		// free(sampleName);
	};

	void init (char * input_bam_file, char * Version, char * SampleName);

    void writeHeader_locuswise(ofstream &out);

    void writeBody_locuswise(vector<VNTR *> &vntrs, ofstream &out, int tid, int nproc);

    void writeHeader_readwise(ofstream &out);

    void writeBody_readwise(vector<VNTR *> &vntrs, ofstream &out, int tid, int nproc);

	void writeNullAnno(vector<VNTR *> &vntrs, ofstream &out_nullAnno, int tid, int nproc);
};

#endif