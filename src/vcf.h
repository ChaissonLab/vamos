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
#include <string>

using namespace std;

class OutWriter
{
public:
	string sampleName;
	string version;
	vector<string> target_names; // reference names
	vector<uint32_t> contigLengths;
	int32_t ncontigs;
	bool set;

	OutWriter () {
		sampleName = "";
		version = "";
		set = 0;
		ncontigs = 0;
	};

	~OutWriter () 
	{
		// free(version);
		// free(sampleName);
	};

	void init (string input_bam_file, string Version, string SampleName);

    void writeHeader_locuswise(ofstream &out);

    void writeBody_locuswise(vector<VNTR *> &vntrs, ofstream &out, int tid, int nproc);

    void writeHeader_readwise(ofstream &out);

    void writeBody_readwise(vector<VNTR *> &vntrs, ofstream &out, int tid, int nproc);

	void writeNullAnno(vector<VNTR *> &vntrs, ofstream &out_nullAnno, int tid, int nproc);
};

#endif
