#ifndef VNTR_H_
#define VNTR_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include "read.h"
using namespace std;

/*
class MOTIF contains:
@seq: the sequence of the motif
@len: the lenght of the motif
*/
class MOTIF
{
public:
	string seq;
	int len;
	MOTIF () {};
	MOTIF(string &Seq) : seq(Seq) { len = Seq.length();};
	~MOTIF() {};
};

/*
class VNTR contains:
@chr: the chromosome 
@start: the chromosome start coordinate
@end: the chromosome end coordinate
@reads: a group of sequences overlapping with the current VNTR locus
@annos: the index of the motifs annotation for each sequence, 
	 with function `commaSeparatedStringannotation`, you can get a string annotation
@concensus: the index of the concensus motifs annotation for the current VNTR locus
*/

class VNTR
{
public: 
	uint32_t ref_start;
	uint32_t ref_end;
	uint32_t len;
	string chr;
	string region;
	vector<MOTIF> motifs;
	vector<READ *> reads; 
	vector<vector<int>> annos; // the motif annotation for each read sequence
	vector<int> concensus_h1; // diploid genome
	vector<int> concensus_h2;

	VNTR () {};

	VNTR (string Chr, uint32_t Start, uint32_t End, uint32_t Len) : chr(Chr), ref_start(Start), ref_end(End), len(Len) 
	{
		string s = ":" + to_string(ref_start);
		string e = "-" + to_string(ref_end);
		region = Chr + s + e;
	};

	~VNTR () {};

	void clear ()
	{
		for (size_t i = 0; i < reads.size(); ++i) 
		{ 
			delete reads[i];
		}
		reads.clear();
		return;
	}

	/* for each sequence, get the annotation of motifs */
	void motifAnnoForOneVNTR ();

	void annoTostring (vector<int> &anno, string &annostr);

	/* for all the sequences at the current VNTR locus, get the concensus annotation */
	void concensusMotifAnnoForOneVNTR ();

	/* for one concensus annotation, output the comma-delimited annotation */
	void commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_rep);
};


#endif
