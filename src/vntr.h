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
@consensus;: the index of the consensus; motifs annotation for the current VNTR locus
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
	vector<vector<uint8_t>> annos; // the motif annotation for each read sequence
	vector<string> annoStrs; 
	vector<uint8_t> consensus_h1; // diploid genome
	vector<uint8_t> consensus_h2;
	vector<uint8_t> consensus;
	int nreads;
	bool het;
	bool naive;

	VNTR () { het = false; nreads = 0; };

	VNTR (string Chr, uint32_t Start, uint32_t End, uint32_t Len) : chr(Chr), ref_start(Start), ref_end(End), len(Len) 
	{
		string s = ":" + to_string(ref_start);
		string e = "-" + to_string(ref_end);
		region = Chr + s + e;
		het = false; 
		nreads = 0;
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
	void motifAnnoForOneVNTR (bool naiveAnnoAlg);

	// string * getAnnoStr (int i);

	size_t getAnnoStrLen (int i);

	void annoTostring ();

	/* for all the sequences at the current VNTR locus, get the consensus; annotation */
	void concensusMotifAnnoForOneVNTR ();

	int hClust (vector<int> &gp1, vector<int> &gp2, double edist []);

	/* for one consensus; annotation, output the comma-delimited annotation */
	void commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_rep);
};


void outputConsensus (vector<uint8_t> &consensus);

#endif
