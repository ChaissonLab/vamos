#ifndef VNTR_H_
#define VNTR_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "htslib/sam.h"
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
@reps: the index of the motifs representation for each sequence, 
	 with function `commaSeparatedStringRepresentation`, you can get a string representation
@concensus: the index of the concensus motifs representation for the current VNTR locus
*/

class VNTR
{
public: 
	char * chr;
	uint32_t ref_start;
	uint32_t ref_end;
	uint32_t len;
	char * region;
	vector<MOTIF> motifs;
	vector<READ *> reads; 
	vector<vector<int>> reps; // the motif representation for each read sequence
	vector<int> concensus_h1; // diploid genome
	vector<int> concensus_h2;

	VNTR () {};

	VNTR (string Chr, uint32_t Start, uint32_t End, uint32_t Len) : ref_start(Start), ref_end(End), len(Len) 
	{
		chr = (char *) malloc(Chr.length() + 1);
		strcpy(chr, Chr.c_str());

		string s = ":" + to_string(ref_start);
		string e = "-" + to_string(ref_end);
		int sz = Chr.length() + s.length() + e.length(); // c string : original string + "\0"
		char * region = (char *) malloc(sz + 1);
		strcpy(region, chr);
		strcat(region, s.c_str());
		strcat(region, e.c_str());

		printf("region %s", region); 
	};

	~VNTR () 
	{
		free(chr);
		free(region);
		for (uint32_t i = 0; i < reads.size(); ++i) 
		{
			delete reads[i];
		}
		reads.clear();
		return;
	};

	/* get the sequences from input_bam_file that overlapping with chr:start-end */
	// void readSeq (string & input_bam_file);

	/* for each sequence, get the representation of motifs */
	void motifRepresentationForOneVNTR ();

	/* 
	for all the sequences at the current VNTR locus, get the concensus representation 
													 get the concensus 
	*/
	void concensusMotifRepForOneVNTR ();

	void commaSeparatedMotifRepForConsensus (bool h1, string &motif_rep);

};


#endif
