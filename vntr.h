#ifndef VNTR_H_
#define VNTR_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include "io.h"

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

class VNTR () 
{
public: 
	char* chr;
	uint32_t ref_start;
	uint32_t ref_end;
	vector<string> motifs;
	vector<READ *> reads; 
	vector<vector<int>> reps;
	vector<int> concensus_h1; // diploid genome
	vector<int> concensus_h2;

	VNTR () {};
	VNTR (string Chr, uint32_t Start, uint32_t End) : ref_start(Start), ref_end(End) {chr = Chr.data()};
	~VNTR () {free();};

	/* get the sequences from input_bam_file that overlapping with chr:start-end */

	void readSeq (string & input_bam_file) 
	{
		reads.push_back(readSeqFromBam(input_bam_file));
		return;
	}

	/* for each sequence, apply Bida's code to get the representation of motifs */
	void motifRepresentationForOneSeq (const vector<string> &motifs, string &vntr) 
	{
		/* modify vector<vector<int>> reps
		   for now, just using dummy code
		*/
		reps.resize(reads.size());
		for (int i = 0; i < reads.size(); ++i)
		{
			reps[i].resize(10, 1);
		}
		return;
		
	}

	/* 
	for all the sequences at the current VNTR locus, get the concensus representation 
													 get the concensus 
	*/
	void concensusMotifRep ()
	{
		/* based on vector<int> rep, get a concensus representation of the current locus 
		   for now, just using some dummy code */
		concensus.resize(10, 1);
		return;
	}

	void commaSeparatedMotifRepForConsensus (bool h1, string & motif_rep)
	{
		if (h1)
			for (const auto &it : concensus_h1) { motif_rep += "VNTR_" + stoi(it) + ",";}
		else
			for (const auto &it : concensus_h2) { motif_rep += "VNTR_" + stoi(it) + ",";}
		if (!motif_rep.isempty())
			motif_ref.pop_back();
	}

	/* return the string representation for the current VNTR locus with comma as separator */
	void commaSeparatedMotifRepForOneSequence (int i)
	{
		vector<string> motif_rep;
		for (auto &it : reps[i])
		{
			motif_rep.push_back(stoi(it));
		} 
		return;
	}

	void commaSeparatedMotifRepForAllSequence (int i)
	{
		vector<vector<string>> motifreps;
		for (auto &seq : reps)
		{
			motifreps.push_back(commaSeparatedStringRepresentationForOneSequence(seq));
		}
		return;
	}

	void free()
	{
		for (uint32_t i = 0; i < reads.size(); ++i) 
		{
			delete reads[i];
		}
		reads.clear();
		return;
	}
};

/* read vntrs coordinates from file `input_vntr_bed`*/
void readVNTRFromBed (const string &input_vntr_bed, vector<VNTR> &vntrs)
{
	vector<vector<string>> items;
	read_tsv(input_vntr_bed, items);
	uint32_t start, end;
	string chr;
	for (const auto &it : items)
	{
		chr = it[0];
		start = stoi(it[1]);
		end = stoi(it[2]);
		vntrs.push_back(VNTR(chr, start, end));
	}
	return;
}


#endif
