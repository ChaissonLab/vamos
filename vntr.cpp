#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "htslib/htslib/sam.h"
#include "vntr.h"

void VNTR::motifRepresentationForOneVNTR () 
{
	/* modify vector<vector<int>> reps
	   for now, just using dummy code
	*/
	reps.resize(reads.size());
	for (int i = 0; i < reads.size(); ++i)
	{
		/* 
			apply Bida's code here
			input: string : reads[i].seq
				  vector<string> : motifs
			output: vector<vector<int>> reps[i]
		*/
		reps[i].resize(10, 1);
	}
	return;
	
}

void VNTR::concensusMotifRepForOneVNTR ()
{
	/* based on vector<int> rep, get a concensus representation of the current locus 
	   for now, just using some dummy code 
	   input: vector<vector<int>> reps
	   output: vector<int> concensus_h1
	           vector<int> concensus_h2
	*/
	concensus_h1.resize(10, 1);
	concensus_h2.resize(10, 1);
	return;
}

void VNTR::commaSeparatedMotifRepForConsensus (bool h1, string &motif_rep)
{
	if (h1)
	{
		for (auto &it : concensus_h1) { motif_rep += "VNTR_" + to_string(it) + ",";}
	}
	else
	{
		for (auto &it : concensus_h2) { motif_rep += "VNTR_" + to_string(it) + ",";}
	}
	if (!motif_rep.empty()) motif_rep.pop_back();
}
