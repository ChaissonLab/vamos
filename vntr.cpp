#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include "vntr.h"
#include "naive_anno.cpp"

void VNTR::motifAnnoForOneVNTR () 
{
	annos.resize(reads.size());
	for (int i = 0; i < reads.size(); ++i)
	{
		/* 
			apply Bida's code here
			input:   string : reads[i].seq; vector<MOTIF> : motifs
			output:  vector<vector<int>> annos[i]
		*/
		anno(annos[i], motifs, reads[i]->seq);
		cerr << "finish for reads: " << i << endl;
	}
	return;
}

void VNTR::concensusMotifAnnoForOneVNTR ()
{
	/* based on vector<vector<int>> annos, get a concensus representation of the current locus 
	   input:   vector<vector<int>> annos
	   output:  vector<int> concensus_h1, vector<int> concensus_h2
	*/
	// concensus_h1.resize(10, 1);
	// concensus_h2.resize(10, 1);


	
	return;
}

void VNTR::commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_rep)
{
	if (h1)
		for (auto &it : concensus_h1) { motif_rep += "VNTR_" + to_string(it) + ",";}
	else
		for (auto &it : concensus_h2) { motif_rep += "VNTR_" + to_string(it) + ",";}
	if (!motif_rep.empty()) motif_rep.pop_back();
}
