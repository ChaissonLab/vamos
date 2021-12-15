#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include "vntr.h"
#include "naive_anno.cpp"
#include "seqan/sequence.h"
#include "seqan/align.h"
#include "seqan/score.h"
#include "edlib.h"
#include "dataanalysis.h"

void VNTR::motifAnnoForOneVNTR () 
{
	nreads = reads.size();
	if (nreads == 0) return;

	annos.resize(nreads);
	for (int i = 0; i < nreads; ++i)
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

void annoTostring_helper (string &annoStr, vector<int> &annoNum)
{
	int i;
	annoStr.clear();
	annoStr.resize(annoNum.size(), '+');
	for (i = 0; i < annoNum.size(); ++i)
	{
		// TODO: needs to remove '-' from encoding
		int num = annoNum[i];
		if (num >= 0 && num <= 222) {
			string tmp_s;
			tmp_s.assign(1, (char) (num + 33));
			annoStr.replace(i, 1, tmp_s);  
		}
	}
	return;
}

void VNTR::annoTostring ()
{
	if (nreads == 0) return;
	annoStrs.resize(nreads);
	int l, r, i, num;
	for (int r = 0; r < nreads; ++r)
	{
		annoTostring_helper(annoStrs[r], annos[r]);
	    cerr << "string: " << r << "  " << annoStrs[r];
	    cerr << endl;
	}
	return;
}

size_t VNTR::getAnnoStrLen (int i)
{
	assert(i < nreads);
	return annoStrs[i].length();
}

int VNTR::hClust (vector<int> &gp1, vector<int> &gp2, double edist [])
{
	/*
	 compute the edit distance matrix for annotation strings 
	 input: vector<string> annoStrs
	 output: double edist[nreads * nreads] = {...}
	*/
	int i, j;
	EdlibAlignResult result;
	for (i = 0; i < nreads; ++i)
	{
		for (j = i + 1; j < nreads; ++j) 
		{
			result = edlibAlign(annoStrs[i].c_str(), getAnnoStrLen(i), 
								annoStrs[j].c_str(), getAnnoStrLen(j), edlibDefaultAlignConfig());
            if (result.status == EDLIB_STATUS_OK) 
            {
                edist[i * nreads + j] = (double) result.editDistance; // [i][j]
                edist[j * nreads + i] = (double) result.editDistance; // [j][i]
            }
            else 
            {
                printf("edlib goes wrong!\n");
                return 1;
            }
		}
	}

	/* do hierarchical clustering based on the edit ditance matrix */
    alglib::clusterizerstate s;
    alglib::ahcreport rep;
    alglib::ae_int_t disttype;
    alglib::clusterizercreate(s);

    alglib::real_2d_array d;
    d.setcontent(nreads, nreads, edist);
    alglib::clusterizersetdistances(s, d, true);
    alglib::clusterizerrunahc(s, rep);

	/* get the cluster results when nclusts == 2 */
    alglib::integer_1d_array cidx;
    alglib::integer_1d_array cz;    
    clusterizergetkclusters(rep, 2, cidx, cz);

    gp1.clear(); gp2.clear();
    for (i = 0; i < nreads; ++i)
    {
    	if (cidx[i] == 0)
    		gp1.push_back(i);
    	else
    		gp2.push_back(i);
    }
    return 0;
}

void outputConsensus (vector<int> &consensusStr)
{
    cerr << "consensus sequence is\n";
    for (auto &it : consensusStr)
		cerr << to_string(it) << "->";
	cerr << endl;
	return;
}

void MSA_helper (int sampleSz, seqan::StringSet<seqan::String<char>> &annoSet, vector<int> &consensus)
{
    seqan::Align<seqan::String<char>> align(annoSet); // Initialize the Align object using a StringSet.
	seqan::Score<int> scoreScheme(1, -1, -1, -1); 
	int score = seqan::globalAlignment(align, scoreScheme);  // Compute a global alingment using the Align object.
    // int score = seqan::globalAlignment(align, seqan::EditDistanceScore());  // Compute a global alingment using the Align object.
    cerr << "score = " << score << endl;
    cerr << "align\n" << align << endl;

    /*
     create the profile string
     profile: row: alignment position
     		  col: reads
    */
    int r, idx; uint32_t i;
    uint32_t alnSz = seqan::length(seqan::row(align, 0));
    seqan::String<seqan::ProfileChar<char>> profile; 
    seqan::resize(profile, alnSz); 

    for (int r = 0; r < sampleSz; ++r) 
    {
        for (i = 0; i < alnSz; ++i) 
            profile[i].count[seqan::getValue(seqan::row(align, r), i) - 33] += 1;
    }

    for (i = 0; i < alnSz; ++i)
    {
        idx = seqan::getMaxIndex(profile[i]);
        consensus.push_back(idx);
    }
    return;
}

void MSA (const vector<int> &gp, const vector<string> &annoStrs, vector<int> &consensus)
{
	if (gp.empty()) return;

    seqan::StringSet<seqan::String<char>> annoSet;
    for (auto &r : gp)
    	seqan::appendValue(annoSet, annoStrs[r]);
    MSA_helper (gp.size(), annoSet, consensus);
	return;
}

void MSA (const vector<string> &annoStrs, vector<int> &consensus)
{
    seqan::StringSet<seqan::String<char>> annoSet;
    for (auto &str : annoStrs)
    	seqan::appendValue(annoSet, str);
    MSA_helper (annoStrs.size(), annoSet, consensus);
	return;
}

double avgSilhouetteCoeff (int nreads, double edist [], vector<int> &gp1, vector<int> &gp2)
{
	double sumSilC = 0;
	EdlibAlignResult result;
	double a1 = 0.0;
	double b1 = 10000000.0;
	double sc1 = 0.0;

	double a2 = 0.0;
	double b2 = 10000000.0;
	double sc2 = 0.0;

	for (auto &i : gp1)
	{
		a1 = 0.0; b1 = 10000000.0;
		for (auto &j : gp1) 
		{
			if (i == j) continue;
			a1 += edist[i * nreads + j];
		}

		if (gp1.size() > 1)
			a1 /= gp1.size() - 1;

		for (auto &j : gp2)
			b1 = min(edist[i * nreads + j], b1);

		sc1 += (b1 - a1) / max(a1, b1);
	}

	for (auto &i : gp2)
	{
		a2 = 0.0; b2 = 10000000.0;
		for (auto &j : gp2) 
		{
			if (i == j) continue;
			a2 += edist[i * nreads + j];
		}

		if (gp2.size() > 1)
			a2 /= gp2.size() - 1;

		for (auto &j : gp1)
			b2 = min(edist[i * nreads + j], b2);

		sc2 += (b2 - a2) / max(a2, b2);
	}

	return (sc1 + sc2) / (gp1.size() + gp2.size());
}

/* based on vector<vector<int>> annos, get a consensus representation of the current locus 
   input:   vector<vector<int>> annos
   output:  vector<int> consensus_h1, vector<int> consensus_h2
*/
void VNTR::concensusMotifAnnoForOneVNTR ()
{
	/* compute the MSA for each cluster when nclusts == 2 */
	vector<int> gp1, gp2;
	double edist [nreads * nreads];
	hClust(gp1, gp2, edist); 

	/* compare which is better: nclusts == 2 or nclusts == 1 */
	double sc_2clusts = avgSilhouetteCoeff(nreads, edist, gp1, gp2);
	if (sc_2clusts > 0.0)
	{
		het = true;
		MSA (gp1, annoStrs, consensus_h1);
		MSA (gp2, annoStrs, consensus_h2);
	}
	else 
	{
		het = false;
		MSA (annoStrs, consensus);
	}
	return;
}

void VNTR::commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_rep)
{
	if (h1)
		for (auto &it : consensus_h1) { motif_rep += "VNTR_" + to_string(it) + ",";}
	else
		for (auto &it : consensus_h2) { motif_rep += "VNTR_" + to_string(it) + ",";}
	if (!motif_rep.empty()) motif_rep.pop_back();
	return;
}