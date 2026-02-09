#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "option.h"
#include <math.h>
using namespace std;



void string_decomposer(vector<uint8_t> &optMotifs, vector<int> &optMotifStarts, vector<int> &optMotifEnds,
		       vector<int> &motifQV, vector<int> &starts, vector<int> &ends,
		       vector<vector<int> > &motifNMatches, vector<MOTIF> &motifs, const char *vntr, int N, const OPTION &opt,
		       SDTables &sdTables, vector<int > &mismatchCI)
{

    string vntrS(vntr, N);
    sdTables.Init(motifs, vntrS);
    int left=0;
    int down=1;
    int diag=2;
    int diagMis=3;    
    int endPath=4;
    int nInf=-99999999;

    // First column is a gap. First row is handled by a separate vector
    for (auto m=0; m < motifs.size(); m++)
    {
        int curGap=-opt.penalty_indel;
        for (auto mi=0; mi < motifs[m].len; mi++)
        {
            sdTables.scoreMat[m][0][mi] = curGap;
            sdTables.pathMat[m][0][mi]  = down;
            curGap -= opt.penalty_indel;
        }
    }

    int globMaxScore=nInf;
    int globMaxScoreRow=-1;
    int globMaxScoreCol=-1;
    int globMaxScoreMatrix=-1;
    for (auto s=0; s < N; s++)
    {
        for (auto m=0; m < motifs.size(); m++)
        {

            // Current column in matrices is s+1, string s.
            // Current row in block m is mi, string mi
            for (auto mi=0; mi < motifs[m].len; mi++)
            {
                int diagScore, insScore;
		bool match=true;
                if (mi == 0)
                {
                    int ms=0;
                    if (motifs[m].seq[0] == vntr[s]) { ms=opt.match;}
		    else { ms =-opt.penalty_mismatch; match=false;}
                    diagScore= sdTables.scoreRow0[s] + ms;
                    insScore = sdTables.scoreRow0[s] - opt.penalty_indel;
                }
                else
                {
                    int ms=0;
                    if (motifs[m].seq[mi] == vntr[s]) { ms=opt.match;}
		    else { ms =-opt.penalty_mismatch; match=false;}    
                    diagScore= sdTables.scoreMat[m][s][mi-1] + ms;
                    insScore = sdTables.scoreMat[m][s+1][mi-1] - opt.penalty_indel;
                }

                int delScore = sdTables.scoreMat[m][s][mi] - opt.penalty_indel; 
                int optScore=0;
                int optPath=0;
                if (diagScore >= insScore and diagScore >= delScore)
                {
                    optScore=diagScore;
		    if (match) 
		      optPath=diag;
		    else
		      optPath=diagMis;
                }
                else if (insScore >= diagScore and insScore >= delScore)
                {
                    optScore=insScore;
                    optPath=down;
                }
                else
                {
                    optScore=delScore;
                    optPath=left;
                }
		if (optScore > globMaxScore) {
		  globMaxScore = optScore;
		  globMaxScoreMatrix=m;
		  globMaxScoreRow=s+1;
		  globMaxScoreCol=mi;
		}
		if (optScore < 0 and opt.local) {
		  optScore = 0;
		  optPath=endPath;
		}
		
                sdTables.scoreMat[m][s+1][mi] = optScore;
                sdTables.pathMat[m][s+1][mi]  = optPath;
            }
        }
    
        // score loop back.

        int maxScore = sdTables.scoreMat[0][s+1][motifs[0].len-1];
        int maxIndex = 0;
        for (auto m=1; m < motifs.size(); m++)
        {
            if (sdTables.scoreMat[m][s+1][motifs[m].len-1] > maxScore )
            {
                maxScore = sdTables.scoreMat[m][s+1][motifs[m].len-1];
                maxIndex = m;
            }
        }

        sdTables.scoreRow0[s+1] = maxScore;
        sdTables.pathRow0[s+1]  = maxIndex;
	//
	// For local alignment, if the maximum score among all entries is zero
	// the alignment should terminate at this motif.
	//
	if (opt.local and maxScore == 0) {
	  sdTables.pathRow0[s+1]  = -1;
	}
	    

    }

    // Trace back;
    int maxScore=-9999999;
    int maxIndex;
    for (auto m=0; m < motifs.size(); m++)
    {
        if (m == 0 or sdTables.scoreMat[m][N][motifs[m].len-1] > maxScore )
        {
            maxScore = sdTables.scoreMat[m][N][motifs[m].len-1];
            maxIndex = m;
        }
    }
    if (opt.local and globMaxScore == 0) {
      return;
    }
    int cm, cs, cmi;
    if (opt.local == false) {
      cmi=motifs[maxIndex].len-1;
      cs=N;
      cm=maxIndex;
    }
    else {
      cm = globMaxScoreMatrix;
      cs = globMaxScoreRow;
      cmi= globMaxScoreCol;
    }
    ends.push_back(cs);
    optMotifs.push_back((uint8_t)cm);
    optMotifEnds.push_back(cmi+1);
    bool endTrace=false;
    int motifNMat=0;
    int motifNIndel=0;
    int motifNMismat=0;
    string motifStr="";
    string refStr="";
    
    string arrows[] = {"gmotif", "gref", "match", "mis", "pathend"};
    while (cs > 0 and cm >= 0 and endTrace == false)
    {
        assert(optMotifs.size() < 100000);
	assert(cmi >= 0);
	assert(cm >= 0);
	assert(cs >= 0);
        int arrow=sdTables.pathMat[cm][cs][cmi];
	//	cout << cm << " " << cs << " " << cmi << "\t" << arrows[arrow] << " " << motifs[cm].seq[cmi] << " " << vntr[cs-1] << endl;
        // At the beginning of a motif, may need to wrap around
        if (cmi == 0)
        {
            if (arrow == left)
            {
                assert(cs > 0);
                cs--;
            }
            else
            {
                cm = sdTables.pathRow0[cs-1];
		cmi=0;
                if (cs > 1 and cm >= 0)
                {
                cmi= motifs[cm].len-1;		  
		  if (motifNMat+motifNIndel+motifNMismat > 0) {
		    //		    cerr << "Motif score: " << motifNMat << " " <<  motifNIndel << " " << motifNMismat << " " << ((float)motifNMat ) / (motifNMat+motifNIndel+motifNMismat) << endl;
		  }
		  motifStr.push_back('+');
		  refStr.push_back('+');
		  
		  motifNMat = motifNIndel=motifNMismat=0;

		  optMotifStarts.push_back(0);
		  optMotifs.push_back((uint8_t)cm);
		  optMotifEnds.push_back(cmi+1);
		  ends.push_back(cs-1);
		  starts.push_back(cs-1);
                }
                cs--;
            }
        }
        else
        {
            if (arrow == diag or arrow== diagMis)
            {
	      motifStr.push_back(motifs[cm].seq[cmi]);
	      refStr.push_back(vntr[cs-1]);
                cs--;
                cmi--;
		if (arrow== diag) motifNMat++;
		if (arrow == diagMis) motifNMismat++;
            }
            else if (arrow == down)
            {
	      motifStr.push_back(motifs[cm].seq[cmi]);
	      refStr.push_back('-');
                cmi--;
		motifNIndel++;
            }
            else if (arrow == left)
            {
	      motifStr.push_back('-');
	      refStr.push_back(vntr[cs-1]);
                cs--;
		motifNIndel++;		
            }
	    else if (arrow == endPath) {
	      //
	      endTrace = true;
	    }
        }
        assert(cmi >= 0);
        assert(cm >=-1);
    }
    optMotifStarts.push_back(cmi); 
    starts.push_back(cs);
    //
    // Change for local alignment
    reverse(refStr.begin(), refStr.end());
    reverse(motifStr.begin(), motifStr.end());
    assert(starts.size() == ends.size());
    reverse(optMotifs.begin(), optMotifs.end());
    reverse(optMotifStarts.begin(), optMotifStarts.end());
    reverse(optMotifEnds.begin(), optMotifEnds.end());    
    reverse(starts.begin(), starts.end());
    reverse(ends.begin(), ends.end());
    reverse(motifNMatches.begin(), motifNMatches.end());
}

