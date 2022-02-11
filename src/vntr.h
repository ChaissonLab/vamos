#ifndef VNTR_H_
#define VNTR_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "read.h"
#include "option.h"

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


class SDTables
{
public:
  vector<vector<vector< int > > > pathMat;
  vector<vector<vector< int > > > nMatchMat;
  vector<vector<vector< int > > > nDelMat;  
  vector<vector<vector< int > > > scoreMat;
  vector<int> pathRow0;
  vector<int> scoreRow0;
  vector<vector<int > > nMatchOnOptPath;
  vector<vector<int > > nDelOnOptPath;  

  void Init(vector<MOTIF> &motifs, string &seq) 
  {
	pathMat.resize(motifs.size());
	nMatchMat.resize(motifs.size());
	nDelMat.resize(motifs.size());    
	scoreMat.resize(motifs.size());

	for (auto m=0; m < motifs.size(); m++) 
	{
		pathMat[m].resize(seq.size()+1);
		nMatchMat[m].resize(seq.size()+1);
		nDelMat[m].resize(seq.size()+1);      
		scoreMat[m].resize(seq.size()+1);
		for (auto s=0; s < seq.size() + 1; s++) 
		{
			pathMat[m][s].resize(motifs[m].len);
			nMatchMat[m][s].resize(motifs[m].len);
			nDelMat[m][s].resize(motifs[m].len);	
			scoreMat[m][s].resize(motifs[m].len);
			fill(pathMat[m][s].begin(), pathMat[m][s].end(), 0);
			fill(nMatchMat[m][s].begin(), nMatchMat[m][s].end(), 0);
			fill(nDelMat[m][s].begin(), nDelMat[m][s].end(), 0);	
			fill(scoreMat[m][s].begin(), scoreMat[m][s].end(), 0);	
		}
	}
    pathRow0.resize(seq.size() + 1);
    scoreRow0.resize(seq.size() + 1);    
  }

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
	int len; // the ref len
	int cur_len; // the sample seq len
	string chr;
	string region;
	vector<MOTIF> motifs;
	vector<READ *> reads; 
	vector<vector<uint8_t>> annos; // the motif annotation for each read sequence
	vector<string> annoStrs; 
	int nreads;

	vector<READ *> clean_reads; 
	vector<vector<uint8_t>> clean_annos; 
	vector<string> clean_annoStrs;
	int ncleanreads; 
	vector<vector<double>> clean_edist;

	vector<vector<uint8_t>> consensus; 
	bool het;
	bool skip;
	bool nullAnno;
	vector<bool> nullAnnos;
   
	VNTR () { het = false; nreads = 0; skip = false; nullAnno = false;};

	VNTR (string Chr, uint32_t Start, uint32_t End, uint32_t Len) : chr(Chr), ref_start(Start), ref_end(End), len(Len) 
	{
		string s = ":" + to_string(ref_start);
		string e = "-" + to_string(ref_end);
		region = Chr + s + e;
		het = false; 
		nreads = 0;
		skip = false;
		nullAnno = false;
	};

	~VNTR () {};

	void clearRead ()
	{
		for (size_t i = 0; i < reads.size(); ++i) 
		{ 
			delete reads[i];
		}
		reads.clear();
		return;
	}

	/* for each sequence, get the annotation of motifs */
        void motifAnnoForOneVNTR (const OPTION &opt, SDTables &sdTables, vector<int > &mismatchCI);
	// string * getAnnoStr (int i);

	size_t getAnnoStrLen (int i);

	void annoTostring (const OPTION &opt);

	/* for all the sequences at the current VNTR locus, get the consensus; annotation */
	void concensusMotifAnnoForOneVNTR (const OPTION &opt);

	void concensusMotifAnnoForOneVNTRByABpoa (const OPTION &opt);

	void concensusMotifAnnoForOneVNTRBySeqan (const OPTION &opt);

	int hClust (vector<int> &gp1, vector<int> &gp2, double dists []);

	/* for one consensus; annotation, output the comma-delimited annotation */
	void commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_rep);

	int cleanNoiseAnno(const OPTION &opt);
};


void outputConsensus (vector<uint8_t> &consensus);


class Order {
public:
	vector<double> *dist;
	vector<int> index;

	Order() {dist = NULL;};
	Order(vector<double> *a): dist(a) {};

	void Update(vector<double> *a) 
	{
		dist = a;
		index.resize(a->size());
		for (int i = 0; i < index.size(); i++) index[i] = i;
		Sort();
	}

	int operator()(const int i, const int j) 
	{
		return (*dist)[i] < (*dist)[j];
	}

	void Sort() 
	{
		std::sort(index.begin(), index.end(), *this);
	}

	double & operator[](int i) 
	{
		return (*dist)[index[i]];
	}

	int size() 
	{
		return index.size();
	}
};


#endif
