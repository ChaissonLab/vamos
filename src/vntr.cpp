#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <numeric>
#include "vntr.h"
#include "bounded_anno_local.cpp"
#include "abpoa.h"
#include "option.h"
#include "msa.h"
#include "FitLengthMixtures.h"

extern int naive_flag;
extern int debug_flag;
extern int hclust_flag;
extern int anno_flag;
extern int subseq_flag;
extern int locuswise_flag;


void outputConsensus(vector<uint8_t> &consensus, const OPTION &opt)
{
    size_t i = 0;
    for (auto &it : consensus)
    {
        if (debug_flag)
        {
            if (i < consensus.size() - 1)
                cerr << to_string(it) << ",";
            else
                cerr << to_string(it);
        }
        i += 1;
    }
    if (debug_flag) cerr << endl;
    return;
}

void SetPloidy(vector<VNTR*> &vntrs, int minChrY, int nChrY) {
  bool isXY = false;
  if (nChrY > minChrY) {
    isXY = true;
  }
  for (auto &v: vntrs) {
    if (isXY == true and (v->chr == "chrX" || v->chr == "X" || v->chr == "chrY" || v->chr == "Y")) {
      v->ploidy=1;
    }
    else {
      v->ploidy=2;
    }
  }
}

void VNTR::consensusReadForHapByABpoa(const OPTION &opt) 
{
    if (skip) return;
    if (reads.size() == 0) return;
    het = false;
    assert(h1_reads.size() == 0);
    assert(h2_reads.size() == 0);
    //    cout << "CONS:\t" << region << "\t";
    for (int i = 0; i < reads.size(); ++i) 
    {
      if (reads[i]->haplotype != 0) het = true;
      if (reads[i]->haplotype == 1) h1_reads.push_back(i);
      else if (reads[i]->haplotype == 2) h2_reads.push_back(i);
    }

    // if all reads to one hap, take as homozygous
    if ( (h1_reads.size() == 0 and h2_reads.size() > 0) or
        (h1_reads.size() > 0 and h2_reads.size() == 0) ) {het = false;}

    // homozygous, get the one consensus
    if (!het)
    {
      vector<int> readLengths;
      for (auto &s: reads) {
	readLengths.push_back(s->seq.size());
      }
      DetectionResult mmResult;
      mmResult.model = PoissonModel::Single;
      if (readLengths.size() > 2) {
	mmResult=detectPoissonModel(readLengths);
      }
      if (mmResult.model == PoissonModel::Mixture) {
	int nl1=0, nl2=0;
	for (auto rl: readLengths) {
	  if (WithinRange(rl, mmResult.mixture.lambda1, 0.9) and
	      !WithinRange(rl, mmResult.mixture.lambda2, 0.9)) {
	    nl1++;
	  }
	  else if (WithinRange(rl, mmResult.mixture.lambda2, 0.9) and
		   !WithinRange(rl, mmResult.mixture.lambda1, 0.9)) {
	    nl2++;
	  }
	}
	if (nl1 > 0.3 * readLengths.size() and nl2 > 0.3*readLengths.size()) {
	  //
	  // This is a heterozygous locus that is autozygous according to snps
	  //
	  het = true;	  
	  for (int i=0; i < reads.size(); i++) {
	    int readLen = reads[i]->seq.size();
	    
	    if (WithinRange(readLen, mmResult.mixture.lambda1, 0.9) and
		!WithinRange(readLen, mmResult.mixture.lambda2, 0.9)) {
	      h1_reads.push_back(i);
	    }
	    else if (WithinRange(readLen, mmResult.mixture.lambda2, 0.9) and
		     !WithinRange(readLen, mmResult.mixture.lambda1, 0.9)) {
	      h2_reads.push_back(i);
	    }
	  }
	}
      }
      //
      // Check again to see if this is heterozygous after checking mixture model partitions.
      // 
      if (het == false) {
        int motifs_size = motifs.size();
        h1_reads.resize(reads.size());
        iota (begin(h1_reads), end(h1_reads), 0); // 0 to n_seqs - 1

        READ * read = new READ();
        char qname_1[] = "h1";
        read->qname = "h1"; //(char *) malloc(2 + 1);
	//        strcpy(read->qname, qname_1);
        Hap_seqs.push_back(read);
        // msa to get the consensus
	int totalReadSize=0;
	for (int rit =0; rit <reads.size(); rit++) {
	  totalReadSize+=reads[rit]->seq.size();
	}
	MSA msa(h1_reads, reads);
	msa.MSA_seq_group(h1_reads, reads, Hap_seqs[0], opt.pruneExtremities);
      }
    }
    // heterozygous, get the two consensuses
    if (het) {
      vector<int> h1rl, h2rl;
      for (auto &s: h1_reads) {
	h1rl.push_back(reads[s]->seq.size());
      }
      for (auto &s: h2_reads) {
	h2rl.push_back(reads[s]->seq.size());
      }

      DetectionResult mmResult;
      mmResult.model = PoissonModel::Single;
      bool usedMM =false;
      if (h1rl.size() > 2 and h2rl.size() > 2) {
	mmResult=detectPoissonModel(h1rl, h2rl);
	usedMM=true;
      }
	if (usedMM == true and mmResult.model == PoissonModel::Mixture) {
	  vector<int> h1RemoveI, h2RemoveI;	  
	  vector<int> h1RemoveLen, h2RemoveLen;
	  vector<int> h1RemoveIdx, h2RemoveIdx;	  
	  vector<int> h1l, h2l;
	  int nReassigned=0;
	  double h1Mean=0, h2Mean=0;
	  for (int i=0; i < h1_reads.size(); i++) {
	    int seqLen=reads[h1_reads[i]]->seq.size();
	    h1Mean+=seqLen;
	    h1l.push_back(seqLen);
	  }
	  if (h1_reads.size() > 0) { h1Mean = h1Mean/h1_reads.size();}
	  for (int i=0; i < h2_reads.size(); i++) {
	    int seqLen=reads[h2_reads[i]]->seq.size();
	    h2Mean+=seqLen;
	    h2l.push_back(seqLen);	    
	  }
	  if (h2_reads.size() > 0) { h2Mean = h2Mean/h2_reads.size();}

	  if (fabs(mmResult.mixture.lambda2 - h1Mean)  < fabs(mmResult.mixture.lambda1 - h1Mean)) {
	    auto temp=mmResult.mixture.lambda2;
	    mmResult.mixture.lambda2 = mmResult.mixture.lambda1;
	    mmResult.mixture.lambda1=temp;
	  }
	  int h1Size=0;
	  int h2Size=0;
	  for (int i=h1_reads.size(); i > 0; i--) {
	    int seqLen=reads[h1_reads[i-1]]->seq.size();
	    double dist1 = fabs(mmResult.mixture.lambda1-seqLen);
	    double dist2 = fabs(mmResult.mixture.lambda2-seqLen);
	    if (WithinRange((float)seqLen, mmResult.mixture.lambda2, 0.9) and
		! WithinRange((float)seqLen, mmResult.mixture.lambda1, 0.7)) {
	      h1RemoveLen.push_back(seqLen);
	      h1RemoveI.push_back(i-1);
	      h1RemoveIdx.push_back(h1_reads[i-1]);
	      ++nReassigned;
	    }
	    else {
	      h1Size++;
	    }
	  }
	  for (int i=h2_reads.size(); i > 0; i--) {
	    int seqLen=reads[h2_reads[i-1]]->seq.size();
	    double dist1 = fabs(mmResult.mixture.lambda1-seqLen);
	    double dist2 = fabs(mmResult.mixture.lambda2-seqLen);	    
	    if (WithinRange((float)seqLen, mmResult.mixture.lambda1, 0.9) and
		! WithinRange((float)seqLen, mmResult.mixture.lambda2, 0.7)) {
	      h2RemoveLen.push_back(seqLen);
	      h2RemoveI.push_back(i-1);
	      h2RemoveIdx.push_back(h2_reads[i-1]);
	      ++nReassigned;	      
	    }
	    else {
	      h2Size++;
	    }
	  }
	  if (nReassigned > 0 and h1Size > 0.5*h1_reads.size() and h2Size > 0.5*h2_reads.size()) {
	    // The clusters are of sufficient size to reassign repeats.
	    for (auto i: h1RemoveI) {
	      h1_reads.erase(h1_reads.begin() + i);
	    }
	    for (auto i: h2RemoveI) {
	      h2_reads.erase(h2_reads.begin() + i);
	    }
	    // Add indexes back to cluster.
	    for (auto idx: h1RemoveIdx) {
	      h2_reads.push_back(idx);
	    }
	    for (auto idx: h2RemoveIdx) {
	      h1_reads.push_back(idx);
	    }
	  }
	}	  
	
        int motifs_size = motifs.size();
        if (h1_reads.size() > 0)
        {
	    int nLoop= 0;
	    //	    if (reads.size() > 0 and reads[0]->seq.size() > 20) {

	    READ * read_h1 = new READ();
	    char qname_1[] = "h1";      
	    read_h1->qname = "h1"; //(char *) malloc(2 + 1);
	    Hap_seqs.push_back(read_h1);
	    // msa to get the consensus
	    
	    MSA msa_1(h1_reads, reads);
            msa_1.MSA_seq_group(h1_reads, reads, Hap_seqs[0], opt.pruneExtremities);	    
	    int totalReadSize=0;
	    //	    vector<int> rls;
	    for (int rit =0; rit <h1_reads.size(); rit++) {
	      totalReadSize+=reads[h1_reads[rit]]->seq.size();
	      //	      rls.push_back(reads[h1_reads[rit]]->seq.size());	      
	    }
	    
        }
	else {
	  if (het) {
	    READ * read_h1 = new READ();
	    char qname_1[] = "h1";
	    read_h1->seq = "";
	    read_h1->len = 0;
	    Hap_seqs.push_back(read_h1);
	  }
	}
        if (h2_reads.size() > 0)
        {
	  
            READ * read_h2 = new READ();
            char qname_2[] = "h2";
	    read_h2->qname = "h2"; //(char *) malloc(2 + 1);
	    //            strcpy(read_h2->qname, qname_2);
	    //            read_h2->qname[2] = '\0';
            Hap_seqs.push_back(read_h2);
	    int totalReadSize=0;
	    //	    vector<int> rls;
	    for (int rit =0; rit <h2_reads.size(); rit++) {
	      totalReadSize+=reads[h2_reads[rit]]->seq.size();
	      //	      rls.push_back(reads[h2_reads[rit]]->seq.size());
	    }
	    /*
	    cout << "Hap2 mean: " << totalReadSize / h2_reads.size() << endl;
	    sort(rls.begin(), rls.end());
	    for (auto r: rls) { cout << " " << r;} cout << endl;
	    */
            // msa to get the consensus
            MSA msa_2(h2_reads, reads);
	    //            assert(1 < Hap_seqs.size());
            msa_2.MSA_seq_group(h2_reads, reads, Hap_seqs[Hap_seqs.size()-1], opt.pruneExtremities);
        }
	else {
	  if (het) {
            READ * read_h2 = new READ();
            char qname_2[] = "h2";
	    read_h2->qname = "h2"; //(char *) malloc(2 + 1);
            Hap_seqs.push_back(read_h2);
	    read_h2->seq = "";
	    read_h2->len = 0;
	  }
	}

        // cerr.write(Hap_seqs[0]->seq, Hap_seqs[0]->len);
        // cerr.write(Hap_seqs[1]->seq, Hap_seqs[1]->len);
        // cerr << endl;
	if (het) {
	  assert(Hap_seqs.size() == 2);
	}
    }

    return;
}



void VNTR::ReconstructTRSequence(vector<uint8_t> &optMotifs, vector<int> &optMotifStarts, vector<int> &optMotifEnds, string &reconstructedSeq) {
  int v=0;
  reconstructedSeq="";
  for (v=0; v < optMotifs.size(); v++) {
    reconstructedSeq+=motifs[optMotifs[v]].seq.substr(optMotifStarts[v], optMotifEnds[v]-optMotifStarts[v]);
  }
}


void VNTR::motifAnnoForOneVNTR(const OPTION &opt, SDTables &sdTables, vector<int> &mismatchCI) 
{
    if (skip) return;

    if (locuswise_flag) 
    {
        // annotate for reads
        if (het) consensus.resize(2);
        else consensus.resize(1);

        int n_consensus = consensus.size();
        nullAnnos.resize(n_consensus, false);

        for (int i = 0; i < n_consensus; ++i)
        {
            // TODO: this statemnet should never be checked.
            if (Hap_seqs[i]->len == 0 or Hap_seqs[i]->len > 30000)
            {
                skip = true;
                return;
            }
            if (naive_flag)
            {
                cerr << "This option is no longer supported" << endl;
                exit(0);
            }
            else
            {
                vector<int> starts, ends, sdAnnos, sdQV;
		vector<int> motifStarts, motifEnds;
                vector<vector<int > > motifNMatch;
		if (Hap_seqs[i]->seq.size() > 0) {
		  string_decomposer(consensus[i], motifStarts, motifEnds, sdQV, starts, ends, motifNMatch,
				    motifs, Hap_seqs[i]->seq.c_str(), Hap_seqs[i]->len, opt, sdTables, mismatchCI);
		  string reconstructedSeq;
		  ReconstructTRSequence(consensus[i], motifStarts, motifEnds, reconstructedSeq);
		  reconstructedTRSeqs.push_back(reconstructedSeq);
		}
            }
            if (consensus[i].size() == 0) 
            {
                skip = true;
                nullAnno = true;
                nullAnnos[i] = true;
            } 
            else if (debug_flag)
            {
                cerr << endl;
		cerr << Hap_seqs[i]->qname;
                cerr << endl;
		cerr << Hap_seqs[i]->seq;
                cerr << endl;
                cerr << Hap_seqs[i]->len << "  " << i << endl;
                for (auto &idx : consensus[i]) cerr << int(idx) << ",";
                cerr << endl;
            }
        }
    }
    else  
    {
        // annotate for consensus read
        annos.resize(reads.size());
	haps.resize(reads.size(), 0);
	reconstructedTRSeqs.resize(reads.size());
        nullAnnos.resize(reads.size(), false);
        if (reads.size() > opt.maxCoverage)
        {
          cerr << "WARNING, skipping locus " << region << " because it may be a centromeric tandem repeat" << endl;
          return;
        }
        annos.resize(reads.size());
        nullAnnos.resize(reads.size(), false);

        for (int i = 0; i < reads.size(); ++i)
        {
   	  haps[i] = reads[i]->haplotype;

            if (reads[i]->len > opt.maxLocusLength)
            {
                skip = true;
                return;
            }     

            if (naive_flag)
            {
                cerr << "Naive annotation is not supported" << endl;
                exit(0);
            }
            else
            {
                vector<int> starts, ends, sdAnnos, sdQV;
                vector<vector<int> > motifNMatch;
		vector<int> motifStarts, motifEnds;		
                string_decomposer(annos[i],  motifStarts, motifEnds, sdQV, starts, ends, motifNMatch,
				  motifs, reads[i]->seq.c_str(),  reads[i]->len, opt, sdTables, mismatchCI);
		string reconstructedSeq;
		ReconstructTRSequence(annos[i], motifStarts, motifEnds, reconstructedSeq);
		reconstructedTRSeqs[i]=reconstructedSeq;
            }
            
            if (annos[i].size() == 0) 
            {
                skip = true;
                nullAnno = true;
                nullAnnos[i] = true;
            } 
            else if (debug_flag)
            {
                cerr << endl;
		cerr << reads[i]->qname;
                cerr << endl;		
                cerr << reads[i]->seq;
                cerr << endl;
                cerr << reads[i]->len << "  " << i << endl;
                for (auto &idx : annos[i]) cerr << int(idx) << ",";
                cerr << endl;
            }
        }
    }

    return;
}


void VNTR::commaSeparatedMotifAnnoForConsensus(bool h1, string &motif_anno)
{
    if (h1)
      {
	if (consensus.size() > 0) {
	  for (auto &it : consensus[0])
	    {
	      motif_anno += to_string(it) + "-";
	      len_h1 += 1;
	    }
	}
	else {
	  motif_anno = "DEL ";
	  len_h1 = 0;
	}
      }
    
    else if (het)
    {
      if (consensus.size() > 1) {
        for (auto &it : consensus[1])
	  {
            motif_anno += to_string(it) + "-";
            len_h2 += 1;
	  }
      }
      else {
	motif_anno="DEL ";
	len_h2 = 0;
      }
    }
    else
      {
	if (consensus.size() > 0) {
	  for (auto &it : consensus[0])
	    {
	      motif_anno += to_string(it) + "-";
	      len_h2 += 1;
	    }
	}
	else {
	  motif_anno = "DEL ";
	  len_h2 = 0;
	}
      }      

    if (!motif_anno.empty()) motif_anno.pop_back();
    return;
}

void SortVNTRIndexByPos(vector<VNTR*> &vntrs, vector<int> &idx) {
  CompareVNTRPos comp;
  comp.vntrs = &vntrs;
  sort(idx.begin(), idx.end(), comp);
}

void ChromToVNTRMap(vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap) {
  for (int i=0; i < vntrs.size(); i++ ){
    if (vntrMap.find(vntrs[i]->chr) == vntrMap.end()) {
      vntrMap[vntrs[i]->chr] = vector<int>();
    }
    vntrMap[vntrs[i]->chr].push_back(i);
  }

  for (auto it = vntrMap.begin(); it != vntrMap.end(); it++ ) {
    SortVNTRIndexByPos(vntrs, it->second);
  }
}
