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
// #include <seqan/align.h>
// #include <seqan/graph_msa.h>

extern int naive_flag;
extern int debug_flag;
extern int hclust_flag;
// extern int seqan_flag;
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


void VNTR::consensusReadForHapByABpoa(const OPTION &opt) 
{
    if (skip) return;
    if (reads.size() == 0) return;
    het = false;
    assert(h1_reads.size() == 0);
    assert(h2_reads.size() == 0);
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
        int motifs_size = motifs.size();
        h1_reads.resize(reads.size());
        iota (begin(h1_reads), end(h1_reads), 0); // 0 to n_seqs - 1

        READ * read = new READ();
        char qname_1[] = "h1";
        read->qname = "h1"; //(char *) malloc(2 + 1);
	//        strcpy(read->qname, qname_1);
        Hap_seqs.push_back(read);
        // msa to get the consensus

	MSA msa(h1_reads, reads);
	msa.MSA_seq_group(h1_reads, reads, Hap_seqs[0]);

	//	cerr << ">h1" << endl;
	//	cerr << read->seq <<endl;

    }
    // heterozygous, get the two consensuses
    else
    {
        int motifs_size = motifs.size();

        if (h1_reads.size() > 0)
        {
	    int nLoop= 0;
	    //	    if (reads.size() > 0 and reads[0]->seq.size() > 20) {

	    READ * read_h1 = new READ();
	    char qname_1[] = "h1";      
	    read_h1->qname = "h1"; //(char *) malloc(2 + 1);
	    //            strcpy(read_h1->qname, qname_1);
	    //            read_h1->qname[2] = '\0';
	    Hap_seqs.push_back(read_h1);
	    // msa to get the consensus
	    MSA msa_1(h1_reads, reads);
	    msa_1.MSA_seq_group(h1_reads, reads, Hap_seqs[0]);
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

            // msa to get the consensus
            MSA msa_2(h2_reads, reads);
	    //            assert(1 < Hap_seqs.size());
            msa_2.MSA_seq_group(h2_reads, reads, Hap_seqs[Hap_seqs.size()-1]);
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


void VNTR::consensusMotifAnnoForOneVNTRByABpoa(const OPTION &opt)
{
    if (skip) return;
    int n_seqs = annos.size();
    int motifs_size = motifs.size();
    if (n_seqs == 0) return;

    consensus.resize(1);
    if (n_seqs == 1)
    {
        consensus[0] = annos[0];
        return;
    }

    MSA msa(annos, motifs_size);
    msa.MSA_anno_group(motifs.size(), annos, consensus[0]);
    assert(consensus[0].size() > 0);

    return;
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
            // naive_anno(consensus[i], motifs, Hap_seqs[i]->seq, Hap_seqs[i]->len);
            else
            {
                vector<int> starts, ends, sdAnnos, sdQV;
                vector<vector<int > > motifNMatch;
		if (Hap_seqs[i]->seq.size() > 0) {
		  string_decomposer(consensus[i], sdQV, starts, ends, motifNMatch,
				    motifs, Hap_seqs[i]->seq.c_str(), Hap_seqs[i]->len, opt, sdTables, mismatchCI);
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
        nullAnnos.resize(reads.size(), false);
        if (reads.size() > 200)
        {
          cerr << "WARNING, skipping locus " << region << " because it may be a centromeric tandem repeat" << endl;
          return;
        }

        for (int i = 0; i < reads.size(); ++i)
        {
   	  haps[i] = reads[i]->haplotype;

            if (reads[i]->len > 20000)
            {
                skip = true;
                cerr << "skip the vntr (length > " << reads[i]->len << ")" << endl;
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
                string_decomposer(annos[i], sdQV, starts, ends, motifNMatch,
				  motifs, reads[i]->seq.c_str(),  reads[i]->len, opt, sdTables, mismatchCI);

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

    if (skip) cerr << "skip the vntr" << endl;
    return;
}


void VNTR::commaSeparatedMotifAnnoForConsensus(bool h1, string &motif_anno)
{
    if (h1)
      {
	if (consensus.size() > 0) {
	  for (auto &it : consensus[0])
	    {
	      motif_anno += to_string(it) + ",";
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
            motif_anno += to_string(it) + ",";
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
	      motif_anno += to_string(it) + ",";
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
