#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <numeric>
#include "vntr.h"
// #include "bounded_anno.cpp"
#include "bounded_anno_local.cpp"
#include "edlib.h"
#include "dataanalysis.h"
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


void outputConsensus (vector<uint8_t> &consensus, const OPTION &opt)
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


void VNTR::consensusReadForHapByABpoa (const OPTION &opt) 
{
    if (skip) return;
    if (nreads == 0) return;

    het = false;
    for (int i = 0; i < nreads; ++i) 
    {
        if (reads[i]->haplotype != 0) het = true;
        if (reads[i]->haplotype == 1) h1_reads.push_back(i);
        else if (reads[i]->haplotype == 2) h2_reads.push_back(i);
    }

    // hom, get the one consensus 
    if (!het) 
    {
        int motifs_size = motifs.size();
        h1_reads.resize(nreads);
        iota (begin(h1_reads), end(h1_reads), 0); // 0 to n_seqs - 1

        READ * read = new READ();
        char qname_1[] = "h1";
        read->qname = (char *) malloc(2 + 1);
        strcpy(read->qname, qname_1);
        Hap_seqs.push_back(read);

        MSA * msa = new MSA(h1_reads, reads);
        msa->MSA_seq_group (h1_reads, reads, Hap_seqs[0]);

        delete msa;

        // cerr.write(Hap_seqs[0]->seq, Hap_seqs[0]->len);
        // cerr << endl;
    }
    else
    {
        int motifs_size = motifs.size();
        assert (h1_reads.size() > 0 and h2_reads.size() > 0);
        READ * read_h1 = new READ();
        READ * read_h2 = new READ();

        char qname_1[] = "h1";
        char qname_2[] = "h2";

        read_h1->qname = (char *) malloc(2 + 1);
        strcpy(read_h1->qname, qname_1);
        read_h2->qname = (char *) malloc(2 + 1);
        strcpy(read_h2->qname, qname_2);        

        Hap_seqs.push_back(read_h1);
        Hap_seqs.push_back(read_h2);

        MSA * msa_1 = new MSA(h1_reads, reads);
        MSA * msa_2 = new MSA(h2_reads, reads);

        msa_1->MSA_seq_group (h1_reads, reads, Hap_seqs[0]);
        msa_2->MSA_seq_group (h2_reads, reads, Hap_seqs[1]);

        delete msa_1;
        delete msa_2;

        // cerr.write(Hap_seqs[0]->seq, Hap_seqs[0]->len);
        // cerr.write(Hap_seqs[1]->seq, Hap_seqs[1]->len);
        // cerr << endl;
    }

    return;
}


void VNTR::consensusMotifAnnoForOneVNTRByABpoa (const OPTION &opt)
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

    MSA * msa = new MSA(annos, motifs_size);   
    msa->MSA_anno_group (motifs.size(), annos, consensus[0]);
    assert(consensus[0].size() > 0);
    delete msa;

    return;
}

void VNTR::motifAnnoForOneVNTR (const OPTION &opt, SDTables &sdTables, vector<int > &mismatchCI) 
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

            if (consensus[i].size() > 20000) {
                skip = true;
                cerr << "skip the vntr" << endl;
                return;
            }    
            if (naive_flag) {
	      cerr << "This option is no longer supported" << endl;
	      exit(0);
	    }
	      //                naive_anno(consensus[i], motifs, Hap_seqs[i]->seq, Hap_seqs[i]->len);
            else {
                vector<int> starts, ends, sdAnnos, sdQV;
                vector<vector<int > > motifNMatch;
                // cout << "decomposing " << region << " " << i << " " << reads[i]->len << " " << motifs.size() << endl;
                string_decomposer(consensus[i], sdQV, starts, ends, motifNMatch,
                        motifs, Hap_seqs[i]->seq, Hap_seqs[i]->len, opt, sdTables, mismatchCI);
                
                /*
                cout << "heuristic: " << annos[i].size() << " sd " << sdAnnos.size() << endl;
                for (auto j=0; j< annos[i].size(); j++) 
                    cout << std::setw(3) << (int)annos[i][j] << ",";
                cout << endl;

                for (auto j=0; j < sdQV.size(); j++) 
                    cout << std::setw(3) << sdQV[j] << ",";
                cout << endl;

                for (auto j=0; j < sdAnnos.size(); j++) 
                    cout << std::setw(3) << sdAnnos[j] << ",";
                cout << endl;
                */
            }
            
            if (consensus[i].size() == 0) 
            {
                skip = true;
                nullAnno = true;
                nullAnnos[i] = true;
            } 
            else if (debug_flag) {
                cerr << endl;
                cerr.write(Hap_seqs[i]->qname, Hap_seqs[i]->l_qname);
                cerr << endl;
                cerr.write(Hap_seqs[i]->seq, Hap_seqs[i]->len);
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
        annos.resize(nreads);
        nullAnnos.resize(nreads, false);
        if (nreads > 200) {
          cerr << "WARNING, skipping locus " << region << " because it may be a centromeric tandem repeat" << endl;
          return;
        }

        for (int i = 0; i < nreads; ++i)
        {
            /* 
                apply Bida's code here
                input:   string : reads[i].seq; vector<MOTIF> : motifs
                output:  vector<vector<int>> annos[i]
            */

            if (reads[i]->len > 20000) {
                skip = true;
                cerr << "skip the vntr" << endl;
                return;
            }     

            if (naive_flag) {
	      cerr << "Naive annotation is not supported" << endl;
	      exit(0);
	    }//                naive_anno(annos[i], motifs, reads[i]->seq, reads[i]->len);
            else {
                // bounded_anno(annos[i], motifs, reads[i]->seq, reads[i]->len, opt);
                vector<int> starts, ends, sdAnnos, sdQV;
                vector<vector<int > > motifNMatch;
                // cout << "decomposing " << region << " " << i << " " << reads[i]->len << " " << motifs.size() << endl;
                string_decomposer(annos[i], sdQV, starts, ends, motifNMatch,
                        motifs, reads[i]->seq,  reads[i]->len, opt, sdTables, mismatchCI);
                
                /*
                cout << "heuristic: " << annos[i].size() << " sd " << sdAnnos.size() << endl;
                for (auto j=0; j< annos[i].size(); j++) 
                    cout << std::setw(3) << (int)annos[i][j] << ",";
                cout << endl;

                for (auto j=0; j < sdQV.size(); j++) 
                    cout << std::setw(3) << sdQV[j] << ",";
                cout << endl;

                for (auto j=0; j < sdAnnos.size(); j++) 
                    cout << std::setw(3) << sdAnnos[j] << ",";
                cout << endl;
                */
            }
            
            if (annos[i].size() == 0) 
            {
                skip = true;
                nullAnno = true;
                nullAnnos[i] = true;
            } 
            else if (debug_flag) {
                cerr << endl;
                cerr.write(reads[i]->qname, reads[i]->l_qname);
                cerr << endl;
                cerr.write(reads[i]->seq, reads[i]->len);
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

void VNTR::commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_anno)
{
	if (h1) {
        for (auto &it : consensus[0]) { 
            motif_anno += to_string(it) + ",";
            len_h1 += 1;
        }        
    }
	else if (het) {
        for (auto &it : consensus[1]) { 
            motif_anno += to_string(it) + ",";
            len_h2 += 1;
        }       
    }
	else {
        for (auto &it : consensus[0]) { 
            motif_anno += to_string(it) + ",";
            len_h2 += 1;
        }        
    }

	if (!motif_anno.empty()) motif_anno.pop_back();
	return;
}
