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
#include "naive_anno.cpp"
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
    /*
    // sort annos
    vector<vector<uint8_t>> annos_cp(annos);
    std::sort(annos_cp.begin(), annos_cp.end(), []( vector<uint8_t> & a, vector<uint8_t> & b){ return a.size() > b.size(); });

    // collect sequence length and sequence
    int i, j;
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs);
    for (i = 0; i < n_seqs; ++i) {
        seq_lens[i] = annos_cp[i].size();
        bseqs[i] = (uint8_t*) malloc(sizeof(uint8_t) * (seq_lens[i] + 1));
        for (j = 0; j < seq_lens[i]; ++j) 
        {
            assert(annos_cp[i][j] < motifs_size); 
            bseqs[i][j] = annos_cp[i][j] + 1; // +1 to avoid 0
        }
        bseqs[i][seq_lens[i]] = '\0';
    }

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();

    // alignment parameters
    abpt->align_mode = 1;  // 0:global 1:local, 2:extension
    abpt->match = 1;       // match score
    abpt->mismatch = 1;    // mismatch penalty
    abpt->gap_mode = ABPOA_AFFINE_GAP; // gap penalty mode
    abpt->gap_open1 = 1;   // gap open penalty #1
    abpt->gap_ext1 = 1;    // gap extension penalty #1
    abpt->gap_open2 = 1;   // gap open penalty #2
    abpt->gap_ext2 = 1;    // gap extension penalty #2
                           // gap_penalty = min{gap_open1 + gap_len * gap_ext1, gap_open2 + gap_len * gap_ext2}
    abpt->is_diploid = 0;
    abpt->out_msa = 0; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->progressive_poa = 1;
    abpt->out_gfa = 0;
    abpt->ret_cigar = 0; // turn off the cigar string. Otherwise, cg_backtrack step crashes! 
    abpt->rev_cigar = 0;
    abpt->out_pog = 0;

    abpoa_post_set_para(abpt);

    // variables to store result
    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n = 0;
    uint8_t **msa_seq; int msa_l = 0;

    // perform abpoa-msa
    if (n_seqs == 2) {
      skip = true;
      return;
    }
    skip = true;
    for (auto i = 0; i < n_seqs; i++) {
      if (seq_lens[i] != 2) {
        skip = false;
      }
    }
    if (skip) {
      return;
    }

    

    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);
    // abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, NULL);//&cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);    
    // cons_n = ab->abc->n_cons;

    assert(cons_n > 0); // cons_n == 0 means no consensus sequence exists
    if (cons_n > 1) cerr << "het" << endl;


    int numHap = cons_n > 1 ? 2 : 1;
    consensus.resize(numHap);
    het = cons_n > 1 ? true : false;
    for (i = 0; i < numHap; ++i)
    {
        for (j = 0; j < cons_l[i]; ++j) 
        {
            assert(cons_seq[i][j] <= motifs_size);
            consensus[i].push_back((uint8_t)(cons_seq[i][j] - 1));            
        }
        // for (j = 0; j < ab->abc->cons_len[i]; ++j)
        // {
        //     assert(ab->abc->cons_base[i][j] <= motifs_size);
        //     consensus[i].push_back((uint8_t)(ab->abc->cons_base[i][j]) - 1);
        // }
	    // cerr << "got consensus of size " << consensus.size() << endl;	      
        assert(consensus[i].size() > 0);        
    }
    if (cons_n) {
        for (i = 0; i < cons_n; ++i) 
        {
            free(cons_seq[i]); 
            free(cons_cov[i]);
        } 
        free(cons_seq); 
        free(cons_cov); 
        free(cons_l);
    }
    if (msa_l) 
    {
        for (i = 0; i < n_seqs; ++i) 
        {
            free(msa_seq[i]); 
            free(msa_seq);
        }
    }
    
    abpoa_free(ab); 
    abpoa_free_para(abpt); 

    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); 
    free(bseqs); 
    free(seq_lens);
    

    if (het)
    {
        if (debug_flag) {cerr << "het! " << endl; cerr << "h1 anno: " << endl;}
        outputConsensus(consensus[0], opt);
        if (debug_flag) cerr << "h2 anno: " << endl;
        outputConsensus(consensus[1], opt);
    }
    else
    {
        if (debug_flag) {cerr << "hom! " << endl; cerr << "h1&h2 anno: " << endl;}
        outputConsensus(consensus[0], opt);     
    }
    */

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
            if (naive_flag)
                naive_anno(consensus[i], motifs, Hap_seqs[i]->seq, Hap_seqs[i]->len);
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
            if (naive_flag)
                naive_anno(annos[i], motifs, reads[i]->seq, reads[i]->len);
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


// // skip "-" (45) for encoding
// static void encode (int motifs_size, vector<uint8_t> &annoNum, string &annoStr)
// {
//     string tmp_s;
//     uint8_t num;
//     for (size_t i = 0; i < annoNum.size(); ++i)
//     {
//         num = annoNum[i];
//         assert(num < motifs_size and motifs_size <= 255);

//         tmp_s.assign(1, (char) num);
//         annoStr.replace(i, 1, tmp_s);  
//         assert(!annoStr.empty());

//         // for outputing string to stdout
//         // if (num < 45)
//         // {
//         //  tmp_s.assign(1, (char) (num + 33));
//         //  annoStr.replace(i, 1, tmp_s);  
//         // }
//         // else if (num >= 45 and num <= 221)
//         // {
//         //  tmp_s.assign(1, (char) (num + 34));
//         //  annoStr.replace(i, 1, tmp_s);       
//         // }        
//     }
//     return;
// }

// void annoTostring_helper (int motifs_size, string &annoStr, vector<uint8_t> &annoNum)
// {
//     annoStr.clear();
//     annoStr.resize(annoNum.size(), '+');
//     encode(motifs_size, annoNum, annoStr);
//     return;
// }

// void VNTR::annoTostring (const OPTION &opt)
// {
//     if (skip) return;
//     annoStrs.resize(nreads);
//     for (int r = 0; r < nreads; ++r)
//     {
//         annoTostring_helper(motifs.size(), annoStrs[r], annos[r]);
//         assert(annoStrs[r].length() > 0);
//         // if (debug_flag) cerr << "string " << r << ":  " << annoStrs[r] << endl;
//     }
//     return;
// }

// size_t VNTR::getAnnoStrLen (int i)
// {
//     assert(i < nreads);
//     return annoStrs[i].length();
// }

// double findMedian(Order &a, int n)
// {
//     if (n % 2 == 0) 
//         return (double)(a[(n - 1) / 2] + a[n / 2]) / 2.0;
//     else 
//         return (double)a[n / 2];
// }
  
// pair<double, double> findQuantile (Order &a, int n)
// {
//     int Q2 = n / 4;
//     int Q3 = (n * 3) / 4;
//     return make_pair(a[Q2], a[Q3]);
// }

// int VNTR::cleanNoiseAnno(const OPTION &opt) 
// {
//     if (skip) return 0;

//     int i, j;
//     vector<vector<double>> edist(nreads);
//     for (i = 0; i < nreads; ++i) edist[i].resize(nreads, 0);

//     EdlibAlignResult result;
//     for (i = 0; i < nreads; ++i)
//     {
//         for (j = i + 1; j < nreads; ++j) 
//         {
//             result = edlibAlign(annoStrs[i].c_str(), getAnnoStrLen(i), annoStrs[j].c_str(), getAnnoStrLen(j), edlibDefaultAlignConfig());
//             if (result.status == EDLIB_STATUS_OK) 
//             {
//                 edist[i][j] = (double) result.editDistance; // [i][j]
//                 edist[j][i] = (double) result.editDistance; // [j][i]
//             }
//             else 
//             {
//                 printf("edlib goes wrong!\n");
//                 return 1;
//             }
//         }
//     }  

//     if (opt.filterNoisy)
//     {
//         vector<Order> edist_order(nreads);
//         for (i = 0; i < nreads; ++i)
//             edist_order[i].Update(&edist[i]);

//         vector<double> median(nreads, 0.0);
//         for (i = 0; i < nreads; ++i) 
//             median[i] = findMedian(edist_order[i], nreads);

//         Order median_order(&median);
//         median_order.Update(&median);
//         auto [q2, q3] = findQuantile(median_order, nreads);
//         auto IQR = q3 - q2;
//         auto outliner = q3 + opt.filterStrength * IQR;

//         vector<bool> remove(nreads, 0);
//         for (i = 0; i < nreads; ++i)
//         {
//             if (median[i] > outliner) remove[i] = 1;
//         }

//         ncleanreads = 0;
//         for (i = 0; i < nreads; ++i) 
//         {
//             if (remove[i] == 0)
//             {
//                 ncleanreads += 1;
//                 clean_reads.push_back(reads[i]);
//                 clean_annoStrs.push_back(annoStrs[i]);
//             }
//         }

//         clean_annos.resize(ncleanreads);
//         j = 0;
//         for (i = 0; i < nreads; ++i)
//         {
//             if (remove[i] == 0)
//             {
//                 clean_annos[j].resize(annos[i].size());
//                 copy (annos[i].begin(), annos[i].end(), clean_annos[j].begin());
//                 j += 1;
//             }
//         }

//         clean_edist.resize(ncleanreads);
//         int new_i = 0, new_j = 0;
//         for (i = 0; i < nreads; ++i)
//         {
//             if (remove[i] == 0)
//             {
//                 for (j = 0; j < nreads; ++j) 
//                 {
//                    if (remove[j] == 0)
//                      clean_edist[new_i].push_back(edist[i][j]); 
//                 }
//                 new_i += 1;
//             }
//         }    

//         for (i = 0; i < ncleanreads; ++i)
//             assert(clean_edist[i].size() == ncleanreads);

//         if (ncleanreads == 0) skip = true;
//         if (debug_flag) cerr << "ncleanreads: " << ncleanreads << endl;       
//     }
//     else
//     {
//         ncleanreads = nreads;
//         for (i = 0; i < nreads; ++i) 
//         {
//             clean_reads.push_back(reads[i]);
//             clean_annoStrs.push_back(annoStrs[i]);
//         }

//         clean_annos.resize(ncleanreads);
//         for (i = 0; i < nreads; ++i)
//             clean_annos[i] = annos[i];

//         clean_edist.resize(ncleanreads);
//         for (i = 0; i < nreads; ++i)
//             clean_edist[i] = edist[i]; 
//     }

//     return 0;    
// }

// int VNTR::hClust (vector<int> &gp1, vector<int> &gp2, double dists [])
// {
//     /*
//      compute the edit distance matrix for annotation strings 
//      input: vector<string> annoStrs
//      output: double edist[nreads * nreads] = {...}
//     */
//     gp1.clear(); gp2.clear();
//     if (nreads == 1)
//     {
//        gp1.push_back(0);
//        return 0;
//     }

//     int i, j;

//     /* do hierarchical clustering based on the edit ditance matrix */
//     alglib::clusterizerstate s;
//     alglib::ahcreport rep;
//     alglib::clusterizercreate(s);

//     alglib::real_2d_array d;
//     d.setcontent(nreads, nreads, dists);
//     alglib::clusterizersetdistances(s, d, true);
//     alglib::clusterizerrunahc(s, rep);

//     /* get the cluster results when nclusts == 2 */
//     alglib::integer_1d_array cidx;
//     alglib::integer_1d_array cz;    
//     clusterizergetkclusters(rep, 2, cidx, cz);

//     gp1.clear(); gp2.clear();
//     for (i = 0; i < nreads; ++i)
//     {
//         if (cidx[i] == 0)
//             gp1.push_back(i);
//         else
//             gp2.push_back(i);
//     }
//     return 0;
// }

// double avgSilhouetteCoeff (int nreads, double edist [], vector<int> &gp1, vector<int> &gp2)
// {
//     double a1 = 0.0;
//     double b1 = 10000000.0;
//     double sc1 = 0.0;

//     double a2 = 0.0;
//     double b2 = 10000000.0;
//     double sc2 = 0.0;

//     for (auto &i : gp1)
//     {
//         a1 = 0.0; b1 = 10000000.0;
//         for (auto &j : gp1) 
//         {
//             if (i == j) continue;
//             a1 += edist[i * nreads + j];
//         }

//         if (gp1.size() > 1)
//             a1 /= gp1.size() - 1;

//         for (auto &j : gp2)
//             b1 = min(edist[i * nreads + j], b1);

//         sc1 += (b1 - a1) / max(a1, b1);
//     }

//     for (auto &i : gp2)
//     {
//         a2 = 0.0; b2 = 10000000.0;
//         for (auto &j : gp2) 
//         {
//             if (i == j) continue;
//             a2 += edist[i * nreads + j];
//         }

//         if (gp2.size() > 1)
//             a2 /= gp2.size() - 1;

//         for (auto &j : gp1)
//             b2 = min(edist[i * nreads + j], b2);

//         sc2 += (b2 - a2) / max(a2, b2);
//     }

//     return (sc1 + sc2) / (gp1.size() + gp2.size());
// }

// /* based on vector<vector<uint8_t>> annos, get a consensus representation of the current locus 
//    input:   vector<vector<uint8_t>> annos
//    output:  vector<vector<uint8_t>> consensus
// */
// void VNTR::consensusMotifAnnoForOneVNTRByClustering (const OPTION &opt)
// {
//     if (skip) return;

//     int n_seqs = annos.size();
//     if (n_seqs == 0) return;

//     if (n_seqs == 1)
//     {
//         consensus.resize(1);
//         consensus[0] = annos[0];
//         return;
//     }

//     /* compute the MSA for each cluster when nclusts == 2 */
//     vector<int> gp1, gp2;
//     int motifs_size = motifs.size();

//     // change from edist to dists
//     int i, j;
//     double dists [nreads * nreads];

//     for (i = 0; i < nreads; ++i)
//     {
//         for (j = 0; j < nreads; ++j)
//            dists[i * nreads + j] = edist[i][j];
//     }        


//     hClust(gp1, gp2, dists); 

//     /* compare which is better: nclusts == 2 or nclusts == 1 */
//     double sc_2clusts = avgSilhouetteCoeff(nreads, dists, gp1, gp2);
//     if (sc_2clusts > 0.0 and (float) gp1.size() / nreads >= 0.4 and (float) gp2.size() / nreads >= 0.4)
//     {
//         het = true;
//         consensus.resize(2);

//         MSA * msa_1 = new MSA(gp1, annos_cp, motifs_size);
//         MSA * msa_2 = new MSA(gp2, annos_cp, motifs_size);
//         msa_1->MSA_anno_group (motifs.size(), gp1, annos_cp, consensus[0]);
//         msa_2->MSA_anno_group (motifs.size(), gp2, annos_cp, consensus[1]);
//         assert(consensus[0].size() > 0 and consensus[1].size() > 0);
//         delete msa_1;
//         delete msa_2;
//     }
//     else 
//     {
//         het = false;
//         consensus.resize(1);

//         MSA * msa = new MSA(gp1, annos_cp, motifs_size);

//         if ((float) gp1.size() / nreads >= (float) gp2.size() / nreads)
//         {

//             msa->MSA_anno_group (motifs.size(), gp1, annos_cp, consensus[0]);
//         }
//         else 
//         {
//             msa->MSA_anno_group (motifs.size(), gp2, annos_cp, consensus[0]);
//         }
//         assert(consensus[0].size() > 0);
//         delete msa;
//     }

//     if (het)
//     {
//         if (debug_flag) {cerr << "het! " << endl; cerr << "h1 anno: " << endl;}
//         outputConsensus(consensus[0], opt);
//         if (debug_flag) cerr << "h2 anno: " << endl;
//         outputConsensus(consensus[1], opt);
//     }
//     else
//     {
//         if (debug_flag) {cerr << "hom! " << endl; cerr << "h1&h2 anno: " << endl;}
//         outputConsensus(consensus[0], opt);     
//     }
//     return;
// }
