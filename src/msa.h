#ifndef MSA_H_
#define MSA_H_

#include "read.h"
#include "abpoa.h"
#include "option.h"


class MSA
{
public: 
    uint8_t **cons_seq; 
    int *cons_l;
    int **cons_cov;
    int cons_n;
    uint8_t **msa_seq;
    int msa_l;
    int n_seqs;
    int *seq_lens;
    uint8_t **bseqs;
    abpoa_t *ab;
    abpoa_para_t *abpt;

    MSA () 
    {
        seq_lens = NULL; 
        bseqs = NULL; 
        cons_seq = NULL; 
        cons_l = NULL; 
        cons_cov = NULL; 
        cons_n = 0; 
        msa_seq = NULL;
        msa_l = 0; 
        n_seqs = 0;
        // initialize variables
        ab = abpoa_init();
        abpt = abpoa_init_para();
    };

    MSA (const vector<vector<uint8_t>> &annos, int motifs_size)
    {
        cons_seq = NULL; 
        cons_l = NULL; 
        cons_cov = NULL; 
        cons_n = 0; 
        msa_seq = NULL;
        msa_l = 0; 
        ab = abpoa_init();
        abpt = abpoa_init_para();

        n_seqs = annos.size();
        seq_lens = new int[n_seqs];
        bseqs = new uint8_t*[n_seqs];

        // sort annos
        vector<vector<uint8_t>> annos_cp(annos);
        sort(annos_cp.begin(), annos_cp.end(), [](vector<uint8_t> & a, vector<uint8_t> & b){return a.size() > b.size();});

        int i, j;
        for (i = 0; i < n_seqs; ++i) {
            seq_lens[i] = annos_cp[i].size();
            bseqs[i] = new uint8_t[seq_lens[i] + 1];
            // bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * (seq_lens[i] + 1));
            for (j = 0; j < seq_lens[i]; ++j)
            {
                assert(annos_cp[i][j] < motifs_size);
                bseqs[i][j] = annos_cp[i][j] + 1;            
            }
            bseqs[i][seq_lens[i]] = '\0';
        }
    };

    MSA (const vector<int> &gp, const vector<READ *> &reads)
    {
        cons_seq = NULL; 
        cons_l = NULL; 
        cons_cov = NULL; 
        cons_n = 0; 
        msa_seq = NULL;
        msa_l = 0; 
        ab = abpoa_init();
        abpt = abpoa_init_para();

        n_seqs = gp.size();
        seq_lens = new int[n_seqs];
        bseqs = new uint8_t*[n_seqs];
        int i, j;
        for (i = 0; i < n_seqs; ++i) {
            seq_lens[i] = reads[gp[i]]->len;
            bseqs[i] = new uint8_t[seq_lens[i] + 1];
            
            for (j = 0; j < seq_lens[i]; ++j)
            {
                bseqs[i][j] = reads[gp[i]]->seq[j];            
            }
            bseqs[i][seq_lens[i]] = '\0';

        }
    };

    // void encodeSeqsFromAnno (const vector<int> &gp, const vector<vector<uint8_t>> &annos, int _n_seqs, int motifs_size);

    // void encodeSeqsFromRaw (const vector<int> &gp, const vector<READ *> &reads);

    void runConsensus ();

    void extractConsensusAnno (vector<uint8_t> &consensus, int motifs_size);

    void extractConsensusRawSeq (READ * Hap_seq);

    void MSA_anno_group (int motifs_size, const vector<vector<uint8_t>> &annos, vector<uint8_t> &consensus);

    void MSA_seq_group (const vector<int> &gp, const vector<READ *> &reads, READ * Hap_seq);

    void cleanMSA ();
};



#endif
