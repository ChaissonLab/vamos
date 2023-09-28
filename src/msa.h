#ifndef MSA_H_
#define MSA_H_

#include "read.h"
#include "abpoa.h"
#include "option.h"


class MSA
{
public:
    uint8_t ** cons_seq;                    // sequence of the best msa consensus sequence
    int * cons_l;                           // length of the best msa consensus sequence
    int ** cons_cov;                        // 
    int cons_n;                             // number of all consensus sequences
    uint8_t ** msa_seq;                     // 
    int msa_l;                              // 
    int n_seqs, n_empty_seqs;                             // total number of read sequences, total number of deletion sequences.
    int * seq_lens;                         // read sequences length
    uint8_t ** bseqs;                       // read sequences in integer 0-4 (AaCcGgTtNn ==> 0,1,2,3,4)
    abpoa_t * ab;                           // struct for abpoa parameters
    abpoa_para_t * abpt;                    // struct for abpoa parameters
    unsigned char ab_nt4_table[256] = {     // table to convert nt characters to integer 0-4 (AaCcGgTtNn ==> 0,1,2,3,4)
        0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };


    MSA() 
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
	n_empty_seqs = 0;
        // initialize variables
    };

  ~MSA() {
    cleanMSA();
  }
    /**
     * @brief Construct a new MSA object for msa of anno sequences - initialize all anno seqs as vector of motif indices
     * @param annos annotations (by motif index vector) of all reads for this VNTR
     * @param motifs_size number of motifs for this VNTR
     */
    MSA(const vector<vector<uint8_t> > &annos, int motifs_size)
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

        // sort annos as annos_cp
        vector<vector<uint8_t> > annos_cp(annos);
        sort(annos_cp.begin(), annos_cp.end(), [](vector<uint8_t> & a, vector<uint8_t> & b)
            {return a.size() > b.size();});

        int i, j;
        for (i = 0; i < n_seqs; ++i) {

            seq_lens[i] = annos_cp[i].size();
            bseqs[i] = new uint8_t[seq_lens[i] + 1];
            // bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * (seq_lens[i] + 1));

            for (j = 0; j < seq_lens[i]; ++j)
            {
                assert(annos_cp[i][j] < motifs_size);
                bseqs[i][j] = annos_cp[i][j] + 1; // +1 is cancelled when the cons_seq is converted back
            }
            bseqs[i][seq_lens[i]] = '\0';
        }
    };


    /**
     * @brief Construct a new MSA object for msa of raw sequences - initialize all raw nt seqs as 0-4 nt encoding
     * @param gp vector of indices of all reads in a haplotype group
     * @param reads vector of all reads
     */
    MSA(const vector<int> &gp, const vector<READ *> &reads)
    {
        cons_seq = NULL;
        cons_l = NULL;
        cons_cov = NULL;
        cons_n = 0;
        msa_seq = NULL;
        msa_l = 0;
        ab = abpoa_init();
        abpt = abpoa_init_para();
	n_seqs = 0;
	n_empty_seqs = 0;
	int numNonZero = 0;
	for (int i = 0; i < gp.size(); i++) {
	  if (reads[gp[i]]->seq.size() > 0) {
	    n_seqs++;
	  }
	  else {
	    n_empty_seqs++;
	  }
	}
	
	if (n_seqs > 0) {
	    seq_lens = new int[n_seqs+1]; // n_seqs+1, one extra read as padding?
	    bseqs = new uint8_t*[n_seqs];
	}
	else {
	  return;
	}
        int i, j;
        seq_lens[n_seqs] = 0; // n_seqs+1, one extra read as padding?
	int cur=0;
        for (i = 0; i < gp.size(); ++i)
        {
            assert(gp[i] < reads.size());
	    if (reads[gp[i]]->seq.size() > 0) {
	      seq_lens[cur] = reads[gp[i]]->seq.size();
	      bseqs[cur] = new uint8_t[seq_lens[cur] + 1];
	      memcpy(bseqs[cur], reads[gp[i]]->seq.c_str(), seq_lens[cur]);
	      
	      for (auto j=0; j < seq_lens[cur]; j++)
		{
		  bseqs[cur][j] = ab_nt4_table[bseqs[cur][j]];
		}
	      bseqs[cur][seq_lens[cur]] = '\0'; // last position of each read as '\0' (NULL character)
	      cur++;	      
	    }

	    /*
	    else {
	      seq_lens[i] = 0;
	      bseqs[i] = new uint8_t[1];
	      bseqs[i][0] = '\0';
	      }*/
        }
    };


    /**
     * @brief Run the abpoa msa
     */
    void runConsensus();


    /**
     * @brief Convert the consensus (by motif indices + 1) back to the original motif indices string
     * @param consensus converted final MSA consensus (by motif indices)
     * @param motifs_size number of motifs
     */
    void extractConsensusAnno(vector<uint8_t> &consensus, int motifs_size);


    /**
     * @brief Convert the consensus (by 0-4 nt encoding) back to the original nt sequence
     * @param Hap_seq converted final MSA consensus (by nt AaCcGgTtNn)
     */
    void extractConsensusRawSeq(READ * Hap_seq);


    /**
     * @brief Generate consensus sequence for a group of motif anno sequences through abpoa
     * @param motifs_size number of motifs
     * @param annos annotations (by motif index vector) of all reads for this VNTR
     * @param consensus converted final MSA consensus (by motif indices)
     */
    void MSA_anno_group(int motifs_size, const vector<vector<uint8_t>> &annos, vector<uint8_t> &consensus);


    /**
     * @brief Generate consensus sequence for a group of read sequences through abpoa
     * @param gp indeces of all input reads for MSA
     * @param reads all input reads
     * @param Hap_seq converted final MSA consensus (by nt AaCcGgTtNn)
     */
    void MSA_seq_group(const vector<int> &gp, const vector<READ *> &reads, READ * Hap_seq);


    /**
     * @brief Clear memory
     */
    void cleanMSA();
};


#endif
