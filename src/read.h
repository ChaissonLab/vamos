#ifndef READ_H_
#define READ_H_
#include <string>
#include <vector>

using namespace std;



/**
 * @brief class for one SNV
 * 
 */
class SNV {
public:
    int pos;        // position of this SNV
    int nuc;        // nucleotide base of this SNV
};




/**
 * @brief class for one aligned read/contig
 */
class READ
{
public: 
    char * qname;                   // query name of this read/contig
    uint16_t l_qname;               // 
    char * chr;                     // aligned reference chromosome of this read/contig
    uint32_t len;                   // length of this read/contig
    char * seq;                     // sequence of this read/contig
    bool rev;                       // if this read/contig is reverse aligned
    int haplotype;                  // haplotype of this read/contig, 0 for homo, 1/2 for heter hap1/hap2
    string upstream, downstream;    // 
    vector<SNV> snvs;               // vector of all SNVs of this read/contig
    uint16_t flag=0;                // sam flag of this read/contig


    /// @brief Construct a new READ object
    READ()
    {
        qname = NULL;
        l_qname = 0;
        chr = NULL;
        seq = NULL;
        rev = false;
        haplotype = 0;
    };
    READ(char * Qname, char * Chr, uint32_t Len) : qname(Qname), chr(Chr), len(Len) { seq = NULL;};


    /// @brief Destroy the READ object
    ~READ() 
    {
        free(seq);
        free(qname);
    };
};


#endif
