#ifndef READ_H_
#define READ_H_
#include <string>
#include <vector>

using namespace std;

/* 
class READ constains
@qname: the read ID
@chr: the chromosome the read is mapped to
@len: the read length
@seq: the read sequence
*/

class SNV {
public:
  int pos;
  int nuc;
};

class READ
{
public: 
	char * qname;
	uint16_t l_qname;
	char * chr;
	uint32_t len;
	char * seq;
	bool rev;
    int haplotype;
    string upstream, downstream;
  	vector<SNV> snvs;
  	uint16_t flag=0;
	READ () 
	{
		qname = NULL;
		l_qname = 0;
		chr = NULL;
		seq = NULL;
		rev = false;
		haplotype = 0;
	};
	READ (char * Qname, char * Chr, uint32_t Len) : qname(Qname), chr(Chr), len(Len) { seq = NULL;};
	~READ () 
	{
		free(seq);
		free(qname);
	};
};

class less_than_key
{
public:
    inline bool operator() (READ * read1, READ * read2)
    {
        return (read1->len < read2->len);
    }
};

#endif
