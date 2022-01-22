#ifndef READ_H_
#define READ_H_
#include <string>
using namespace std;
/* 
class READ constains
@qname: the read ID
@chr: the chromosome the read is mapped to
@len: the read length
@seq: the read sequence
*/
class READ
{
public: 
	char * qname;
	uint16_t l_qname;
	char * chr;
	uint32_t len;
	char * seq;
	bool rev;
	READ () 
	{
		qname = NULL;
		l_qname = 0;
		chr = NULL;
		seq = NULL;
		rev = false;
	};
	READ (char * Qname, char * Chr, uint32_t Len) : qname(Qname), chr(Chr), len(Len) { seq = NULL;};
	~READ () 
	{
		free(seq);
		free(qname);
	};
};

#endif
