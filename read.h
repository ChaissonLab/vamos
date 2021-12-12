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
	char * chr;
	uint32_t len;
	char * seq;
	READ () {};
	READ (char * Qname, char * Chr, uint32_t Len) : qname(Qname), chr(Chr), len(Len) {};
	~READ () 
	{
		free(qname);
		free(chr);
		free(seq);
	};
};

#endif
