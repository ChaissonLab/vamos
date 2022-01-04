#ifndef OPTION_H_
#define OPTION_H_

using namespace std;

class OPTION
{
public:
	bool fasterAnnoAlg;
	bool debug;
	int nproc;

	OPTION ()
	{
		fasterAnnoAlg = true;
		debug = false;
		nproc = 1;
	};

	~OPTION () {};
};

#endif
