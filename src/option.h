#ifndef OPTION_H_
#define OPTION_H_

using namespace std;

class OPTION
{
public:
	bool fasterAnnoAlg;
	bool debug;

	OPTION ()
	{
		fasterAnnoAlg = true;
		debug = false;
	};

	OPTION (bool FasterAnnoAlg, bool Debug) : fasterAnnoAlg(FasterAnnoAlg), debug(Debug) {};
	~OPTION () {};
};

#endif
