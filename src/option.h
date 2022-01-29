#ifndef OPTION_H_
#define OPTION_H_

// using namespace std;

class OPTION
{
public:
	int nproc;
	bool filterNoisy;
	double filterStrength;

	OPTION ()
	{
		nproc = 1;
		filterNoisy = false;
		filterStrength = 0.0;
	};

	~OPTION () {};
};

#endif
