#ifndef THREADS_H_
#define THREADS_H_

#include <stdio.h>     
#include <stdlib.h>   
#include "io.h"
#include "vntr.h"
#include "option.h"
#include <mutex>
#include <map>
#include <sys/time.h>



class ProcInfo
{
public:
	vector<VNTR *> * vntrs;
        map<string, vector<int > > *vntrMap;
        vector<bool> *procChrom;
        map<string, vector<char > > *processed;
        map<string, vector<int > > *bucketEndPos;
        map<string, vector<int > > *bucketStartIndex;
        map<string, vector<int > > *bucketEndIndex;  
	int thread;
	OPTION * opt;
	IO * io;
	ofstream * out;
	ofstream * out_nullAnno;
	mutex *mtx; 
	int * numOfProcessed;
	struct timeval start_time, stop_time, elapsed_time;
        vector< int > *mismatchCI;
        vector<Pileup > *pileups;
        bool refine;
        int oneOffset;
        int forcePhase;
	ProcInfo () {};
	~ProcInfo () {};
};

#endif
