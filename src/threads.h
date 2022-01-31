#ifndef THREADS_H_
#define THREADS_H_

#include <stdio.h>     
#include <stdlib.h>   
#include "io.h"
#include "vntr.h"
#include "option.h"
#include <mutex>   
#include <sys/time.h>



class ProcInfo
{
public:
	vector<VNTR *> * vntrs;
	int thread;
	OPTION * opt;
	IO * io;
	ofstream * out;
	ofstream * out_nullAnno;
	mutex *mtx; 
	int * numOfProcessed;
	struct timeval start_time, stop_time, elapsed_time;

	ProcInfo () {};
	~ProcInfo () {};
};

#endif
