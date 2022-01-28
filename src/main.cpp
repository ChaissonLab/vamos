#include <stdio.h>     
#include <stdlib.h>   
#include <getopt.h>
#include <string>
#include <unistd.h>
#include "vntr.h"
#include "vcf.h"
#include "option.h"
#include <sys/time.h>
#include "threads.h"
#include "io.h"
#include <mutex>   

using namespace std;

struct timeval start_time, stop_time, elapsed_time;

/* vamos -in read.bam -vntr vntrs.bed -motif motifs.csv -o out.vcf */

void printUsage(IO &io) 
{
	printf("Usage: vamos [-h] [-i in.bam] [-v vntrs.bed] [-m motifs.csv] [-o output.vcf] [-s sample_name] [-t threads] [-n] [-d] [-c] [-f]\n");
	printf("Version: %s\n", io.version);
	printf("Options:\n");
	printf("       -i  FILE      input alignment file (bam format), bam file needs to be indexed \n");
	printf("       -v  FILE      the tab-delimited coordinate of each VNTR locus - `chrom\tstart\tend`, each row represents a VNTR locus\n");
	printf("       -m  FILE      the comma-delimited motif sequences list for each VNTR locus, each row represents a VNTR locus\n");
	printf("       -o  FILE      output vcf file\n");
	printf("       -s  CHAR      the sample name\n");
	printf("       -t  INT       number of threads\n");
	printf("       -n            specify the naive version of code to do the annotation, default is faster implementation\n");
	printf("       -d            print out debug information\n");
	printf("       -c            use hierarchical clustering to judge if a VNTR locus is het or hom\n");
	printf("       -f            filter noisy read annotations\n");
	printf("       -h            print out help message\n");
} 

void ProcVNTR (int s, VNTR * it, const OPTION &opt) 
{
	if (it->nreads == 0 or it->motifs.size() > 255) 
	{
		it->skip = true;
		if (opt.debug) cerr << "skip one vntr due to > 255 motifs" << endl;
		return;
	}

	if (opt.debug) cerr << "start to do the annotation: " << s << endl;

	it->motifAnnoForOneVNTR(opt); 
	it->annoTostring(opt);
	it->cleanNoiseAnno(opt);
	if (opt.hc)
		it->concensusMotifAnnoForOneVNTR(opt);
	else
		it->concensusMotifAnnoForOneVNTRUsingABpoa(opt);
	it->clearRead();
	return;
}

void *ProcVNTRs (void *procInfoValue)
{
	ProcInfo *procInfo = (ProcInfo *)procInfoValue;
	cerr << "start thread: " << procInfo->thread << endl;
	int i, s;
	int sz = (procInfo->vntrs)->size();

	// read bam 
	(procInfo->io)->readSeqFromBam ((*(procInfo->vntrs)), (procInfo->opt)->nproc, procInfo->thread, sz);

	for (i = procInfo->thread, s = 0; i < sz; i += (procInfo->opt)->nproc, s += 1)
	{
		if (procInfo->opt->debug) cerr << "processing vntr: " << i << endl;
		// procInfo->numOfProcessed += ProcVNTR (s, (*(procInfo->vntrs))[i], *(procInfo->opt));
		ProcVNTR (s, (*(procInfo->vntrs))[i], *(procInfo->opt));
	}

	procInfo->mtx->lock();

	cerr << "outputing vcf" << endl;
	(procInfo->io)->writeVCFBody(*(procInfo->out), (*(procInfo->vntrs)), procInfo->thread, (procInfo->opt)->nproc);	
	// (procInfo->io->vcfWriter).writeNullAnno((*(procInfo->vntrs)), *(procInfo->out_nullAnno), procInfo->thread, (procInfo->opt)->nproc);	

	procInfo->mtx->unlock();
	cerr << "finish thread: " << procInfo->thread << endl;
	pthread_exit(NULL);     /* Thread exits (dies) */	
}

int main (int argc, char **argv)
{
	int c;
	IO io;
	OPTION opt;

	const struct option long_options[] =
	{
		/* These options don’t set a flag. We distinguish them by their indices. */
		{"input",         required_argument,       0, 'i'},
		{"vntr",          required_argument,       0, 'v'},
		{"motif",         required_argument,       0, 'm'},
		{"output",        required_argument,       0, 'o'},
		{"sampleName",    required_argument,       0, 's'},
		{"numThreads",    required_argument,       0, 't'},
		{"naiveAnnoAlg",  no_argument,             0, 'n'},
		{"hclust",        no_argument,             0, 'c'},
		{"filterNoisy",   required_argument,       0, 'f'},
		{"help",          no_argument,             0, 'h'},
		{NULL, 0, 0, '\0'}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;

	while ((c = getopt_long (argc, argv, "i:v:m:o:s:t:f:hndc", long_options, &option_index)) != -1)
	{
		switch (c)
		{
		case 'i':
			printf ("option -input with `%s'\n", optarg);
			io.input_bam = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.input_bam, optarg);
			break;

		case 'v':
			printf ("option -vntr with `%s'\n", optarg);
			io.vntr_bed = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.vntr_bed, optarg);
			break;

		case 'm':
			printf ("option -motif with `%s'\n", optarg);
			io.motif_csv = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.motif_csv, optarg);
			break;

		case 'o':
			printf ("option -output with `%s'\n", optarg);
			io.out_vcf = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.out_vcf, optarg);
			break;

		case 's':
			printf ("option -sampleName with `%s'\n", optarg);
			io.sampleName = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.sampleName, optarg);
			break;

		case 't':
			printf ("option -numThreads with `%s'\n", optarg);
			opt.nproc = atoi(optarg);
			break;

		case 'n':
			printf ("option -naiveAnnoAlg");
			opt.fasterAnnoAlg = false;
			break;

		case 'd':
			printf ("option -debug\n");
			opt.debug = true;
			break;

		case 'c':
			printf ("option -hclust\n");
			opt.hc = true;
			break;

		case 'f':
			printf ("option -filterNoisy\n");
			opt.filterStrength = stod(optarg);
			opt.filterNoisy = true;
			break;

		case 'h':
			printUsage(io);
			exit(EXIT_SUCCESS);

		case '?':
			printf("Unknown option: %c\n", optopt);
			exit(EXIT_FAILURE);

		case ':':
			cerr << "[ERROR] missing option argument" << endl;
			exit(EXIT_FAILURE);

		default:
			printUsage(io);
			exit(EXIT_FAILURE);
		}		
	}
	
	/* Check mandatory parameters */
	bool missingArg = false;
	if (io.input_bam == NULL) 
	{
		printf("-i is mandatory!\n");
		missingArg = true;
	}
	if (io.vntr_bed == NULL)
	{
		printf("-v is mandatory!\n");
		missingArg = true;
	}
	if (io.motif_csv == NULL)
	{
		printf("-m is mandatory!\n");
		missingArg = true;
	}
	if (io.out_vcf == NULL)
	{
		printf("-o is mandatory!\n");
		missingArg = true;
	}
	if (io.sampleName == NULL)
	{
		printf("-s is mandatory!\n");
		missingArg = true;
	}

	if (missingArg)
	{
		printUsage(io);
		exit(EXIT_FAILURE);
	}

	vector<VNTR *> vntrs;

	gettimeofday(&start_time,NULL);

	/* read VNTR bed file */
	io.readVNTRFromBed(vntrs);

	cerr << "finish reading vntrs.bed" << endl;

	/* read motif csv file */
	io.readMotifsFromCsv(vntrs);

	cerr << "finish reading motifs.csv" << endl;

	/* process each VNTR */
	// io.readSeqFromBam(vntrs); 

	/* set up out stream and write VCF header */
    ofstream out(io.out_vcf);
    if (out.fail()) 
    {
        cerr << "ERROR: Unable to open file " << io.out_vcf << endl;
        exit(EXIT_FAILURE);
    } 	
	io.writeVCFHeader(out);
	cerr << "finishing reading!" << endl;

	gettimeofday(&stop_time,NULL);
	timersub(&stop_time, &start_time, &elapsed_time); 
	printf("reading time was %f sec\n",
	elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

	/* debug code */
	// ofstream out_nullAnno("/project/mchaisso_100/cmb-16/jingwenr/trfCall/vamos/src/test_pipeline/sub_50000/nullAnno.bed");

	/* Create threads */
	int i;
	if (opt.nproc > 1)
	{
		pthread_t *tid = new pthread_t[opt.nproc];
		mutex mtx; 
		int numOfProcessed = 0;
		vector<ProcInfo> procInfo(opt.nproc);		
		for (i = 0; i < opt.nproc; i++){ 
			procInfo[i].vntrs = &vntrs;
			procInfo[i].thread = i;
			procInfo[i].opt = &opt;
			procInfo[i].io = &io;
			procInfo[i].mtx = &mtx;
			procInfo[i].numOfProcessed = &numOfProcessed;
			procInfo[i].out = &out;
			// procInfo[i].out_nullAnno = &out_nullAnno;

			if (pthread_create(&tid[i], NULL, ProcVNTRs, (void *) &procInfo[i]) )
			{
				cerr << "ERROR: Cannot create thread" << endl;
				exit(EXIT_FAILURE);
			}
		}

		for (i = 0; i < opt.nproc; i++) {
			pthread_join(tid[i], NULL);
		}
		
		delete[] tid;
		// for (i = 1; i < opt.nproc; i++) 
		// 	procInfo[0].timing.Add(procnfo[i].timing);
		
		// if (opt.timing != "") 
		// 	procInfo[0].timing.Summarize(opt.timing);
		
	}
	else 
	{
		io.readSeqFromBam (vntrs, 1, 0, vntrs.size());
		int s = 0;
		for (auto &it: vntrs) 
		{
			ProcVNTR (s, it, opt);
			s += 1;
		}	
		cerr << "outputing vcf" << endl;
		io.writeVCFBody(out, vntrs, -1, 1);
	}

	for (size_t i = 0; i < vntrs.size(); ++i) 
		delete vntrs[i];
	
	out.close();
	// out_nullAnno.close();
	exit(EXIT_SUCCESS);
}

