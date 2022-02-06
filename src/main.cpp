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
#include <unistd.h>
#include <iomanip>
#include <errno.h>

// using namespace std;
int naive_flag = false;
int debug_flag = false;
int hclust_flag = false;
int seqan_flag = false;
int output_read_anno_flag = false;
int conseq_flag = false;
int conseq_anno_flag = false;
int raw_anno_flag = false;

struct timeval pre_start_time, pre_stop_time, pre_elapsed_time;
struct timeval single_start_time, single_stop_time, single_elapsed_time;

/* vamos -in read.bam -vntr vntrs.bed -motif motifs.csv -o out.vcf */

void process_mem_usage(double &vm_usage, double &resident_set)
{
    vm_usage = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    string ignore;
    ifstream ifs("/proc/self/stat", ios_base::in);
    ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
        >> ignore >> ignore >> vsize >> rss;

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = (rss * page_size_kb) / (1024.0 * 1024.0);
}

void ProcVNTR (int s, VNTR * it, const OPTION &opt) 
{
	if (it->nreads == 0 or it->motifs.size() > 255) 
	{
		it->skip = true;
		if (debug_flag and it->nreads == 0) cerr << "no reads" << endl;
		if (debug_flag and it->motifs.size() > 255) cerr << "skip one vntr due to > 255 motifs" << endl;

		return;
	}

	if (debug_flag) cerr << "start to do the annotation: " << s << endl;

	it->motifAnnoForOneVNTR(opt); 
	it->annoTostring(opt);
	it->cleanNoiseAnno(opt);
	if (hclust_flag)
		it->concensusMotifAnnoForOneVNTR(opt);
	else if (seqan_flag)
		it->concensusMotifAnnoForOneVNTRBySeqan(opt);
	else
		it->concensusMotifAnnoForOneVNTRByABpoa(opt);
	return;
}

void *ProcVNTRs (void *procInfoValue)
{
	ProcInfo *procInfo = (ProcInfo *)procInfoValue;
	gettimeofday(&(procInfo->start_time), NULL);
	cerr << "start thread: " << procInfo->thread << endl;
	int i, s;
	int sz = (procInfo->vntrs)->size();

	// read bam 
	(procInfo->io)->readSeqFromBam ((*(procInfo->vntrs)), (procInfo->opt)->nproc, procInfo->thread, sz);

	for (i = procInfo->thread, s = 0; i < sz; i += (procInfo->opt)->nproc, s += 1)
	{
		if (debug_flag) cerr << "processing vntr: " << i << endl;
		// procInfo->numOfProcessed += ProcVNTR (s, (*(procInfo->vntrs))[i], *(procInfo->opt));
		ProcVNTR (s, (*(procInfo->vntrs))[i], *(procInfo->opt));
	}

	procInfo->mtx->lock();
	cerr << "outputing vcf" << endl;
	(procInfo->io)->writeVCFBody(*(procInfo->out), (*(procInfo->vntrs)), procInfo->thread, (procInfo->opt)->nproc);	

	procInfo->mtx->unlock();
	cerr << "finish thread: " << procInfo->thread << endl;
	gettimeofday(&(procInfo->stop_time), NULL);
	timersub(&(procInfo->stop_time), &(procInfo->start_time), &(procInfo->elapsed_time)); 
	pthread_exit(NULL);     /* Thread exits (dies) */	
}

void printUsage(IO &io) 
{
	printf("Usage: vamos [subcommand] [options] [-i in.bam] [-v vntrs.bed] [-m motifs.csv] [-o output.vcf] [-s sample_name] [-x subsequence.fa]\n");
	printf("Version: %s\n", io.version);
	printf("subcommand:\n");
	printf("vamos --conseq      [-i in.bam] [-v vntrs.bed] [-o output.fa] ");
	printf("vamos --conseq_anno [-i in.fa]  [-v vntrs.bed] [-m motifs.csv] [-o output.vcf] [-s sample_name] (ONLY FOR SINGLE LOCUS!!)");
	printf("vamos --raw_anno    [-i in.bam] [-v vntrs.bed] [-m motifs.csv] [-o output.vcf] [-s sample_name]");
	printf("   Input:\n");
	printf("       -i   FILE         input alignment/read file (bam/fa), bam file needs to be indexed \n");
	printf("       -v   FILE         the tab-delimited coordinate of each VNTR locus - `chrom\tstart\tend`, each row represents a VNTR locus\n");
	printf("       -m   FILE         the comma-delimited motif sequences list for each VNTR locus, each row represents a VNTR locus\n");
	printf("       -o   FILE         output vcf/fa file\n");
	printf("       -s   CHAR         the sample name\n");
	printf("   Dynamic Programming:\n");
	printf("       -pi  DOUBLE       penalty of indel in dynamic programming (double) DEFAULT: 1.0\n");
	printf("       -pm  DOUBLE       penalty of mismatch in dynamic programming (double) DEFAULT: 1.0\n");
	printf("       --naive           specify the naive version of code to do the annotation, DEFAULT: faster implementation\n");
	printf("   Aggregate Annotation:\n");
	printf("       -f   DOUBLE       filter noisy read annotations, DEFAULT: 0.0 (no filter)\n");
	printf("       --clust           use hierarchical clustering to judge if a VNTR locus is het or hom\n");
	printf("       --seqan           use seqan lib to do MSA (haploid only), DEFAULT: abPoa\n");
	printf("       --readanno        output read annotation in VCF and output vntr sequences to stdin\n");
	printf("   General Setting:\n");
	printf("       -t   INT          number of threads, DEFAULT: 1\n");
	printf("       --debug           print out debug information\n");
	printf("       -h                print out help message\n");
} 

int main (int argc, char **argv)
{
	int c;
	IO io;
	OPTION opt;

	const struct option long_options[] =
	{
		/* These options set a flag. */
		{"naive",         no_argument,             &naive_flag,                    1},
		{"debug",         no_argument,             &debug_flag,                    1},
		{"clust",         no_argument,             &hclust_flag,                   1},
		{"seqan",         no_argument,             &seqan_flag,                    1},
		{"readanno",      no_argument,             &output_read_anno_flag,         1},
		{"raw_anno",      no_argument,             &raw_anno_flag,                 1},
		{"conseq_anno",   no_argument,             &conseq_anno_flag,              1},
		{"conseq",        no_argument,             &conseq_flag,                   1},

		/* These options don’t set a flag. We distinguish them by their indices. */
		{"input",           required_argument,       0, 'i'},
		{"vntr",            required_argument,       0, 'v'},
		{"motif",           required_argument,       0, 'm'},
		{"output",          required_argument,       0, 'o'},
		{"out_fa",          required_argument,       0, 'x'},
		{"sampleName",      required_argument,       0, 's'},
		{"numThreads",      required_argument,       0, 't'},
		{"filterNoisy",     required_argument,       0, 'f'},
		{"penlaty_indel",   required_argument,       0, 'd'},
		{"penlaty_mismatch",required_argument,       0, 'c'},
		{NULL, 0, 0, '\0'}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;
	while ((c = getopt_long (argc, argv, "i:v:m:o:s:t:f:d:c:x:h", long_options, &option_index)) != -1)
	{
		switch (c)
		{

        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0) break;
          fprintf (stderr, "option %s", long_options[option_index].name);
          if (optarg) fprintf (stderr, " with arg %s", optarg);
          break;

		case 'i':
			fprintf (stderr, "option -input with `%s'\n", optarg);
			io.input_bam = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.input_bam, optarg);
			break;

		case 'v':
			fprintf (stderr, "option -vntr with `%s'\n", optarg);
			io.vntr_bed = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.vntr_bed, optarg);
			break;

		case 'm':
			fprintf (stderr, "option -motif with `%s'\n", optarg);
			io.motif_csv = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.motif_csv, optarg);
			break;

		case 'o':
			fprintf (stderr, "option -output with `%s'\n", optarg);
			io.out_vcf = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.out_vcf, optarg);
			break;

		case 's':
			fprintf (stderr, "option -sampleName with `%s'\n", optarg);
			io.sampleName = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.sampleName, optarg);
			break;

		case 't':
			fprintf (stderr, "option -numThreads with `%s'\n", optarg);
			opt.nproc = atoi(optarg);
			break;

		case 'd':
			opt.penalty_indel = stod(optarg);
			fprintf (stderr, "option -penlaty_indel with `%f'\n", opt.penalty_indel);
			break;

		case 'c':
			opt.penalty_mismatch = stod(optarg);
			fprintf (stderr, "option -penlaty_mismatch with `%f'\n", opt.penalty_mismatch);
			break;

		case 'f':
			fprintf (stderr, "option -filterNoisy\n");
			opt.filterStrength = stod(optarg);
			opt.filterNoisy = true;
			break;

		case 'h':
			printUsage(io);
			exit(EXIT_SUCCESS);

		case '?':
			fprintf(stderr, "Unknown option: %c\n", optopt);
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
		fprintf(stderr, "-i is mandatory!\n");
		missingArg = true;
	}
	if (io.vntr_bed == NULL)
	{
		fprintf(stderr, "-v is mandatory!\n");
		missingArg = true;
	}
	if (io.motif_csv == NULL and (conseq_anno_flag or raw_anno_flag))
	{
		fprintf(stderr, "-m is mandatory!\n");
		missingArg = true;
	}
	if (io.sampleName == NULL)
	{
		fprintf(stderr, "-s is mandatory!\n");
		missingArg = true;
	}
	if (io.out_vcf == NULL)
	{
		fprintf(stderr, "-o is mandatory!\n");
		missingArg = true;
	}

	if (missingArg)
	{
		printUsage(io);
		exit(EXIT_FAILURE);
	}

  	/* Instead of reporting ‘--verbose’
     and ‘--brief’ as they are encountered,
     we report the final status resulting from them. */
  	if (naive_flag) fprintf(stderr, "naive_flag is set");
  	if (debug_flag) fprintf(stderr, "debug_flag is set");
   	if (hclust_flag) fprintf(stderr, "hclust_flag is set");
  	if (seqan_flag) fprintf(stderr, "seqan_flag is set");
   	if (output_read_anno_flag) fprintf(stderr, "output_read_anno_flag is set");
  	if (conseq_flag) fprintf(stderr, "conseq_flag is set");
  	if (conseq_anno_flag) fprintf(stderr, "conseq_anno_flag is set");
  	if (raw_anno_flag) fprintf(stderr, "raw_anno_flag is set");

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		fprintf (stderr, "non-option ARGV-elements: ");
		while (optind < argc) fprintf (stderr, "%s ", argv[optind++]);
		putchar ('\n');
	}

	vector<VNTR *> vntrs;

	gettimeofday(&pre_start_time, NULL);

	/* read VNTR bed file */
	io.readVNTRFromBed(vntrs);
	cerr << "finish reading vntrs.bed" << endl;

	/* read motif csv file */
	if (!conseq_flag)
	{
		io.readMotifsFromCsv(vntrs);
		cerr << "finish reading motifs.csv" << endl;		
	}


	ofstream out;
	/* set up out stream and write fa/VCF header */
	out.open(io.out_vcf, ofstream::out);
	if (out.fail()) 
	{
	    cerr << "ERROR: Unable to open file " << io.out_vcf << endl;
	    exit(EXIT_FAILURE);
	} 	
	if (!conseq_flag)
		io.writeVCFHeader(out);

	gettimeofday(&pre_stop_time, NULL);
	timersub(&pre_stop_time, &pre_start_time, &pre_elapsed_time); 
	long threads_elapsed_time = 0; 

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

		for (i = 1; i < opt.nproc; i++) 
			threads_elapsed_time += procInfo[i].elapsed_time.tv_sec + procInfo[i].elapsed_time.tv_usec/1000000.0;
	}
	else 
	{
		gettimeofday(&single_start_time, NULL);

		if (conseq_anno_flag)
			io.readSeqFromFasta(vntrs);
		
		if (conseq_flag or raw_anno_flag)
			io.readSeqFromBam (vntrs, 1, 0, vntrs.size());

		if (conseq_flag)
			io.writeFa(out, vntrs);
		else 
		{
			int s = 0;
			for (auto &it: vntrs) 
			{
				ProcVNTR (s, it, opt);
				s += 1;
			}	
			cerr << "outputing vcf" << endl;
			io.writeVCFBody(out, vntrs, -1, 1);			
		}

		gettimeofday(&single_stop_time, NULL);
		timersub(&single_stop_time, &single_start_time, &single_elapsed_time); 
	}
	out.close();

	/* output vntr sequences to stdin */
	if (output_read_anno_flag)
	{
	    for (auto &vntr : vntrs)
	    {
	        for (auto &read : vntr->reads)
	        {
	            cout << ">"; 
	            cout.write(read->qname, read->l_qname);
	            cout << "\n";
	            cout.write(read->seq, read->len);
	            cout << "\n";
	        }
	    }
	}
	for (size_t i = 0; i < vntrs.size(); ++i) 
	{
		vntrs[i]->clearRead();
		delete vntrs[i];
	}
	
	if (opt.nproc > 1) {
		printf("[CPU time: %.2f sec, ", threads_elapsed_time + pre_elapsed_time.tv_sec + pre_elapsed_time.tv_usec/1000000.0);
	}
	else {
		printf("[CPU time: %.2f sec, ", single_elapsed_time.tv_sec + single_elapsed_time.tv_usec/1000000.0 + 
										 pre_elapsed_time.tv_sec + pre_elapsed_time.tv_usec/1000000.0);
	}

	double vm, rss;
	process_mem_usage(vm, rss);
	printf("RSS: %.2f G]\n", rss);

   	exit(EXIT_SUCCESS);
}

