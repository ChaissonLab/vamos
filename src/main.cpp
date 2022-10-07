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
#include "acc_lookup_table.h"
#include <mutex>   
#include <unistd.h>
#include <iomanip>
#include <errno.h>

// using namespace std;
int naive_flag = false;
int debug_flag = false;
int hclust_flag = false;
int output_read_anno_flag = false;
int readwise_anno_flag = false;
int liftover_flag = false;
int single_seq_flag = false;
int locuswise_prephase_flag = false;
int locuswise_flag = false;
int num_processed = 0;
int download_motifs=false;
// int seqan_flag = false;
// int output_read_flag = false;

struct timeval pre_start_time, pre_stop_time, pre_elapsed_time;
struct timeval single_start_time, single_stop_time, single_elapsed_time;


void PrintDownloadMotifs(OPTION &opts) {
  if (opts.download == "original") {  
    fprintf(stdout, "please run:\ncurl \"https://zenodo.org/record/7155334/files/processed_vntrs.tsv?download=1\" > original_motifs.bed\n");
  }
  if (opts.download == "q10") {
    fprintf(stdout,"please run:\ncurl \"https://zenodo.org/record/7155329/files/vntrs_motifs_delta_0.1.bed?download=1\"  > vntrs_motifs_delta_0.1.bed\n");
  }
  if (opts.download == "q20") {
    fprintf(stdout,"please run:\ncurl \"https://zenodo.org/record/7155329/files/vntrs_motifs_delta_0.2.bed?download=1\"  > vntrs_motifs_delta_0.2.bed\n");
  }
  if (opts.download == "q30") {
    fprintf(stdout,"please run:\ncurl \"https://zenodo.org/record/7155329/files/vntrs_motifs_delta_0.3.bed?download=1\"  > vntrs_motifs_delta_0.3.bed\n");
  }
  else {
    fprintf(stderr, "download option '%s' is not supported. Please use one of:\n"
	    "   q10   q=0.10, 64-haplotypes (Ebert 2021).\n" 
	    "   q20   q=0.20, 64-haplotypes (Ebert 2021).\n" 
	    "   q30   q=0.30, 64-haplotypes (Ebert 2021).\n", opts.download.c_str());      
    exit(1);
  }
}

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


void ProcVNTR (int s, VNTR * it, const OPTION &opt, SDTables &sdTables, vector< int > &mismatchCI) 
{
  
	if (it->nreads == 0 or it->motifs.size() > 255) 
	{
		it->skip = true;
		if (debug_flag and it->nreads == 0) cerr << "no reads" << endl;
		if (debug_flag and it->motifs.size() > 255) cerr << "skip one vntr due to > 255 motifs" << endl;
		return;
	}

	if (debug_flag) cerr << "start to do the annotation: " << s << endl;

	if (!locuswise_flag) 
	{
		it->motifAnnoForOneVNTR(opt, sdTables, mismatchCI); 
		
		if (!readwise_anno_flag) 	
			it->consensusMotifAnnoForOneVNTRByABpoa(opt);	

		// it->annoTostring(opt);
		// it->cleanNoiseAnno(opt); // (TODO: remove this)
		// if (hclust_flag)
		// 	it->consensusMotifAnnoForOneVNTRByClustering(opt);
		// else if (!readwise_anno_flag) 
		// 	it->consensusMotifAnnoForOneVNTRByABpoa(opt);
	}
	else 
	{
		it->consensusReadForHapByABpoa(opt);
		it->motifAnnoForOneVNTR(opt, sdTables, mismatchCI); 
		// it->annoTostring(opt);
		// it->cleanNoiseAnno(opt);		
	}
	
	return;
}

void *ProcVNTRs (void *procInfoValue)
{

	ProcInfo *procInfo = (ProcInfo *)procInfoValue;
	gettimeofday(&(procInfo->start_time), NULL);
	cerr << "start thread: " << procInfo->thread << endl;
	int i, s;
	int sz = (procInfo->vntrs)->size();
	SDTables sdTables;

	// read bam 
	(procInfo->io)->readSeqFromBam ((*(procInfo->vntrs)), (procInfo->opt)->nproc, procInfo->thread, sz);

	for (i = procInfo->thread, s = 0; i < sz; i += (procInfo->opt)->nproc, s += 1)
	{
		if (debug_flag) cerr << "processing vntr: " << i << endl;
		// procInfo->numOfProcessed += ProcVNTR (s, (*(procInfo->vntrs))[i], *(procInfo->opt));
		ProcVNTR (s, (*(procInfo->vntrs))[i], *(procInfo->opt), sdTables, *(procInfo->mismatchCI));
	}		
	
	procInfo->mtx->lock();
	cerr << "outputing vcf" << endl;	
	if (locuswise_prephase_flag or locuswise_flag) 
		(procInfo->io)->writeVCFBody_locuswise(*(procInfo->out), (*(procInfo->vntrs)), procInfo->thread, (procInfo->opt)->nproc);	
	else if (readwise_anno_flag) 
		(procInfo->io)->writeBEDBody_readwise(*(procInfo->out), (*(procInfo->vntrs)), procInfo->thread, (procInfo->opt)->nproc);	
	procInfo->mtx->unlock();

	cerr << "finish thread: " << procInfo->thread << endl;
	gettimeofday(&(procInfo->stop_time), NULL);
	timersub(&(procInfo->stop_time), &(procInfo->start_time), &(procInfo->elapsed_time)); 
	pthread_exit(NULL);     /* Thread exits (dies) */	
}

void printUsage(IO &io) 
{

	printf("Usage: vamos [subcommand] [options] [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] \n");
	printf("Version: %s\n", io.version);
	printf("subcommand:\n");
    
	// printf("vamos --liftover   [-i in.bam] [-v vntrs.bed] [-o output.fa] [-s sample_name] (ONLY FOR SINGLE LOCUS!!) \n");
	// printf("vamos --readwise   [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.bed] [-s sample_name] [-t threads] \n");	
	// printf("vamos --locuswise_prephase  [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] \n");
	// printf("vamos --locuswise [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] [-p phase_flank]\n");
	// printf("vamos --single_seq [-b in.fa]  [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] (ONLY FOR SINGLE LOCUS!!) \n");
	printf("vamos --contig [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] \n");
	printf("vamos --read [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] [-p phase_flank] \n");
	printf("\n");
	printf("   Input: \n");
	printf("       -b   FILE         Input indexed bam file. \n");	
	printf("       -r   FILE         File containing region coordinate and motifs of each VNTR locus. \n");
	printf("                         The file format: columns `chrom,start,end,motifs` are tab-delimited. \n");
	printf("                         Column `motifs` is a comma-separated (no spaces) list of motifs for this VNTR. \n");
	printf("       -s   CHAR         Sample name. \n");
	printf("   Output: \n");
	printf("       -o   FILE         Output vcf file. \n");
	printf("   Dynamic Programming: \n");
	printf("       -d   DOUBLE       Penalty of indel in dynamic programming (double) DEFAULT: 1.0. \n");
	printf("       -c   DOUBLE       Penalty of mismatch in dynamic programming (double) DEFAULT: 1.0. \n");
	printf("       -a   DOUBLE       Global accuracy of the reads. DEFAULT: 0.98. \n");	
	printf("       --naive           Specify the naive version of code to do the annotation, DEFAULT: faster implementation. \n");
	// printf("   Aggregate Annotation: \n");
	// printf("       -f   DOUBLE       Filter out noisy read annotations, DEFAULT: 0.0 (no filter). \n");
	// printf("       --clust           use hierarchical clustering to judge if a VNTR locus is het or hom. \n");
	printf("   Phase reads: \n");
	printf("       -p   INT       	 Range of flanking sequences which is used in the phasing step. DEFAULT: 3000 bps. \n");
	printf("   Downloading motifs:\n");
	printf("       -m  MOTIF         Prints a command to download a particular motif set. Current supported motif set is: d10e32. \n"
	       "                         This motif set is selected at a level of Delta=10 from 32 haplotype-resolvd assemblies (Ebert et al., 2021)\n"
	       "                         This may be copied and pasted in the command line, or executed as: vamos -m d10e32\n");
	printf("   Others: \n");
	printf("       -t   INT          Number of threads, DEFAULT: 1. \n");
	printf("       --debug           Print out debug information. \n");
	printf("       -h                Print out help message. \n");
	

	// printf("       --seqan           use seqan lib to do MSA (haploid only), DEFAULT: abPoa\n");
	// printf("       --readanno        output read annotation in VCF and output vntr sequences to stdout. \n");
} 

int main (int argc, char **argv)
{
  
	int c;
	IO io;
	OPTION opt;

	const struct option long_options[] =
	{
		/* These options set a flag. */
		{"naive",               no_argument,             &naive_flag,                    1},
		{"debug",               no_argument,             &debug_flag,                    1},
		{"readanno",            no_argument,             &output_read_anno_flag,         1},
		{"locuswise_prephase",  no_argument,             &locuswise_prephase_flag,       1},
		{"contig",            no_argument,               &locuswise_prephase_flag,       1},		
		{"locuswise",           no_argument,             &locuswise_flag,                1},
		{"read",               no_argument,              &locuswise_flag,                1},		
		{"single_seq",          no_argument,             &single_seq_flag,               1},
		{"readwise",            no_argument,             &readwise_anno_flag,            1},
		{"liftover",            no_argument,             &liftover_flag,                 1},
		// {"clust",               no_argument,             &hclust_flag,                   1},
		// {"seqan",         no_argument,             &seqan_flag,                    1},
		// {"output_read",   no_argument,             &output_read_flag,              1},		


		/* These options don’t set a flag. We distinguish them by their indices. */
		{"bam",             required_argument,       0, 'b'},
		{"region",          required_argument,       0, 'r'},						
		{"vntr",            required_argument,       0, 'v'},
		{"output",          required_argument,       0, 'o'},
		{"out_fa",          required_argument,       0, 'x'},
		{"sampleName",      required_argument,       0, 's'},
		{"numThreads",      required_argument,       0, 't'},
		{"filterNoisy",     required_argument,       0, 'f'},
		{"penlaty_indel",   required_argument,       0, 'd'},
		{"penlaty_mismatch",required_argument,       0, 'c'},
		{"accuracy"        ,required_argument,       0, 'a'},
		{"phase_flank"     ,required_argument,       0, 'p'},
        {"download_db"     ,required_argument,       0, 'm'},
		// {"input",           required_argument,       0, 'i'},
		// {"motif",           required_argument,       0, 'm'},

		{NULL, 0, 0, '\0'}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;
	while ((c = getopt_long (argc, argv, "b:r:a:o:s:t:f:d:c:x:v:m:p:h", long_options, &option_index)) != -1)
	{
		switch (c)
		{

        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0) break;
          fprintf (stderr, "option %s", long_options[option_index].name);
          if (optarg) fprintf (stderr, " with arg %s", optarg);
          break;

		case 'b':
			fprintf (stderr, "option -input with `%s'\n", optarg);
			io.input_bam = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.input_bam, optarg);
			io.input_bam[strlen(optarg)] = '\0';						
			break;

		case 'v':
			fprintf (stderr, "option -vntr with `%s'\n", optarg);
			io.vntr_bed = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.vntr_bed, optarg);
			io.vntr_bed[strlen(optarg)] = '\0';			
			break;

		case 'r':
			fprintf (stderr, "option -region with `%s'\n", optarg);
			io.region_and_motifs = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.region_and_motifs, optarg);
			io.region_and_motifs[strlen(optarg)] = '\0';
			io.region_and_motifs[strlen(optarg)] = '\0';						
			break;


		// case 'i':
		// 	fprintf (stderr, "option -input with `%s'\n", optarg);
		// 	io.input_bam = (char *) malloc(strlen(optarg) + 1);
		// 	strcpy(io.input_bam, optarg);
		// 	break;

		// case 'm':
		// 	fprintf (stderr, "option -motif with `%s'\n", optarg);
		// 	io.motif_csv = (char *) malloc(strlen(optarg) + 1);
		// 	strcpy(io.motif_csv, optarg);
		// 	break;

		case 'o':
			fprintf (stderr, "option -output with `%s'\n", optarg);
			io.out_vcf = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.out_vcf, optarg);
			io.out_vcf[strlen(optarg)] = '\0';						
			break;

		case 's':
			fprintf (stderr, "option -sampleName with `%s'\n", optarg);
			io.sampleName = (char *) malloc(strlen(optarg) + 1);
			strcpy(io.sampleName, optarg);
			io.sampleName[strlen(optarg)] = '\0';						
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

		case 'a':
 		    opt.accuracy = stod(optarg);
			fprintf (stderr, "option -accuracy with `%f'\n", opt.accuracy);
			break;
			
		case 'p':
            opt.phaseFlank = atoi(optarg);
            io.phaseFlank = opt.phaseFlank;
            fprintf (stderr, "option -phase_flank with `%d'\n", opt.phaseFlank);
            break;

        case 'm':
            opt.download=optarg;
			download_motifs=true;
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

	if (locuswise_flag) io.phaseFlank = opt.phaseFlank;

	if (download_motifs) {
	  	PrintDownloadMotifs(opt);
        exit(0);
	}

	/* Check mandatory parameters */
	bool missingArg = false;
	if (io.input_bam == NULL) 
	{
		fprintf(stderr, "ERROR: -b must be specified!\n");
		missingArg = true;
	}
	if (!liftover_flag and io.region_and_motifs == NULL)
	{
		fprintf(stderr, "ERROR:-r must be specified!\n");
		missingArg = true;
	}
	// if (io.motif_csv == NULL and (conseq_anno_flag or locuswise_prephase_flag))
	// {
	// 	fprintf(stderr, "-m is mandatory with conseq and locuswise!\n");
	// 	missingArg = true;
	// }
	if (io.sampleName == NULL)
	{
		fprintf(stderr, "ERROR: -s must be specified!\n");
		missingArg = true;
	}
	if (io.out_vcf == NULL)
	{
		fprintf(stderr, "ERROR: -o must be specified!\n");
		missingArg = true;
	}

	if (missingArg)
	{
		printUsage(io);
		exit(EXIT_FAILURE);
	}

  	if (naive_flag) fprintf(stderr, "naive_flag is set. \n");
  	if (debug_flag) fprintf(stderr, "debug_flag is set. \n");
   	// if (hclust_flag) fprintf(stderr, "hclust_flag is set. \n");
  	// if (seqan_flag) fprintf(stderr, "seqan_flag is set\n");
   	if (output_read_anno_flag) fprintf(stderr, "output_read_anno_flag is set. \n");
  	if (liftover_flag) fprintf(stderr, "liftover_flag is set\n");
  	if (single_seq_flag) fprintf(stderr, "single_seq_flag is set\n");
  	if (readwise_anno_flag) fprintf(stderr, "readwise_anno_flag is set. \n");
  	if (locuswise_prephase_flag) fprintf(stderr, "locuswise_prephase_flag is set. \n");
  	if (locuswise_flag) fprintf(stderr, "locuswise_flag is set. \n");
	

	/* Print any remaining command line arguments (not options). */
	if (optind < argc)
	{
		fprintf (stderr, "non-option ARGV-elements: ");
		while (optind < argc) fprintf (stderr, "%s ", argv[optind++]);
		putchar ('\n');
	}

	vector<VNTR *> vntrs;
	gettimeofday(&pre_start_time, NULL);

	// read input_bam and region_and_motifs
	vector< int> mismatchCI;	
	if (io.region_and_motifs != NULL) {
		io.readRegionAndMotifs(vntrs);
		CreateAccLookupTable(vntrs, opt.accuracy, mismatchCI, 0.999);	  
	}
	if (liftover_flag) {
		io.readVNTRFromBed(vntrs);
	}
	cerr << "finish reading " << vntrs.size() << " vntrs" << endl;
	
	/* set up out stream and write VCF header */
	ofstream out;
	out.open(io.out_vcf, ofstream::out);
	if (out.fail()) 
	  {
	    cerr << "ERROR: Unable to open file " << io.out_vcf << endl;
	    exit(EXIT_FAILURE);
	  } 	
	if (readwise_anno_flag)
		io.writeBEDHeader_readwise(out);
	else if (locuswise_prephase_flag or locuswise_flag or single_seq_flag)
		io.writeVCFHeader_locuswise(out);

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
			procInfo[i].mismatchCI = &mismatchCI;
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

		// read input bam/fasta
		if (single_seq_flag)
			io.readSeqFromFasta(vntrs);
		else 
			io.readSeqFromBam (vntrs, 1, 0, vntrs.size());

		if (!liftover_flag) {
			int s = 0;
			SDTables sdTables;
			for (auto &it: vntrs) 
			{
				ProcVNTR (s, it, opt, sdTables, mismatchCI);
				s += 1;
			}				
		}
	
		// output vcf or bed or fasta
		if (readwise_anno_flag) 
		 	io.writeBEDBody_readwise(out, vntrs, -1, 1);
		else if (locuswise_prephase_flag or locuswise_flag or single_seq_flag) 
			io.writeVCFBody_locuswise(out, vntrs, -1, 1);
		else if (liftover_flag)
			io.writeFa(out, vntrs);
		
		gettimeofday(&single_stop_time, NULL);
		timersub(&single_stop_time, &single_start_time, &single_elapsed_time); 
	}

	out.close();

	/* output vntr sequences to stdout */
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

