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
int read_flag = false;
int somatic_flag = false;
int contig_flag = false;
int num_processed = 0;
int download_motifs=false;
int output_read_seq_flag=false;
// int seqan_flag = false;
// int output_read_flag = false;

struct timeval pre_start_time, pre_stop_time, pre_elapsed_time;
struct timeval single_start_time, single_stop_time, single_elapsed_time;


void PrintDownloadMotifs() {
  cout << "Most up-to-date motif set:" << endl

       << "For annotation on GRCh38:" << endl
       << " curl \"https://zenodo.org/records/11625069/files/vamos.motif.hg38.v2.1.e0.1.tsv.gz?download=1\" > vamos.motif.hg38.v2.1.e0.1.tsv.gz; gunzip vamos.motif.hg38.v2.1.e0.1.tsv.gz" << endl
       << "For annotation on CHM13" << endl
       << " curl \"https://zenodo.org/records/11625069/files/vamos.motif.CHM13.v2.1.orig.tsv.gz?download=1\" > vamos.motif.CHM13.v2.1.orig.tsv.gz; gunzip vamos.motif.CHM13.v2.1.orig.tsv.gz" << endl
       << " These motif sets have ~1.2M sites. They are an increase from v2.0 by adding loci " << endl
       << " that are in mobile elements." << endl    
       << "Previous versions, as well as the unfiltered motif sets, you can navigate to: https://zenodo.org/records/11625069" << endl;
  
    cout << "Motif sets from Ren, Gu, and Chaisson, Genome Biology, 2023 (for backwards compatibility):" << endl
       << "Original (no filtering) " <<endl
       << "   curl \"https://zenodo.org/records/13263615/files/vamos.oriMotifs.GRCh38.tsv.gz?download=1\" > vamos.oriMotifs.GRCh38.tsv.gz; gunzip vamos.oriMotifs.GRCh38.tsv.gz" << endl
       << "q10:" << endl
       << "   curl \"https://zenodo.org/records/13263615/files/vamos.effMotifs-0.1.GRCh38.tsv.gz?download=1\" > vamos.effMotifs-0.1.GRCh38.tsv.gz; gunzip vamos.effMotifs-0.1.GRCh38.tsv.gz" << endl
       << "q20:" << endl
       << "   curl \"https://zenodo.org/records/13263615/files/vamos.effMotifs-0.2.GRCh38.tsv.gz?download=1\" > vamos.effMotifs-0.2.GRCh38.tsv.gz; gunzip vamos.effMotifs-0.2.GRCh38.tsv.gz" << endl
       << "q30:" << endl
       << "   curl \"https://zenodo.org/records/13263615/files/vamos.effMotifs-0.3.GRCh38.tsv.gz?download=1\" > vamos.effMotifs-0.3.GRCh38.tsv.gz; gunzip vamos.effMotifs-0.3.GRCh38.tsv.gz" << endl << endl;
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


void ProcVNTR(int s, VNTR * it, const OPTION &opt, SDTables &sdTables, vector<int> &mismatchCI) 
{

  if (it->reads.size() == 0)
    {
        it->skip = true;
        if (debug_flag and it->reads.size() == 0) cerr << "no reads" << endl;
        return;
    }

    if (debug_flag) cerr << "start to do the annotation: " << s << endl;

    if (contig_flag) {
        it->motifAnnoForOneVNTR(opt, sdTables, mismatchCI); 
    }
    else if (read_flag) {
        it->consensusReadForHapByABpoa(opt);
        it->motifAnnoForOneVNTR(opt, sdTables, mismatchCI);      
    }
    else if (somatic_flag) {
        it->motifAnnoForOneVNTR(opt, sdTables, mismatchCI);
	it->clearReads();
    }
    
    return;
}

void SetEndPos(vector<VNTR*> &vntrs, map<string, vector<int> > &vntrMap, map<string, vector<int> > &endPos, int bufferNLoc) {
  for (auto mapIt : vntrMap) {
    string chrom=mapIt.first;
    //
    // Add the chromosome to the map of buckets.
    //
    if (endPos.find(chrom) == endPos.end()) {
      endPos[chrom] = vector<int>(1,1);
    }
    //
    // Skip empty buckets.
    //
    if (mapIt.second.size() == 0) {
      continue;
    }

    //
    // Create windows with bufferNLoc vntr loci per window.
    //
    int start=vntrs[mapIt.second[0]]->ref_start;
    endPos[chrom][0] = start;
    int cur = 0;
    while (cur < mapIt.second.size()) {
      int next = min(((int) mapIt.second.size())-1, cur+bufferNLoc);
      int nextEndPos = vntrs[mapIt.second[next]]->ref_end;
      endPos[chrom].push_back(nextEndPos);
      cur+=bufferNLoc;
    }
  }
}


void InitProcessed(map<string, vector<int> > &endPos, map<string, vector<char > > &processed) {
  for (auto epIt : endPos) {
    processed[epIt.first] = vector<char>(epIt.second.size()-1, false);
  }
}

bool GetNextUnprocessedRegion(int &curRegion, string &curChrom, map<string, vector<char> > &processed) {
  auto curChromIt = processed.find(curChrom);
  bool foundNext = false;
  while (curChromIt != processed.end() and foundNext == false) {
    //
    // First make sure this isn't pointing past the end of a region.
    // Exit early if so.
    while(curRegion < curChromIt->second.size() and
	  curChromIt->second[curRegion] == 1 ) {
      ++curRegion;
    }
    if (curRegion < curChromIt->second.size() and
	curChromIt->second[curRegion] == 0) {
      //
      // Mark this as being processed.
      //
      curChromIt->second[curRegion] = 1;
      curChrom = curChromIt->first;
      return true;
    }
    else {
      //
      // Start search at the beginning of the next region
      //
      curChromIt++;
      curRegion=0;
    }
  }
  return false;
}

void *CallSNVs (void *procInfoValue) {
  ProcInfo *procInfo = (ProcInfo *) procInfoValue;
  gettimeofday(&(procInfo->start_time), NULL);

  int curRegion = procInfo->thread;
  if (  (*(procInfo->processed)).size() == 0) {
    pthread_exit(NULL);    
  }
  string curChrom = (*(procInfo->processed)).begin()->first;
  bool readsArePhased = false;
  while (true) {
    procInfo->mtx->lock();
    bool foundNextRegion = false;
    foundNextRegion = GetNextUnprocessedRegion(curRegion, curChrom, (*(procInfo->processed)));
    if (foundNextRegion) {
      (*(procInfo->processed))[curChrom][curRegion] = 1;
    }
    procInfo->mtx->unlock();
    if (foundNextRegion == false) {
      pthread_exit(NULL);
    }
    Pileup pileup;
    (procInfo->io)->curChromosome = curChrom;
    if (readsArePhased == false) {
    
      (procInfo->io)->CallSNVs(curChrom,
			       (*(procInfo->bucketEndPos))[curChrom][curRegion],
			       (*(procInfo->bucketEndPos))[curChrom][curRegion+1],
			       *(procInfo->vntrs),
			       *(procInfo->vntrMap),
			       pileup, readsArePhased, *(procInfo->opt));
    }


    (procInfo->io)->StoreReadsOnChrom(curChrom,				      
				      (*(procInfo->bucketEndPos))[curChrom][curRegion],
				      (*(procInfo->bucketEndPos))[curChrom][curRegion+1],
				      *(procInfo->vntrs),
				      *(procInfo->vntrMap),
				      pileup, procInfo->thread, readsArePhased);
    vector<int>::iterator start, it,end;
     int regionStart = (*(procInfo->bucketEndPos))[curChrom][curRegion];
     int regionEnd   = (*(procInfo->bucketEndPos))[curChrom][curRegion+1];
     GetItOfOverlappingVNTRs(*(procInfo->vntrs),
			     (*(procInfo->vntrMap)), curChrom,
			     regionStart, regionEnd, start, end);

     SDTables sdTables;
     it = start;
     while (it != end and it != (*(procInfo->vntrMap))[curChrom].end()) {
       VNTR* vntr = (*procInfo->vntrs)[*it];
       ProcVNTR (vntr->index, vntr, *(procInfo->opt), sdTables, *(procInfo->mismatchCI));
       vntr->clearReads();
       it++;
     }
     cout << "Done processing " << curChrom << ":" << regionStart <<"-" << regionEnd << endl;
  }
  cerr << "Done reading " << procInfo->thread << endl;
  pthread_exit(NULL);     /* Thread exits (dies) */      
}


void *ProcVNTRs (void *procInfoValue)
{
    ProcInfo *procInfo = (ProcInfo *) procInfoValue;
    gettimeofday(&(procInfo->start_time), NULL);
    cerr << "start thread: " << procInfo->thread << endl;
    int i, s;
    int sz = (procInfo->vntrs)->size();
    SDTables sdTables;

    // read bam
    // procInfo->io->initializeBam();

    for (i = procInfo->thread, s = 0; i < sz; i += (procInfo->opt)->nproc, s += 1)
    {
        if (debug_flag) cerr << "processing vntr: " << i << endl;
        ProcVNTR (s, (*(procInfo->vntrs))[i], *(procInfo->opt), sdTables, *(procInfo->mismatchCI));
        (*procInfo->vntrs)[i]->clearReads();        
    }        
    
    gettimeofday(&(procInfo->stop_time), NULL);
    timersub(&(procInfo->stop_time), &(procInfo->start_time), &(procInfo->elapsed_time));
    pthread_exit(NULL);     /* Thread exits (dies) */    
}

void printUsage(IO &io, OPTION &opt) 
{

    printf("Usage: vamos [subcommand] [options] [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] \n");
    printf("Version: %s\n", io.version.c_str());
    printf("subcommand:\n");
    
    printf("vamos --contig [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] \n");
    printf("vamos --read [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] [-p phase_flank] \n");
    printf("vamos --somatic [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] [-p phase_flank] \n");    
    printf("vamos -m [verison of efficient motif set]\n");
    printf("\n");
    printf("   Input: \n");
    printf("       -b   FILE         Input indexed bam file. \n");    
    printf("       -r   FILE         File containing region coordinate and motifs of each VNTR locus. \n");
    printf("                         The file format: columns `chrom,start,end,motifs` are tab-delimited. \n");
    printf("                         Column `motifs` is a comma-separated (no spaces) list of motifs for this VNTR. \n");
    //    printf("       -R   FILE         Reference sequence where regions are on.\n");    
    printf("       -s   CHAR         Sample name. \n");
    printf("   Input handling:\n");
    printf("       -C   INT          Maximum coverage to call a tandem repeat.\n");
    printf("   Output: \n");
    printf("       -o   FILE         Output vcf file. \n");
    printf("       -S                Output assembly/read consensus sequence in each call.\n");
    printf("   Dynamic Programming: \n");
    printf("       -d   DOUBLE       Penalty of indel in dynamic programming (double) DEFAULT: 1.0. \n");
    printf("       -c   DOUBLE       Penalty of mismatch in dynamic programming (double) DEFAULT: 1.0. \n");
    printf("       --naive           Specify the naive version of code to do the annotation, DEFAULT: faster implementation. \n");
    printf("   Phase reads: \n");
    printf("       -p   INT            Range of flanking sequences which is used in the phasing step. DEFAULT: 15000 bps. \n");
    printf("       -M   INT            Minimum total coverage to allow a SNV to be called (6). \n");
    printf("       -a   INT            Minimum alt coverage to allow a SNV to be called (3). \n");        
    printf("   Downloading motifs:\n");
    PrintDownloadMotifs();
    printf("   Others: \n");
    printf("       -L   INT          Maximum length locus to compute annotation for (%d)\n", opt.maxLocusLength);
    printf("       -p   INT          Phase flank- how many bases on each side of a VNTR to collect SNVs to phase (default=15000)\n");
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
        {"naive",               no_argument,        &naive_flag,                    1},
        {"debug",               no_argument,        &debug_flag,                    1},
        {"readanno",            no_argument,        &output_read_anno_flag,         1},
        {"locuswise_prephase",  no_argument,        &locuswise_prephase_flag,       1},
        {"contig",              no_argument,        &contig_flag            ,       1}, // was locuswise_prephase_flag
        {"locuswise",           no_argument,        &locuswise_flag,                1},
        {"read",                no_argument,        &read_flag,                     1}, // was locuswise
        {"somatic",             no_argument,        &somatic_flag,                  1}, // was locuswise	
        {"single_seq",          no_argument,        &single_seq_flag,               1},
        {"readwise",            no_argument,        &readwise_anno_flag,            1},
        {"liftover",            no_argument,        &liftover_flag,                 1},
        {"output_seq",          no_argument,        &output_read_seq_flag,          1},	
        // {"clust",               no_argument,        &hclust_flag,                   1},
        // {"seqan",               no_argument,        &seqan_flag,                    1},
        // {"output_read",         no_argument,        &output_read_flag,              1},

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
        {"max_length",      required_argument,       0, 'L'},
        {"max_coverage",    required_argument,       0, 'C'},
        {"min_snv_coverage",   required_argument,       0, 'M'},
        {"min_alt_coverage",   required_argument,       0, 'a'},	
        {"penlaty_mismatch",required_argument,       0, 'c'},
        {"accuracy"        ,required_argument,       0, 'a'},
        {"phase_flank"     ,required_argument,       0, 'p'},
        {"download_db"     ,required_argument,       0, 'm'},
	//	{"reference"       ,required_argument,       0, 'R'},
	{"chry"            ,required_argument,       0, 'y'},	
        // {"input",           required_argument,       0, 'i'},
        // {"motif",           required_argument,       0, 'm'},

        {NULL, 0, 0, '\0'}
    };
    /* getopt_long stores the option index here. */
    int option_index = 0;
    while ((c = getopt_long (argc, argv, "Sb:r:a:o:C:s:t:f:d:c:x:v:m:R:y:p:hL:", long_options, &option_index)) != -1)
    {
        switch (c)
        {

            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0) break;
                fprintf (stderr, "option %s", long_options[option_index].name);
                if (optarg) fprintf (stderr, " with arg %s", optarg);
                break;
            case 'L':
	      fprintf(stderr, "option --max_length '\n", optarg);
	      io.maxLength=atoi(optarg);
	      opt.maxLocusLength = io.maxLength;
	      break;
            case 'C':
	      fprintf(stderr, "option --max_coverage '\n", optarg);
	      opt.maxCoverage=atoi(optarg);
	      break;	      
            case 'S':
	        fprintf (stderr, "option --output_seq '\n", optarg);
	        output_read_seq_flag = true;
		break;
            case 'b':
                fprintf (stderr, "option --bam with `%s'\n", optarg);
		io.input_bam = optarg;
                break;
		/*	    case 'R':
                fprintf (stderr, "option --reference with `%s'\n", optarg);
                io.reference = optarg;
                break;
		*/
            case 'v':
                fprintf (stderr, "option --vntr with `%s'\n", optarg);
                io.vntr_bed = optarg;
                break;
            case 'r':
                fprintf (stderr, "option --region with `%s'\n", optarg);
                io.region_and_motifs = optarg;
                break;

            // case 'i':
            //     fprintf (stderr, "option -input with `%s'\n", optarg);
            //     io.input_bam = (char *) malloc(strlen(optarg) + 1);
            //     strcpy(io.input_bam, optarg);
            //     break;

            // case 'm':
            //     fprintf (stderr, "option -motif with `%s'\n", optarg);
            //     io.motif_csv = (char *) malloc(strlen(optarg) + 1);
            //     strcpy(io.motif_csv, optarg);
            //     break;

            case 'o':
                fprintf (stderr, "option --output with `%s'\n", optarg);
                io.out_vcf = optarg;
                break;

            case 's':
                fprintf (stderr, "option --sampleName with `%s'\n", optarg);
                io.sampleName = optarg;
                break;

            case 't':
                fprintf (stderr, "option --numThreads with `%s'\n", optarg);
                opt.nproc = atoi(optarg);
                break;

            case 'y':
                fprintf (stderr, "option --chry with `%s'\n", optarg);
                opt.minChrY = atoi(optarg);
		io.minChrY = opt.minChrY;
                break;
            case 'd':
                opt.penalty_indel = stod(optarg);
                fprintf (stderr, "option --penlaty_indel with `%f'\n", opt.penalty_indel);
                break;

            case 'c':
                opt.penalty_mismatch = stod(optarg);
                fprintf (stderr, "option --penlaty_mismatch with `%f'\n", opt.penalty_mismatch);
                break;

            case 'f':
                fprintf (stderr, "option --filterNoisy\n");
                opt.filterStrength = stod(optarg);
                opt.filterNoisy = true;
                break;

            case 'a':
                opt.minAltCoverage = atoi(optarg);
                fprintf (stderr, "option -min_alt_coverage with `%d'\n", opt.minAltCoverage);
                break;
                    
            case 'M':
                opt.minSNVCoverage = atoi(optarg);
                fprintf (stderr, "option -min_snv_coverage with `%d'\n", opt.minSNVCoverage);
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
	      printUsage(io, opt);
                exit(EXIT_SUCCESS);

            case '?':
                fprintf(stderr, "Unknown option: %c\n", optopt);
                exit(EXIT_FAILURE);

            case ':':
                cerr << "[ERROR] missing option argument" << endl;
                exit(EXIT_FAILURE);

            default:
	      printUsage(io, opt);
                exit(EXIT_FAILURE);
        }        
    }
    if ( argc == 1) {
      printUsage(io, opt);
      exit(0);
    }
    if (read_flag) {
        locuswise_flag = true;
    }
    else if (somatic_flag) {
      locuswise_flag = false;
    }
    else if (contig_flag) {
      locuswise_flag = true;
      locuswise_prephase_flag = true;
    }
    else {
        cerr << "Either --contig or --read must be specified for input that is aligned contigs or reads" << endl;
        exit(1);
    }

    if (contig_flag) {
        opt.inputType=by_contig;
    }

    if (locuswise_flag) io.phaseFlank = opt.phaseFlank;

    if (download_motifs) {
        PrintDownloadMotifs();
        exit(0);
    }

    /* Check mandatory parameters */
    bool missingArg = false;
    if (io.input_bam == "") 
    {
        fprintf(stderr, "ERROR: -b must be specified!\n");
        missingArg = true;
    }
    if (io.region_and_motifs == "")
    {
        fprintf(stderr, "ERROR:-r must be specified!\n");
        missingArg = true;
    }
    // if (io.motif_csv == NULL and (conseq_anno_flag or locuswise_prephase_flag))
    // {
    //     fprintf(stderr, "-m is mandatory with conseq and locuswise!\n");
    //     missingArg = true;
    // }
    if (io.sampleName == "")
    {
        fprintf(stderr, "ERROR: -s must be specified!\n");
        missingArg = true;
    }
    if (io.out_vcf == "")
    {
        fprintf(stderr, "ERROR: -o must be specified!\n");
        missingArg = true;
    }

    if (missingArg)
    {
      printUsage(io, opt);
        exit(EXIT_FAILURE);
    }

    if (naive_flag) fprintf(stderr, "naive_flag is set. \n");
    if (debug_flag) fprintf(stderr, "debug_flag is set. \n");
    // if (hclust_flag) fprintf(stderr, "hclust_flag is set. \n");
    // if (seqan_flag) fprintf(stderr, "seqan_flag is set\n");
    if (output_read_anno_flag) fprintf(stderr, "output_read_anno_flag is set. \n");
    if (single_seq_flag) fprintf(stderr, "single_seq_flag is set\n");
    if (readwise_anno_flag) fprintf(stderr, "readwise_anno_flag is set. \n");
    if (locuswise_prephase_flag) fprintf(stderr, "locuswise_prephase_flag is set. \n");
    

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
    std::map<string, vector<int> > vntrMap;
    io.vntrMap = &vntrMap;
    if (io.region_and_motifs != "")
    {
        io.readRegionAndMotifs(vntrs);
	ChromToVNTRMap(vntrs, vntrMap);
        CreateAccLookupTable(vntrs, 0.98, mismatchCI, 0.999);      
    }
    cerr << "Read " << vntrs.size() << " target loci" << endl;

    
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


    //
    // Read input.
    //
    if (contig_flag) {
      io.StoreAllContigs(vntrs, vntrMap);
      int i;
      if (opt.nproc > 1)    {
        pthread_t *tid = new pthread_t[opt.nproc];
	
        mutex mtx;
        mutex ioLock;
        vector<ProcInfo> procInfo(opt.nproc);
	
        for (i = 0; i < opt.nproc; i++) { 
	  procInfo[i].vntrs = &vntrs;
	  procInfo[i].thread = i;
	  procInfo[i].opt = &opt;
	  // procInfo[i].io = &io;
	  procInfo[i].io = new IO;
	  procInfo[i].io->ioLock = &ioLock;
	  procInfo[i].io->bai = io.bai;
	  procInfo[i].io->fp_in = io.fp_in;
	  procInfo[i].io->bamHdr = io.bamHdr;
	  procInfo[i].io->idx = io.idx;
	  procInfo[i].io->numProcessed = &num_processed;
	  procInfo[i].io->minChrY = io.minChrY;
	  procInfo[i].mtx = &mtx;
	  procInfo[i].mismatchCI = &mismatchCI;
	  procInfo[i].out = &out;
	  // procInfo[i].out_nullAnno = &out_nullAnno;

	  if (pthread_create(&tid[i], NULL, ProcVNTRs, (void *) &procInfo[i]) )
            {
	      cerr << "ERROR: Cannot create thread" << endl;
	      exit(EXIT_FAILURE);
            }
        }

        for (i = 0; i < opt.nproc; i++)
	  {
            pthread_join(tid[i], NULL);
            delete procInfo[i].io;
	  }
        delete[] tid;
      
        for (i = 1; i < opt.nproc; i++) 
	  threads_elapsed_time += procInfo[i].elapsed_time.tv_sec + procInfo[i].elapsed_time.tv_usec/1000000.0;
      }
      else {
	gettimeofday(&single_start_time, NULL);
	
	// read input bam/fasta
	if (single_seq_flag)
	  io.readSeqFromFasta(vntrs);
	
	int s = 0;
	SDTables sdTables;
	for (auto i=0; i < vntrs.size(); i++)  { 
	  ProcVNTR (s, vntrs[i], opt, sdTables, mismatchCI);
	    vntrs[i]->clearReads();	    
	}
          
        gettimeofday(&single_stop_time, NULL);
        timersub(&single_stop_time, &single_start_time, &single_elapsed_time);
      }
    }
    else {
      //
      // Need to read all reads, do parallel processing of contigs.
      //
      int maxIOThread=opt.nproc;
      if (maxIOThread > 1) {
	  pthread_t *tid = new pthread_t[maxIOThread];
	  vector<ProcInfo> procInfo(maxIOThread);	  
	  vector<bool> procChrom(io.chromosomeNames.size(), false);
	  vector<Pileup> pileups(io.chromosomeNames.size());

	  map<string, vector<int> > bucketEndPos;
	  map<string, vector<char> > processed;
	  
	  SetEndPos(vntrs, vntrMap, bucketEndPos, 1000);
	  InitProcessed(bucketEndPos, processed);
	  
	  mutex ioLock, mtx;
	  for (int i = 0; i < maxIOThread; i++) { 
            procInfo[i].vntrs = &vntrs;
            procInfo[i].vntrMap = &vntrMap;
	    procInfo[i].pileups = &pileups;
            procInfo[i].thread = i;
            procInfo[i].opt = &opt;
            procInfo[i].io = new IO;
            procInfo[i].io->input_bam = io.input_bam;
	    //	    procInfo[i].io->reference = io.reference;
	    procInfo[i].io->initializeBam();
	    procInfo[i].io->initializeRefFasta();	    
	    procInfo[i].io->chromosomeNames = io.chromosomeNames;
            procInfo[i].io->ioLock = &ioLock;
            procInfo[i].io->numProcessed = &num_processed;
	    procInfo[i].io->thread = i;
	    procInfo[i].io->minChrY = opt.minChrY;
            procInfo[i].mtx = &mtx;
            procInfo[i].mismatchCI = &mismatchCI;
            procInfo[i].out = &out;
	    procInfo[i].processed = &processed;
	    procInfo[i].bucketEndPos = &bucketEndPos;	    
            // procInfo[i].out_nullAnno = &out_nullAnno;
            if (pthread_create(&tid[i], NULL, CallSNVs, (void *) &procInfo[i]) )
            {
                cerr << "ERROR: Cannot create thread" << endl;
                exit(EXIT_FAILURE);
            }
	    
	  }
	  for (int i = 0; i < maxIOThread; i++)
	    {
	      pthread_join(tid[i], NULL);
	      delete procInfo[i].io;
	    }	  
      }
      else {
        io.initializeBam();
	io.initializeRefFasta();
	bool readsArePhased = false;
	for (int i=0; i < io.chromosomeNames.size(); i++ ) {
	  Pileup pileup;
	  if ( vntrMap.find(io.chromosomeNames[i] ) != vntrMap.end() ) {
	    int start = vntrs[vntrMap[io.chromosomeNames[i]][0]]->ref_start;
	    int end   = vntrs[vntrMap[io.chromosomeNames[i]][vntrMap[io.chromosomeNames[i]].size()-1]]->ref_end;
	    if (readsArePhased == false) {
	      io.CallSNVs(io.chromosomeNames[i], start, end, vntrs, vntrMap, pileup, readsArePhased, opt);
	    }
	    io.StoreReadsOnChrom(io.chromosomeNames[i], start, end, vntrs, vntrMap, pileup, 1, readsArePhased);

	    SDTables sdTables;
	    for (int j=0; j < vntrMap[io.chromosomeNames[i]].size(); j++) {
	      int v=vntrMap[io.chromosomeNames[i]][j];
	      ProcVNTR (vntrs[v]->index, vntrs[v], opt, sdTables, mismatchCI);
	      vntrs[v]->clearReads();
	    }	    
	}
      }
    }
    }
    
      /* Create threads */

    //    io.clear();
        // output vcf or bed or fasta
    if (readwise_anno_flag or somatic_flag) 
      io.writeBEDBody_readwise(out, vntrs, -1, 1);
    else if (locuswise_prephase_flag or locuswise_flag or single_seq_flag) 
      io.writeVCFBody_locuswise(out, vntrs, -1, 1);
    
    out.close();

    /* output vntr sequences to stdout */
    if (output_read_anno_flag)
    {
        for (auto &vntr : vntrs)
        {
            for (auto &read : vntr->reads)
            {
	      cout << ">" << read->qname << endl;
	      cout << read->seq << endl;
            }
        }
    }

    if (opt.nproc > 1)
    {
        printf("[CPU time: %.2f sec, ", threads_elapsed_time + pre_elapsed_time.tv_sec + pre_elapsed_time.tv_usec/1000000.0);
    }
    else
    {
        printf("[CPU time: %.2f sec, ", single_elapsed_time.tv_sec + single_elapsed_time.tv_usec/1000000.0 + 
            pre_elapsed_time.tv_sec + pre_elapsed_time.tv_usec/1000000.0);
    }

    double vm, rss;
    process_mem_usage(vm, rss);
    printf("RSS: %.2f G]\n", rss);

    exit(EXIT_SUCCESS);
}

       
