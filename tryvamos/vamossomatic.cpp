#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string.hpp>  // for boost::split
#include <iomanip>

namespace bio = boost::iostreams;
using namespace std;

void process_line(const std::string& line) {
    std::vector<std::string> columns;
    boost::split(columns, line, boost::is_any_of("\t"));

    for (const auto& col : columns) {
        std::cout << col << " ";
    }
    std::cout << std::endl;
}


std::vector<std::string> split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string item;

    // Use std::getline to extract substrings delimited by 'delimiter'
    while (std::getline(ss, item, delimiter)) {
        tokens.push_back(item);
    }

    return tokens;
}

class SampleLocus {
public:
  string chrom="";
  int start=-1;
  int end=-1;
  vector<float> means, vars, iqr1s, iqr2s;

  vector<vector< float> > iqrExtremes;
  int phase;

  int MatchesRegion(string qchrom, int qstart, int qend) {
    return chrom == qchrom and start == qstart and end == qend;
  }
  
  int ParseLine(string &line, string startOn="") {

    stringstream strm(line);
    strm >> chrom >> start >> end;
    if (chrom == "") { return 0;}
    if (startOn != "" and chrom!= startOn) {
      return 1;
    }
    string motifStr;
    strm >> motifStr;
    if (motifStr == "") { return 0;}

    string annoStr;
    strm >> annoStr;
    if (annoStr == "") { return 0;}
    vector<string> annoStrs;
    annoStrs = split(annoStr, ';');
    vector<vector<int> > lengths(3);
    for (auto i=0; i < annoStrs.size()-1; i++) {
      vector<string> readStr = split(annoStrs[i], ':');
      int hap=atoi(readStr[0].c_str());
      int len=atoi(readStr[1].c_str());
      assert(hap >= 0 and hap <= 2);
      lengths[hap].push_back(len);
    }

    if (lengths[0].size() < 3 and lengths[1].size() > 3 and lengths[2].size() > 3) {
      lengths[0] = lengths[1];
      lengths[1] = lengths[2];
      phase = 1;
    }
    else if (lengths[0].size() > 3 and lengths[1].size() < 3 and lengths[2].size() < 3) {
      lengths[1] = lengths[0];
      phase = 0;
    }
    else {
      //      cout << "SHORT " << lengths[0].size() << "\t" << lengths[1].size() << "\t" << lengths[2].size() << endl;
      lengths.resize(0);
      phase = -1;
      return 1;
    }
    lengths.resize(2);
    for (auto &vect : lengths ) {
      sort(vect.begin(), vect.end());
    }
    means.resize(2); vars.resize(2); iqr1s.resize(2); iqr2s.resize(2);
    iqrExtremes.resize(2);
    iqrExtremes[0].resize(0);
    iqrExtremes[1].resize(0);
    for (auto i =0; i < lengths.size(); i++) {
      if (lengths[i].size() == 0) {
	means[i] = vars[i] = iqr1s[i] = iqr2s[i] = 0;
	continue;
      }
      int sum=0;
      for (auto v: lengths[i]) {
	sum+=v;
      }
      means[i] = ((float)sum)/lengths[i].size();
      //
      // Mean is normalized by length.
      //
      vector<float> lenNorm(lengths[i].size());
      float normSum=0, normSumSq=0;
      
      for (auto j = 0; j < lenNorm.size(); j++) {
	lenNorm[j] = lengths[i][j]/means[i];
	normSum+=lenNorm[j];
	normSumSq += lenNorm[j] * lenNorm[j];
      }
      float normMean = normSum / lengths[i].size();
      vars[i] = normSumSq/lengths[i].size() - normMean * normMean;

      //
      // Store iqrs;
      //
      int q1 = (int)lengths[i].size()*0.25;
      int q3 = (int)lengths[i].size()*0.75;
      iqr1s[i] = lengths[i][q1];
      iqr2s[i] = lengths[i][q3];
      for (auto v: lengths[i]) {
	float iqrSpan = iqr2s[i] - iqr1s[i];
	if (iqrSpan > 0.001 and v > iqr2s[i] + (iqrSpan*3)) {
	  iqrExtremes[i].push_back(v);
	}
      }
    }	
    //    cout << "LENGTHS:\t" << lengths[0].size() << "\t" << lengths[1].size() << endl;
    return 1;
  }
};

class Sample {
public:
  string filename;
  bio::filtering_streambuf<bio::input> fileIn;
  std::istream *inStream;
  SampleLocus locus;
  int Init(string fn) {
    filename = fn;
    fileIn.push(bio::gzip_decompressor());
    fileIn.push(bio::file_source(filename, std::ios_base::in | std::ios_base::binary));
    
    inStream= new std::istream(&fileIn);
    return inStream->good();
  }
  int GetNext(string startOn="") {
    string line;
    getline(*inStream, line);
    int retval;
    if (line == "") { return 0;}
    retval = locus.ParseLine(line, startOn);
    return retval;
  }
  int MatchesRegion(string qchrom, int qstart, int qend) {
    return locus.MatchesRegion(qchrom, qstart, qend);
  }
  
};

int main(int argc, char* argv[]) {
  if (argc < 5) {
    cout << "Usage vamossomatic statsFile iqrFile regions.tsv tab1.gz tab2.gz ..." << endl;
    exit(0);
  }
  string statsFileName=argv[1];
  string iqrFileName=argv[2];
  ofstream statsFile(statsFileName.c_str());
  ofstream iqrFile(iqrFileName.c_str());
  
  std::string regionFileName = argv[3];  // Your gzip-compressed file

  // Create a file source and gzip filter
  ifstream regionIn(regionFileName.c_str());
  vector<string> regionChrom;
  vector<int>    regionStart, regionEnd;
  cout << "Reading regions" << endl;
  while(regionIn) {
    string line;
    string chr;
    int start, end;
    regionIn >> chr >> start >> end;
    if (chr != "") {
      regionChrom.push_back(chr);
      regionStart.push_back(start);
      regionEnd.push_back(end);
    }
    else {
      break;
    }
    getline(regionIn, line);
  }
  cout << "Initializing samples" << endl;
  vector<Sample> samples(argc - 4);
  for (auto i=4; i < argc; i++) {
    samples[i-4].Init(argv[i]);
  }
  int numOutput=0;
  int locusIndex=0;
  for (auto s=0; s < samples.size(); s++) {
    //
    // Initialize each sample.
    //
    int retval;
    retval = samples[s].GetNext();
    if (retval) {
      //
      // Allow for the case that the table starts in the middle of the vamos output.
      //
      int nSkipped=0;
      while(!samples[s].MatchesRegion(regionChrom[0], regionStart[0], regionEnd[0])) {
	retval = samples[s].GetNext(regionChrom[0]);
	if (retval == 0) {
	  break;
	}
	nSkipped++;
      }
      if (nSkipped > 0) {
	cout << "Sample " << s << " skipped " << nSkipped << " starting on " << samples[s].locus.chrom << ":" << samples[s].locus.start << "-" << samples[s].locus.end << endl;
      }
    }
  }
  
  for (locusIndex=0; locusIndex < regionChrom.size(); locusIndex++) {
    //
    // Advance all samples until the next region, possibly.
    //
    statsFile << regionChrom[locusIndex] << "\t" <<  regionStart[locusIndex] << "\t" <<  regionEnd[locusIndex] << endl;
    iqrFile   << regionChrom[locusIndex] << "\t" <<  regionStart[locusIndex] << "\t" <<  regionEnd[locusIndex] << endl;    
    for (auto s=0; s < samples.size(); s++) {
      //
      // Advance past last locus if this was a match. If not a match, this
      // locus must be one later in the table.
      //
      if (locusIndex > 0 and 
	  samples[s].MatchesRegion(regionChrom[locusIndex-1], regionStart[locusIndex-1], regionEnd[locusIndex-1])) {
	samples[s].GetNext();
      }
      /*      
      else {
	cout << "SAMPLE " << s << " is waiting for others to catch up." << endl;
      }
      */
    }
    //
    // Output vars
    if ((locusIndex + 1) % 1000 == 0) {
      cout << locusIndex +1 << "\t" << regionChrom[locusIndex] << "\t" <<  regionStart[locusIndex] << "\t" << numOutput << endl;
    }
    int numNotNA=0;
    for (auto s=0; s < samples.size(); s++) {
      if (samples[s].MatchesRegion(regionChrom[locusIndex], regionStart[locusIndex], regionEnd[locusIndex]) and
	  samples[s].locus.phase >= 0) {
	statsFile << std::setprecision(4) << samples[s].locus.means[0] << "\t" << std::setprecision(4) << samples[s].locus.means[1];
	++numNotNA;
      }
      else {
	statsFile << "NA\tNA";
      }
      if (s+1 < samples.size()) { statsFile << "\t";}
      
    }
    statsFile<< endl;
    for (auto s=0; s < samples.size(); s++) {    
      if (samples[s].MatchesRegion(regionChrom[locusIndex], regionStart[locusIndex], regionEnd[locusIndex])and
	  samples[s].locus.phase >= 0) {
	statsFile << std::setprecision(4) << samples[s].locus.vars[0] << "\t" << std::setprecision(4) << samples[s].locus.vars[1];
      }
      else {
	statsFile << "NA\tNA";
      }
      if (s+1 < samples.size()) { statsFile << "\t";}      
    }
    statsFile<< endl;
    if (numNotNA > 0.5 * samples.size()) {
      numOutput++;
    }

    //
    // Output IQR outliers.
    //

    for (auto s=0; s < samples.size(); s++) {    
      if (samples[s].MatchesRegion(regionChrom[locusIndex], regionStart[locusIndex], regionEnd[locusIndex])and
	  samples[s].locus.phase >= 0) {
	if (s > 0) { iqrFile << "\t"; }
	for (auto h=0; h < 2; h++ ) {
	  if (samples[s].locus.iqrExtremes[h].size() > 0) {
	    assert(samples[s].locus.iqr1s[h] != samples[s].locus.iqr2s[h]);
	    iqrFile << h << ":"<<samples[s].locus.iqr1s[h] << ":" << samples[s].locus.iqr2s[h] << ":";
	    for (auto iqi =0; iqi < samples[s].locus.iqrExtremes[h].size(); iqi++) {
	      iqrFile << samples[s].locus.iqrExtremes[h][iqi];
	      if (iqi + 1 < samples[s].locus.iqrExtremes[h].size()) {
		iqrFile << ",";
	      }
	    }
	  }
	  else {
	    iqrFile << h << ":NA:NA:";
	  }
	  iqrFile << ";";
	}
      }
      else {
	iqrFile << "NA\tNA";
      }
      if (s+1 < samples.size()) { iqrFile << "\t";}      
    }
    iqrFile<< endl;

    
  }
  
  return 0;
}
