#ifndef OPTION_H_
#define OPTION_H_

// using namespace std;
enum InputType {
  by_read,
  by_contig
};
class OPTION
{
public:
  int nproc;
  bool filterNoisy;
  double filterStrength;
  int penalty_indel;
  int penalty_mismatch;
  int match;
  int phaseFlank;
  string download;
  string referenceFilename;
  int maxCoverage;
  int maxLocusLength;
  int minAltCoverage;
  int minSNVCoverage;
  InputType inputType;
  bool hapChrX;
  int minChrY;
  bool refine;
  bool local;
  string reference;
  int oneOffset;
  int forcePhase;
  OPTION ()
  {
    nproc = 1;
    filterNoisy = false;
    filterStrength = 0.0;
    match=2;
    penalty_indel = 1;
    penalty_mismatch = 1;
    phaseFlank=15000;
    download="";
    inputType=by_read;
    maxCoverage=200;
    maxLocusLength=10000;
    minSNVCoverage = 6;
    minAltCoverage = 3;
    referenceFilename = "";
    minChrY = 0;
    reference="";
    hapChrX= false;
    refine=false;
    local=false;
    oneOffset=1;
    forcePhase=false;
  };

  ~OPTION () {};
};

#endif
