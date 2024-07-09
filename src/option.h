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
  double penalty_indel;
  double penalty_mismatch;
  double accuracy;
  int phaseFlank;
  string download;
  int maxCoverage;
  int maxLocusSize;
  InputType inputType;
  OPTION ()
  {
    nproc = 1;
    filterNoisy = false;
    filterStrength = 0.0;
    penalty_indel = 1.0;
    penalty_mismatch = 1.0;
    accuracy=0.98;
    phaseFlank=15000;
    download="";
    inputType=by_read;
    maxCoverage=200;
    maxLocusSize=20000;
  };

  ~OPTION () {};
};

#endif
