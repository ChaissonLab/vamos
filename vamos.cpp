#include <iostream>
#include <fstream>
#include <iomanip>
#include <assert.h>
#include <vector>
using namespace std;


int LEFT=1;
int DIAG=2;
int DOWN=3;
int INIT=4;

int MATCH=1;
int MISMATCH=-1;
int GAP=-1;
using namespace std;

template <typename T>
class Matrix {
public:
  int nRow;
  int nCol;
  vector<T> data;
  T Get(int r, int c) {
    return data[r*nCol+c];
  }
  T& Set(int r, int c, T val) {
    data[r*nCol+c] = val;
    return data[r*nCol+c] ;    
  }
  void Init(int r, int c) {
    nRow=r;
    nCol=c;
    data.resize(nRow*nCol);
    fill(data.begin(), data.end(), 0);
  }
  T* Ptr(int r, int c) {
    return &data[r*nCol+c];
  }
  void Print(int w=2) {
    for (int i=0; i < nRow; i++) {
      for (int j=0; j < nCol; j++ ) {
	cout << std::setw(w) << Get(i,j);
	if (i + 1 < nCol) { cout << " "; }
      }
      cout << endl;
    }
  }
};


int TracebackMotif(vector<int> &scoreVntr,
		   Matrix<int> &scoreMat,
		   Matrix<int> &pathMat,
		   int refPos,
		   int &optMotifStart, int &optMotifEnd, int &optVntrPrev) {
  int optTraceScore=-1;
  int optTraceStartRow=-1;
  int optTraceEndRow=-1;
  vector<int> optPath;
  vector<int> vntrMatch, motifMatch, optScore;
  int maxRowScore=0;
  int maxRow=0;
  for (int motifEnd=scoreMat.nRow-1; motifEnd > 0; motifEnd--) {
    int score=scoreMat.Get(motifEnd, refPos);
    cout << "motif: " << motifEnd << " " << refPos << " " << score << endl;
    if (score > maxRowScore) {
      maxRow =  motifEnd;
      maxRowScore = score;
    }
  }
  if (maxRowScore == 0) {
    //
    // Nothing matches at this position.
    //
    optMotifStart=-1;
    optMotifEnd=-1;
    return 0;
  }
  

  int vPos=refPos;
  int motifPos=maxRow;
  int optMotifPrev=0;
  // Trace back from (motifPos, VPos) to a start.
  optPath.resize(0);
  vntrMatch.resize(0);
  motifMatch.resize(0);
  optScore.resize(0);
  optVntrPrev = refPos-1;
  while(pathMat.Get(motifPos, vPos) != INIT)  {
    assert(motifPos >= 0);
    assert(vPos >= 0);
    int arrow=pathMat.Get(motifPos, vPos);

    if (arrow == DIAG) {
      vntrMatch.push_back(vPos);
      motifMatch.push_back(motifPos);
      optScore.push_back(scoreMat.Get(motifPos, vPos));
      motifPos--;
      vPos--;
    }
    else if (arrow == LEFT) {
      vPos--;
      assert(vPos >= 0);
    }
    else {
      motifPos--;
    }
  }
  cout <<"traceback at " << vPos << "\t" << vntrMatch.size() << endl;
  int maxMotifScore=0;
  for (int i=0; i < optScore.size(); i++ ) {      
    //    cout << vntrMatch[i] << " " << motifMatch[i] << " " << optScore[i] << " " <<  scoreVntr[vntrMatch[i]] + maxRowScore - optScore[i] << endl;
    int motifScore=maxRowScore-optScore[i] + scoreVntr[vntrMatch[i]];
    if (motifScore > maxMotifScore) {
      maxMotifScore = motifScore;
      optMotifStart = motifMatch[i];
      optMotifEnd   = motifMatch[motifMatch.size()-1];
      optVntrPrev   = vntrMatch[i];
    }  
  }
  return maxMotifScore;
}
  
template <typename T>
void InitM(string &s1, string &s2, Matrix<T> &mat) {
  mat.Init(s1.size()+1, s2.size() + 1);
}

int SWAlign(string &s1, string &s2, Matrix<int> &scoreMat, Matrix<int> &pathMat) {
  InitM(s1, s2, scoreMat);
  InitM(s1, s2, pathMat);
  for (int i=0; i < s1.size(); i++) {
    pathMat.Set(i+1,0, INIT);
  }
  for (int j=0; j < s2.size(); j++) {
    pathMat.Set(0, j+1, INIT);
  }
  
  for (int i=0; i < s1.size(); i++) {
    int *left=scoreMat.Ptr(i+1,0);
    int *down=scoreMat.Ptr(i,1);
    int *diag=scoreMat.Ptr(i,0);
    int *score=scoreMat.Ptr(i+1,1);
    int *path=pathMat.Ptr(i+1,1);
    char *c1=&s1[i];
    char *c2=&s2[0];
    for (int j=0; j < s2.size(); j++, left++,down++,diag++, c2++, score++, path++) {
      int matchScore=0;
      if (*c1 == *c2) {
	matchScore=MATCH + *diag;
      }
      else {
	matchScore=MISMATCH + *diag;
      }
      int leftScore = *left + GAP;
      int downScore = *down + GAP;
      int maxScore = max(matchScore, max(leftScore, downScore));
      if (matchScore < 0) {
	*score=0;
	*path=INIT;
      }
      else {
	*score = maxScore;
	if (maxScore == matchScore) {
	  *path=DIAG;
	}
	else if (maxScore == leftScore) {
	  *path=LEFT;
	}
	else {
	  *path=DOWN;
	}
      }      
    }
  }
  return 0;
}


int main(int argc, const char* argv[]) {
  if (argc != 3) {
    cerr << "Usage: val motif_file vntr" << endl;
  }
  
  ifstream f1=ifstream(argv[1]);
  ifstream f2=ifstream(argv[2]);
  string vntr;
  vector<string> motifs;
  while (f1) {
    string motif;
    getline(f1,motif);
    if (motif != "") {
      motifs.push_back(motif);
    }
  }

  getline(f2,vntr);
  vector<Matrix<int> > scoreMat, pathMat;
  scoreMat.resize(motifs.size());
  pathMat.resize(motifs.size());
  for (int m=0; m < motifs.size(); m++) {
    SWAlign(motifs[m], vntr, scoreMat[m], pathMat[m]);
  }
  vector<int> scoreVntr(vntr.size()+1, 0), motifVntrIdx(vntr.size()+1, -1), optVntrPrevIdx(vntr.size()+1, -1);
  for (int vp=1; vp < vntr.size()+1; vp++) {
    int optMotifScore=0;
    int optMotifIndex=0;
    int optMotifStart=0, optMotifEnd=0;
    int optVntrPrev=vp-1;
    for (int m=0; m < motifs.size(); m++) {
      int curMotifStart, curMotifEnd;
      int motifScore = 0;
      int curVntrPrev;
      motifScore = TracebackMotif(scoreVntr, scoreMat[m], pathMat[m], vp, curMotifStart, curMotifEnd, curVntrPrev);
      cout << "vp: " << vp << " ms: " << motifScore << " opt " << optMotifScore << "\t" << curVntrPrev << endl;
      if (motifScore > optMotifScore) {
	optMotifScore = motifScore;
	optMotifIndex = m;
	optMotifStart = curMotifStart;
	optMotifEnd   = curMotifEnd;
	optVntrPrev   = curVntrPrev;
      }
    }
    scoreVntr[vp] = optMotifScore;
    motifVntrIdx[vp] = optMotifIndex;
    optVntrPrevIdx[vp] = optVntrPrev;
  }
  for (int i=0; i < scoreVntr.size(); i++) {
    cout << i << "\t" << scoreVntr[i] << "\t" << motifVntrIdx[i] << "\t" << optVntrPrevIdx[i] << endl;
  }
}

