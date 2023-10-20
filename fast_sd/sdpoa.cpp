#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <assert.h>  
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <stack>
#include <algorithm>

using namespace std;
int ReadSeq(ifstream &in, string &header, string &seq) {
  seq="";
  header="";
  string tmp;
  while (getline(in,tmp)) {
    if (tmp.size() == 0) { return 0;}
    if (tmp[0] == '>') {
      header=tmp;
    }
    else {
      seq += tmp;
    }
    if (in.peek() == '>') {
      return 1;
    }
  }
  if ( seq.size() == 0) {
    return 0;
  }
  else {
    return 1;
  }
}

// Split a string into tokens.
vector<string> split(string str, char delimiter) {
  vector<string> tokens;

  // Create a string stream from the string.
  stringstream ss(str);

  // Iterate over the string stream.
  string token;
  while (getline(ss, token, delimiter)) {
    tokens.push_back(token);
  }

  return tokens;
}

void PrintMat(vector<vector< int > > &mat, string seq="") {
  for (int i=0; i< mat.size(); i++ ) {
    if (seq != "" and i > 0) {
      cerr << seq[i-1];
    }
    if (i == 0 ) {
      cerr << " ";
    }
    cerr << std::setw(4) << i;
    for (int j=0; j < mat[i].size(); j++ ) {

      cerr << std::setw(4) << mat[i][j];
      if (j < mat[i].size() - 1 ) {
	cerr << " ";
      }
      else {
	cerr << endl;
      }
    }
  }
}
// A struct to store a node in the sequence graph.
void MakeMat(vector<vector< int > > &mat, int r, int c ) {
  mat.resize(r);
  for (int i=0; i < r; i++ ) {
    if (mat[i].size() > 0) {
      mat[i].clear();
    }
    mat[i].resize(c, 0);
    fill(mat[i].begin(), mat[i].end(), 0);
  }
}

class AlnMat {
public:
  vector<vector< int> > scoreMat, pathMat;
  
};

class Path {
public:
  vector<int> nodes;
  int Last() {return nodes[nodes.size()-1];}
};

// A struct to store a link in the sequence graph.
class Link {
public:
  string name;
  int src;
  int dest;
  int orientation;
};
const int leftArr=0;
const int diagArr=1;
const int upArr=2;
const int termArr=3;
const int wrapArr=4;
const int gap=-3;
const int match=4;
const int misMatch=-2;
const long int negInf=-99999999;

bool IsMatHop(int p) {
  return p>=6;
}

int GetPrev(int p, int &arrow, int &index) {
  if (p < 6 ) {
    return 0;
  }
  else {
    p-=6;
    if (p % 3 == 0 ) {
      arrow=leftArr;
      index=p/3;
    }
    else if (p % 3 == 1) {
      arrow=diagArr;
      index=p/3;
    }
    else {
      assert(p % 3 == 2);
      arrow=wrapArr;
      index=p/3;
    }
  }
  return p;
}

int SetArrow(int arrow, int index ) {
  int p;
  if (arrow == leftArr) {
    p = 6+(index*3);
  }
  else if (arrow == diagArr) {
    p = 7+(index*3);
  }
  else {
    p = 8 +(index*3);
  }
  return p;
}
    



    
class Node {
public:
  int index;
  string seq;
  bool aligned;
  int last;
  string starting, terminal;
  vector<int> in, out;
  //  vector<int> in, outNodes;
  Node() { index=-1; aligned=false; terminal="0"; starting="0";}
  vector<vector< int> > score, path;

  bool IsTerminal(int i, int j) {
    if (i > 0) return false;
    if (starting[j-1] == '1') {
      return true;
    }
    else {
      return false;
    }
  }
  

  void Alloc(string &query) {
    if (index == -1) {
      return;
    }
    int nRow = query.size()+1;
    int nCol = seq.size()+1;
    last=nCol-1;
    MakeMat(score, nRow, nCol);
    MakeMat(path, nRow, nCol);
  }

  int GetAlignScore(int pos, int row=-1) {
    //
    // By default run for whole grid.
    //
    if (row == -1) {
      row = score.size()-1;
    }
    assert(pos+1 < score[row].size());
    return score[row][pos+1];
  }

  void RealignRow(string &query, vector<Node> &nodes, int row, int maxWrapScore, int maxWrapNode) {
    int i=row;
    //
    // Handle first cell in the table, which may either be a restart or transition.
    //
    if (in.size() == 0) {
      if (score[i+1][0] < maxWrapScore) {
	score[i+1][0] = maxWrapScore;
	path[i+1][0]  = SetArrow(wrapArr, maxWrapNode);
      }
    }

    for (int j=0; j < seq.size(); j++) {
      if (j == 0 and in.size() > 0) {
	int maxPrev=negInf;
	int maxPrevNode=-1;
	for (int ini=0; ini < in.size(); ini++ ) {
	  int p=in[ini];
	  int inLeft = nodes[p].score[i+1][nodes[p].last];	  
	  if (inLeft > maxPrev) {
	    inLeft = maxPrev;
	    path[i+1][1] = SetArrow(leftArr, ini);
	  }
	  int inDiag;
	  if (query[i] == seq[0]) {
	    inDiag = nodes[p].score[i][nodes[p].last] + match;
	  }
	  else {
	    inDiag = nodes[p].score[i][nodes[p].last] + misMatch;
	  }
	  if (inDiag >= maxPrev) {
	    path[i+1][1] = SetArrow(diagArr, ini);
	    maxPrev = inDiag;
	  }
	  score[i+1][1] = maxPrev;
	  //
	  // Now check to see if this instead should start from a wraparound score
	  //
	  if (starting[0] == '1' and score[i+1][1] < maxWrapScore ) {
	    path[i+1][1] = SetArrow(wrapArr, maxWrapNode);
	    score[i+1][1] = maxWrapScore;
	  }	    
	}
      }
      else {
	int matScore, leftScore, upScore;
	if (query[i] == seq[j]) {
	  matScore = score[i][j] + match;
	}
	else {
	  matScore = score[i][j] + misMatch;
	}
	leftScore = score[i+1][j] + gap;
	upScore   = score[i][j] + gap;
	
	int maxRecScore = max(max(matScore, leftScore), upScore);

	/*	if (maxRecScore < maxWrapScore) {
	  path[i+1][j+1] = SetArrow(wrapArr, maxWrapNode);
	  score[i+1][j+1] = maxWrapScore;
	}
	else {*/
	  if (maxWrapScore > maxRecScore and starting[j] == '1' ) {
	    path[i+1][j+1] = SetArrow(wrapArr, maxWrapNode);
	    score[i+1][j+1] = maxWrapScore;
	  }
	  else {
	    if (matScore >= leftScore and matScore >= upScore) {
	      path[i+1][j+1]  = diagArr;
	      score[i+1][j+1] = matScore;
	    }
	    else if (leftScore >= matScore and leftScore >= upScore ) {
	      path[i+1][j+1] = leftArr;
	      score[i+1][j+1] = leftScore;
	    }
	    else {
	      path[i+1][j+1] = upArr;
	      score[i+1][j+1] = upScore;
	    }
	  }
	}
      //      }
    }
  }
  
  void Align(string &query, vector<Node> &nodes, int row=-1) {
    // rows are from the graph, cols are query.
    // An entire alignment to the graph is computed successively for each
    // position in the query.

    //
    // If there are no in-edges, initialize a boundary gap column.
    //
    if (row <= 0) {
      if (in.size() == 0) {
	for (int i=0; i < query.size(); i++) {
	  path[i+1][0] = upArr;
	  score[i+1][0] = score[i][0] + gap;
	}
      }

      //
      // Now initialize the first row columns 1-. The first column is special since
      // if this is a node with in-deg > 0 the gap score comes from the prev node.
      //
      if (in.size() > 0) {
	int inLeft, maxLeft=negInf, maxLeftI=-1;
	
      
	for (int ini=0; ini < in.size(); ini++ ) {
	  int p=in[ini];
	  assert(nodes[p].aligned);
	  inLeft = nodes[p].score[0][nodes[p].last] + gap;
	  if (inLeft > maxLeft) {
	    maxLeft = inLeft;
	    maxLeftI = ini;
	  }
	}
	path[0][1] = SetArrow(leftArr, maxLeftI);
	score[0][1] = maxLeft;
      }
      else {
	path[0][1] = leftArr;
	score[0][1] = score[0][0] + gap;
      }
      if (in.size() > 0 and (starting[0] == '1' and score[0][1] < 0)) {
	path[0][1] = termArr;
	score[0][1] = 0;
      }
      
      // Now initialize the rest of the first row.
      for (int i=1; i < seq.size(); i++ ) {
	if (starting[i] == '1') {
	  path[0][i+1] = termArr;
	  score[0][i+1] = 0;
	}
	else {
	  path[0][i+1] = leftArr;
	  score[0][i+1] = score[0][i] + gap;
	}
      }
#ifdef DEBUG
      cerr << "Init " << index << endl;
      for (int i=0; i < score[0].size(); i++ ) {
	cerr << score[0][i] << " ";
      }
      cerr << endl;
      for (int i=0; i < path[0].size(); i++ ) {
	cerr << path[0][i] << " ";
      }
      cerr << endl;
      cerr << endl;      
      
#endif      
      
    }


    //
    // Now process each character in the query.
    //
    int startRow, endRow;
    if (row < 0) {
      startRow = 0;
      endRow=query.size();
    }
    else {
      startRow = row;
      endRow = row+1;
    }
    int i=startRow;
    while (i < endRow) {
      //
      // For the first iterate overall in from other nodes.
      //
      for (int j=0; j < seq.size(); j++) {
      
	assert(seq.size() > 0);
	if (j == 0 and in.size() > 0) {
	  int inLeft, maxLeft=negInf, maxLeftI=-1;
	  int inDiag;
	  int maxScore=negInf;
	  for (int ini=0; ini < in.size(); ini++ ) {
	    int p=in[ini];
	    assert(nodes[p].aligned);
	    inLeft = nodes[p].score[i+1][nodes[p].last] + gap;
	    
	    if (query[i] == seq[0]) {
	      inDiag = nodes[p].score[i][nodes[p].last] + match;
	    }
	    else {
	      inDiag = nodes[p].score[i][nodes[p].last] + misMatch;
	    }
	    if (inLeft >= maxScore) {
	      path[i+1][1] = SetArrow(leftArr, ini);
	      maxScore = inLeft;
	    }
	    if (inDiag >= maxScore) {
	      path[i+1][1] = SetArrow(diagArr, ini);
	      maxScore = inDiag;
	    }
	  }
	  int upScore = score[i][1] + gap;
	  if (upScore > maxScore) {
	    path[i+1][1] = upArr;
	    maxScore = upScore;
	  }
	  score[i+1][1] = maxScore;	  
	}
	else {
	  int matScore, leftScore, upScore;
	  if (query[i] == seq[j]) {
	    matScore = score[i][j] + match;
	  }
	  else {
	    matScore = score[i][j] + misMatch;
	  }
	  leftScore = score[i+1][j] + gap;
	  upScore   = score[i][j] + gap;

	  int maxScore = max(max(matScore, leftScore), upScore);
	  if (matScore >= leftScore and matScore >= upScore) {
	    path[i+1][j+1]  = diagArr;
	    score[i+1][j+1] = matScore;
	  }
	  else if (leftScore >= matScore and leftScore >= upScore ) {
	    path[i+1][j+1] = leftArr;
	    score[i+1][j+1] = leftScore;
	  }
	  else {
	    path[i+1][j+1] = upArr;
	    score[i+1][j+1] = upScore;
	  }
	}
      }

      //
      // Process end of the while loop
      //

      i++;
	
    }

    aligned=true;
  }
};


void RecTopoSort(vector<Node> &nodes, vector<bool> &visited, stack<int> &nodeStack, int cur) {
  visited[cur] = true;
  for (int i=0; i < nodes[cur].out.size(); i++ ){
    int dest = nodes[cur].out[i];
    assert(nodes[dest].index >0);
    if (visited[dest] == false) {
      RecTopoSort(nodes, visited, nodeStack, dest);
    }
  }
  nodeStack.push(cur);
}

void TopologicalSort(vector<Node> &nodes, vector<int> &order) {
  stack<int> nodeStack;
  vector<bool> visited(nodes.size(), false);
  for (int i=0; i < nodes.size(); i++ ) {
    if (nodes[i].index >= 0 and visited[i]  == false) {
      RecTopoSort(nodes, visited, nodeStack, i);
    }
  }
  order.resize(nodeStack.size());
  for (int i=0; i < order.size(); i++ ) {
    order[i]=nodeStack.top();
    nodeStack.pop();
  }
}

void ReadGraph(ifstream &gfa_file, string &header, vector<Node> &nodes, vector<Link> &links, vector<Path> &paths) {
  string line;
  while (getline(gfa_file, line)) {
    // Split the line into tokens.
    vector<string> tokens = split(line, '\t');
    map<int,int>   indexToRank;
    // Check the type of the record.
    if (tokens[0] == "H") {
      header = line;
    }
    if (tokens[0] == "S") {
      // Create a node.
      Node node;
      node.index = stoi(tokens[1]);
      node.seq = tokens[2];
      if (node.index >= nodes.size()) {
	nodes.resize(node.index+1);
      }
      nodes[node.index] = node;
    } else if (tokens[0] == "L") {
      // Create a link.
      Link link;
      link.src = atoi(tokens[1].c_str());
      link.dest = atoi(tokens[3].c_str());
      links.push_back(link);
    }
    else if (tokens[0] == "P") {
      vector<string> pathNodes = split(tokens[2], ',');
      paths.push_back(Path());
      int p=paths.size()-1;
      for (int i=0; i < pathNodes.size(); i++ ) {
	paths[p].nodes.push_back(atoi(pathNodes[i].c_str()));
      }
    }
  }
  for (int l=0; l < links.size(); l++ ){    
    nodes[links[l].src].out.push_back(l);
    nodes[links[l].dest].in.push_back(l);
  }  
}

void MergeSimplePaths(vector<Node> &nodes, vector<Link> &links) {
  vector<char> removed(nodes.size(), 0);
  for (int cur=0; cur < nodes.size(); cur++) {
    while (nodes[cur].index >= 0 and nodes[cur].out.size() == 1 and nodes[links[nodes[cur].out[0]].dest].in.size() == 1) {
      int dest    = links[nodes[cur].out[0]].dest;
      //      cerr << "Merging " << cur << " and " << dest << endl;
      
      int curLink = nodes[cur].out[0];

      assert(dest >= 0);
      
      nodes[cur].seq = nodes[cur].seq + nodes[dest].seq;
      // Concatenate terminal/starting strings
      nodes[cur].terminal = nodes[cur].terminal + nodes[dest].terminal;
      nodes[cur].starting = nodes[cur].starting + nodes[dest].starting;

      int i=0;

      //
      // Propagate links to cur
      nodes[cur].out = nodes[dest].out;
      for (int l=0; l < nodes[dest].out.size(); l++ ){
	links[nodes[dest].out[l]].src = cur;
      }
      removed[dest] = true;
      nodes[dest].index=-1;
      
      /*
      //
      //
      for (int i=0; i < nodes[src].out.size(); i++ ) {
	// Relink src to dest
	if (links[nodes[src].out[i]].dest == nodes[n].index) {
	  links[nodes[src].out[i]].dest = dest;
	  break;
	}
      }
      */
    }     
  }
}


void PrintGraph(ostream &out, string &header, vector<Node> &nodes, vector<Link> &links, vector<Path> &paths ) {
  out << header << endl;
  for (int n=0; n < nodes.size(); n++ ) {
    if (nodes[n].index > 0) {
      out<< "S\t" << nodes[n].index << "_" << nodes[n].seq << "\t" << nodes[n].seq << endl;
    }
  }
  for (int l=0; l < links.size(); l++ ) {
    if (links[l].src != -1) {
      out << "L\t" << links[l].src << "_" << nodes[links[l].src].seq << "\t+\t" << links[l].dest << "_" << nodes[links[l].dest].seq  << "\t+\t0M"<< endl;
    }
  }  
}

void SetTerminal(vector<Node> &nodes, vector<Path> &paths) {
  for (int p=0; p < paths.size(); p++ ) {
    int n=paths[p].Last();
    nodes[n].terminal="1";
    n=paths[p].nodes[0];
    nodes[n].starting="1";
  }
  for (int n=0; n < nodes.size(); n++ ) {
    if (nodes[n].out.size() == 0) {
      nodes[n].terminal="1";
    }
    if (nodes[n].in.size() == 0) {
      nodes[n].starting="1";
    }
  }
}

class Align {
public:
  string read, str, graph;
};

void TraceAlign(vector<Node> &nodes, string &read, int maxNode, int maxPos, Align &aln) {


  int curI, curJ, curN;
  curI = read.size();
  curJ = maxPos;
  curN = maxNode;
  aln.read="";
  aln.str="";
  aln.graph="";
  while ((curI > 0 or curJ > 0) and  nodes[curN].IsTerminal(curI, curJ) == false ) {
    int arrow = nodes[curN].path[curI][curJ];
    if (IsMatHop(arrow)) {
      //      assert(curJ == 1 or curJ == nodes[curN].last or (arrow > 5 and );
      int hopArrow, hopNode;
      GetPrev(arrow, hopArrow, hopNode);
      if (hopArrow == diagArr) {
#ifdef DEBUG	
	cerr << "Mat hop match " << curN << "_" << nodes[curN].seq<< " " << curI << " " << curJ << " " << read[curI-1] << " " << nodes[curN].seq[curJ - 1]  << endl;
#endif

	aln.read+= read[curI-1];
	aln.graph+= nodes[curN].seq[curJ - 1];
	if (read[curI-1] != nodes[curN].seq[curJ - 1]) {
	  aln.str += ':';
	}
	else {
	  aln.str += '|';
	}
	
	curN = nodes[curN].in[hopNode];
	curI = curI-1;
	//
	// If this path can end at this position, do so now. We cannot store a terminal
	// arrow in the previous matrix because paths from other nodes could not be terminal
	// and reach there.
	//
	if (curI == 0 and nodes[curN].starting[0] == '1') {
	  break;
	}
	curJ = nodes[curN].last;
      }
      else if (hopArrow == leftArr)  {
	assert(hopArrow == leftArr);
#ifdef DEBUG
	cerr << "Mat hop gap graph " << " " << curN << "_" << nodes[curN].seq << " " << curI << " " << curJ << endl;
#endif
	aln.read += '-';
	aln.str  += ' ';
	assert(nodes[curN].seq[curJ- 1] != '\0');
	aln.graph += nodes[curN].seq[curJ-1];
	curN = nodes[curN].in[hopNode];
	curJ = nodes[curN].last;
      }
      else if (hopArrow == wrapArr) {
	//
	// This operation doesn't add any gaps or matches, it just jumps to a different
	// part of the graph.
	//
	curN = hopNode;
#ifdef DEBUG
	cerr << "Wrapping back to " << hopNode  << "_" << nodes[curN].seq << " " << nodes[curN].last << endl;
#endif
	curJ = nodes[curN].last;
      }
      else {
	cerr << "Should not hit this hop case." << endl;
	assert(0);
      }
    }
    else {
      if (arrow == diagArr ) {
#ifdef DEBUG
	cerr << "Same match " << curN << "_" << nodes[curN].seq << " " << curI << " " << curJ <<  " " << read[curI-1] << " " << nodes[curN].seq[curJ- 1] << endl;
#endif
	aln.read+= read[curI-1];
	aln.str += '|';
	aln.graph+= nodes[curN].seq[curJ-1];
	curI-=1;
	curJ-=1;
      }
      else if (arrow == leftArr ) {
#ifdef DEBUG
	cerr << "Same gap graph " << curN << "_" << nodes[curN].seq <<" " << curI << " " << curJ << endl;	
#endif
	aln.read += '-';
	aln.str += ' ';
	aln.graph += nodes[curN].seq[curJ - 1];
	curJ-=1;
      }
      else if (arrow == upArr) {
#ifdef DEBUG
	cerr << "Same gap read " << curN << "_" << nodes[curN].seq<< " " << curI << " " << curJ << endl;		
#endif
	aln.read += read[curI-1];
	aln.str += ' ';
	aln.graph+= '-';
	curI-=1;
      }
      else {
	assert(arrow == termArr);
	break;
      }
    }
  }
  std::reverse(aln.read.begin(), aln.read.end());
  std::reverse(aln.str.begin(), aln.str.end());
  std::reverse(aln.graph.begin(), aln.graph.end());
}
int FindMax(vector<Node> &nodes, int &maxTerminalNode, int &maxTerminalPos, int row=-1) {
  int maxTerminal=negInf;
  maxTerminalNode=-1;
  maxTerminalPos=-1;
  for (int i=0; i < nodes.size(); i++ ) {
    if (nodes[i].index != -1) {
      for (int p=0; p < nodes[i].terminal.size(); p++ ) {
	if (nodes[i].terminal[p] == '1') { 
	  int alnScore = nodes[i].GetAlignScore(p, row);
#ifdef DEBUG	 
	  cerr << "Node " << i << " pos " << p << " score " << alnScore << endl;
#endif	  
	  if (alnScore > maxTerminal) {
	    maxTerminal=alnScore;
	    maxTerminalNode= i;
	    maxTerminalPos = p + 1;
	  }
	}
      }
    }    
  }
  return maxTerminal;
}
int RunAlign(vector<Node> &nodes, string &read, vector<int> &order) {
  for (int n=0; n < nodes.size(); n++) {
    nodes[n].Alloc(read);
  }
  for (int j=0; j < read.size(); j++ ) {
    //
    // Run alignment in standard dynamic programming
    //
    for (int i = 0; i < order.size(); i++ ) {
      nodes[order[i]].Align(read, nodes, j);
    }
    int maxNode, maxPos, maxScore;
#ifdef DEBUG
    cerr << "iter " << j << endl;
#endif
    //
    // Rerun alignment for this row after copying over
    // maximum scoring end pos.
    //
    maxScore = FindMax(nodes, maxNode, maxPos, j+1);
    for (int i = 0; i < order.size(); i++ ) {
      nodes[order[i]].RealignRow(read, nodes, j, maxScore, maxNode);
    }
  }

#ifdef DEBUG
  for (int i=0; i < order.size(); i++) {
    int n=order[i];
    cerr << n << "_" << nodes[n].seq << endl;
    cerr << "Aligned " << nodes[n].index << " score:" << endl;
    PrintMat(nodes[n].score, read);
    cerr << "Aligned " << nodes[n].index << " path: " << endl;

    cerr << "         ";

    for (int i=0; i < nodes[n].seq.size(); i++ ) {
      cerr << setw(4) << nodes[n].seq[i];
      if (i < nodes[n].seq.size()-1) { cerr << " "; }
      else {
	cerr << endl;
      }
    }
    PrintMat(nodes[n].path, read);
  }
#endif

  return 0;
}

void LinksToOutEdges(vector<Node> &nodes, vector<Link> &links) {
  for (int n=0; n < nodes.size(); n++ ) {
    for (int i=0; i < nodes[n].in.size(); i++ ) {
      nodes[n].in[i] = links[nodes[n].in[i]].src;
    }
    for (int i=0; i < nodes[n].out.size(); i++ ) {
      nodes[n].out[i] = links[nodes[n].out[i]].dest;
    }
  }
}

int main(int argc, char* argv[]) {
  // Open the GFA file.
  if (argc < 3) {
    cerr << "Usage sdpoa graph.gfa query" << endl;
    exit(0);
  }
  ifstream gfa_file(argv[1]);
  ifstream faFile(argv[2]);
  // Declare a vector to store the nodes.
  vector<Node> nodes;

  // Declare a vector to store the links.
  vector<Link> links;
  vector<Path> paths;
  // Read the GFA file line by line.
  string line;
  string headerLine;

  ReadGraph(gfa_file, headerLine, nodes, links, paths);
  SetTerminal(nodes, paths);
  MergeSimplePaths(nodes, links);
  /*
  for (int i=0; i < nodes.size(); i++) {
    cerr << "node: " << i << endl;
    cerr << nodes[i].seq << endl;
  }
  */
  if (argc == 4) {
    ofstream gfaOut(argv[3]);
    PrintGraph(gfaOut, headerLine, nodes, links, paths);
    exit(0);
  }
  string read, name;
  //
  // Stop using links, those are lame. I'm sure this will cause a bug later on.
  //
  LinksToOutEdges(nodes, links);  
  vector<int> order;
  TopologicalSort(nodes, order);


  while (ReadSeq(faFile, name, read)) {
    RunAlign(nodes, read, order);
    int maxNode, maxPos;
    FindMax(nodes, maxNode, maxPos);
    Align aln;
    TraceAlign(nodes, read, maxNode, maxPos, aln);
    cerr << name << endl;
    cerr << aln.read << endl;
    cerr << aln.str << endl;
    cerr << aln.graph << endl;
  }
  return 0;
}

