#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include "vntr.h"
#include "naive_anno.cpp"
#include "seqan/sequence.h"
#include "seqan/align.h"
#include "seqan/score.h"

void VNTR::motifAnnoForOneVNTR () 
{
	annos.resize(reads.size());
	for (size_t i = 0; i < reads.size(); ++i)
	{
		/* 
			apply Bida's code here
			input:   string : reads[i].seq; vector<MOTIF> : motifs
			output:  vector<vector<int>> annos[i]
		*/
		anno(annos[i], motifs, reads[i]->seq);
		cerr << "finish for reads: " << i << endl;
	}
	return;
}

void VNTR::annoTostring (vector<int> &anno, string &annostr)
{
	if (anno.size() == 0) return;

    annostr.clear();
	for (auto &num : anno)
	{
		// TODO: needs to remove '-' from encoding
		if (num >= 0 && num <= 222) {
			string tmp_s;
			tmp_s.assign(1, (char) (num + 33));
			annostr += tmp_s;			
		}

	}
    cerr << "string: " << annostr;
    cerr << endl;
	return;
}

// ----------------------------------------------------------------------------
// Helper Function getMaxIndex()
// ----------------------------------------------------------------------------

/*!
 * @fn ProfileChar#getMaxIndex
 * @brief Return number of dominating entry in ProfileChar.
 *
 * @signature TSize getMaxIndex(c);
 *
 * @param[in] c     ProfileChar to query for its dominating entry.
 * @return    TSize index (with the @link FiniteOrderedAlphabetConcept#ordValue @endlink) of the dominating character
 *                  in <tt>c</tt>
 */
template <typename TSourceValue, typename TSourceCount, typename TSourceSpec>
pair <typename Size<seqan::ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const>::Type, typename Size<seqan::ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const>::Type>
getFirst_Second_MaxIndex(seqan::ProfileChar<TSourceValue, TSourceCount, TSourceSpec> const & source)
{
    typedef seqan::ProfileChar<TSourceValue, TSourceCount, TSourceSpec> TProfileChar;
    typedef typename Size<TProfileChar>::Type TSize;


    TSize secmaxIndex = 0;
    TSourceCount secmaxCount = source.count[0];

    TSize MaxIndex = 0;
    TSize scdMaxIndex = 0;
    TSourceCount MaxCount = source.count[0];
    TSourceCount scdMaxCount = source.count[0];

    if (seqan::ValueSize<TProfileChar>::VALUE <= 1)
    	return make_pair(MaxIndex, scdMaxIndex);

   if(source.count[0] < source.count[1])
   { 
      MaxCount = source.count[1];
      scdMaxCount = source.count[0];
      MaxIndex = 1;
      scdMaxIndex = 0;
   }
   else
   { 
      MaxCount = source.count[0];
      scdMaxCount = source.count[1];
      MaxIndex = 0;
      scdMaxIndex = 1;
   }

   for (TSize i = 2; i < seqan::ValueSize<TProfileChar>::VALUE; i++) 
   {
      if (source.count[i] > MaxCount) 
      {
         scdMaxCount = MaxCount;
         scdMaxIndex = MaxIndex;
         MaxCount = source.count[i];
         MaxIndex = i;
      }
      else if (source.count[i] > scdMaxCount && source.count[i] != MaxCount) 
      {
         scdMaxCount = source.count[i];
         scdMaxIndex = i;
      }
   }
   
   return make_pair(MaxIndex, scdMaxIndex);
}

void VNTR::concensusMotifAnnoForOneVNTR ()
{
	/* based on vector<vector<int>> annos, get a concensus representation of the current locus 
	   input:   vector<vector<int>> annos
	   output:  vector<int> concensus_h1, vector<int> concensus_h2
	*/
    seqan::StringSet<seqan::String<char>> annoSet;
    string tmp;
    for (auto &it : annos)
    {
    	annoTostring (it, tmp);
    	seqan::appendValue(annoSet, tmp);
    }

    seqan::Align<seqan::String<char>> align(annoSet); // Initialize the Align object using a StringSet.

	seqan::Score<int> scoreScheme(1, -1, -1, -1); 
	int score = seqan::globalAlignment(align, scoreScheme);  // Compute a global alingment using the Align object.
    // int score = seqan::globalAlignment(align, seqan::EditDistanceScore());  // Compute a global alingment using the Align object.

    cerr << "score = " << score << endl;
    cerr << "align\n" << align << endl;

    /*
     create the profile string
     profile: row: alignment position
     		  col: reads
    */
    seqan::String<seqan::ProfileChar<char>> profile; 
    seqan::resize(profile, seqan::length(row(align, 0))); 
    for (int rowNo = 0; rowNo < annos.size(); ++rowNo) 
    {
        for (uint32_t i = 0; i < seqan::length(row(align, rowNo)); ++i) 
            profile[i].count[seqan::getValue(row(align, rowNo), i) - 33] += 1;
    }

    seqan::String<string> consensus_one, consensus_two_f, consensus_two_s;
    for (int i = 0; i < seqan::length(profile); ++i)
    {
        int [idx_f, idx_s] = getFirst_Second_MaxIndex(profile[i]);
        seqan::appendValue(consensus, to_string(idx_f));
    }

    cerr << "consensus sequence is\n";
    for (auto &it : consensus)
		cerr << it << "->";
	cerr << endl;

	return;
}

void VNTR::commaSeparatedMotifAnnoForConsensus (bool h1, string &motif_rep)
{
	if (h1)
		for (auto &it : concensus_h1) { motif_rep += "VNTR_" + to_string(it) + ",";}
	else
		for (auto &it : concensus_h2) { motif_rep += "VNTR_" + to_string(it) + ",";}
	if (!motif_rep.empty()) motif_rep.pop_back();
}
