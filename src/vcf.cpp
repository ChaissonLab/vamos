#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <ostream>
#include "vcf.h"

void VcfWriter::writeHeader(ostream& out)
{
    out << "##fileformat=VCFv4.1\n";
    out << "##source=vamos" << version << '\n';
	for (int32_t i = 0; i < ncontigs; ++i) 
	{
		out << "##contig=<ID=" << target_names[i] << ",length=" << to_string(contigLengths[i]) << ">" << "\n";
	}
	out	<< "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">" << "\n"
		<< "##INFO=<ID=RU,Number=1,Type=String,Description=\"Comma separated motif sequences list in the reference orientation\">" << "\n"
		<< "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << "\n"
		<< "##INFO=<ID=ALTANNO_H1,Number=1,Type=String,Description=\"\"Motif representation for the h1 alternate allele>" << "\n"
		<< "##INFO=<ID=ALTANNO_H2,Number=1,Type=String,Description=\"\"Motif representation for the h2 alternate allele>" << "\n"

		<< "##FILTER=<ID=PASS,Description=\"All filters passed\">" << "\n"
		<< "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << "\n"
		<< "##ALT=<ID=VNTR,Description=\"Allele comprised of VNTR repeat units\">" << "\n";

    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
    return;
}


void VcfWriter::writeBody(vector<VNTR *> &vntrs, ostream& out)
{
	for (auto &it : vntrs)
	{
		if (it->nreads == 0) continue;
		string motif_list, motif_anno_h1, motif_anno_h2, GT; 
		for (auto &piece : it->motifs) 
		{ 
			motif_list += piece.seq + ',';
		}
	    if (!motif_list.empty()) motif_list.pop_back();

	    it->commaSeparatedMotifAnnoForConsensus(1, motif_anno_h1);
	    it->commaSeparatedMotifAnnoForConsensus(0, motif_anno_h2);

	    if (motif_anno_h1 == motif_anno_h2) GT = "1/1";
	    else GT = "1/2";

		out << it->chr << "\t"
			<< to_string(it->ref_start) << "\t"
			<< ".\t"
			<< "N\t"
			<< "<VNTR>\t"
			<< ".\t"
			<< "PASS\t"
			<< "END=" + to_string(it->ref_end) + ";" 
			<< "RU=" + motif_list + ";"
			<< "SVTYPE=VNTR;"
			<< "ALTANNO_H1=" + motif_anno_h1 + ";";

		if (GT == "1/2")
			out << "ALTANNO_H2=" + motif_anno_h2 + ";\t";

		out	<< "PASS\t"
			<< "GT\t"
			<< GT + "\n";
	}
	return;
}

