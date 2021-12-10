#ifndef VCF_H_
#define VCF_H_
#include "htslib/htslib/sam.h"
#include <string>
#include <vector>

class VcfWriter
{
public:
	string sampleName;
	string version;
	char * target_names;
	uint32_t * contigLengths;
	int32_t ncontigs;

	VcfWriter (string input_bam_file, string Version, string SampleName) : version(Version), sampleName(SampleName)
	{
	    bamFile * fp_in = bam_open(input_bam_file, "r"); //open bam file
	    bam_hdr_t * bamHdr = bam_hdr_read(fp_in); //read header

	    int32_t ncontigs = bamHdr->n_targets;
	    target_names = (char *) malloc(ncontigs);
	    contigLengths = (uint32_t *) malloc(ncontigs);

	    for (int32_t i = 0; i < ; ++i)
	    {
	    	contigLengths[i] = bamHdr->target_len[i];
	    	target_names[i] = bamHdr->target_name[i];
	    }

	    bam_hdr_destroy(bamHdr);
	    bam_close(fp_in); 
	};

	~VcfWriter () {};

	void free ()
	{
		free (target_names);
		free (contigLengths)
	}

    void writeHeader(ostream& out);
    void writeBody(ostream& out);
};

void VcfWriter::writeHeader(ostream& out)
{
    out << "##fileformat=VCFv4.1\n";
        << "##source=vamos" << version << '\n'
	for (int32_t i = 0; i < ncontigs; ++i) 
	{
		out << "##contig=<ID=" << target_names[i] << ",length=" << contigLengths[i] << ">" << "\n";
	}
	out	<< "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">" << "\n"
		<< "##INFO=<ID=RU,Number=1,Type=String,Description=\"Comma separated motif sequences list in the reference orientation\">" << "\n"
		<< "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << "\n"
		<< "##INFO=<ID=ALTREP_H1,Number=1,Type=String,Description=\"\"Motif representation for the h1 alternate allele>" << "\n"
		<< "##INFO=<ID=ALTREP_H2,Number=1,Type=String,Description=\"\"Motif representation for the h2 alternate allele>" << "\n"

		<< "##FILTER=<ID=PASS,Description=\"All filters passed\">" << "\n"
		<< "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << "\n"
		<< "##ALT=<ID=VNTR,Description=\"Allele comprised of VNTR repeat units\">" << "\n";

    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
    return;
}


void VcfWriter::writeBody(const vector<VNTR> &vntrs, ostream& out)
{
	for (auto &it : vntrs)
	{
		string motif_list, motif_rep_h1, motif_rep_h2, GT; 
		for (const auto &piece : motifs) 
		{ 
			motif_list += piece + ',';
		}
	    if (!motif_list.empty()) motif_list.pop_back();

	    it.commaSeparatedMotifRepForConsensus(1, motif_rep_h1);
	    it.commaSeparatedMotifRepForConsensus(0, motif_rep_h2);

	    if (motif_rep_h1 == motif_rep_h2)
	    	GT = "1/1"
	    else
	    	GT = "1/2"

		out << it.chr[0] << "\t"
			<< to_string(it.ref_start) << "\t"
			<< ".\t"
			<< "N\t"
			<< "<VNTR>\t"
			<< ".\t"
			<< "PASS\t"
			<< "END=" + to_string(it.ref_end) + ";" 
			<< "RU=" + motif_list + ";"
			<< "SVTYPE=VNTR;"
			<< "ALTREP_H1=" + motif_rep_h1 + ";\t";

		if (GT == "1/2")
			out << "ALTREP_H2=" + motif_rep_h2 + ";\t";

		out	<< "PASS\t"
			<< "GT\t"
			<< GT + "\n";
	}
	return;
}

#endif