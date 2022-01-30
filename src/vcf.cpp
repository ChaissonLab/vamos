#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <istream>
#include <fstream>
#include <ostream>
#include "assert.h"
#include "vcf.h"

void VcfWriter::init (char * input_bam_file, char * Version, char * SampleName)
{
	// version = (char *) malloc(strlen(Version) + 1);
	// strcpy(version, Version);
	// sampleName = (char *) malloc(strlen(SampleName) + 1);
	// strcpy(sampleName, SampleName);
	version = Version;
	sampleName = SampleName;

    samFile * fp_in = hts_open(input_bam_file, "r"); //open bam file
    bam_hdr_t * bamHdr = sam_hdr_read(fp_in); //read header
    ncontigs = bamHdr->n_targets;

    for (int32_t i = 0; i < ncontigs; ++i)
    {
    	contigLengths.push_back(bamHdr->target_len[i]);
    	target_names.push_back(string(bamHdr->target_name[i]));
    }

    assert((uint32_t) ncontigs == target_names.size() and (uint32_t) ncontigs == contigLengths.size());
    bam_hdr_destroy(bamHdr);
    sam_close(fp_in); 
    set = 1;
};


void VcfWriter::writeHeader(ofstream &out)
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
		<< "##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"\"Length of the vntr sequence>" << "\n"

		<< "##FILTER=<ID=PASS,Description=\"All filters passed\">" << "\n"
		<< "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << "\n"
		<< "##ALT=<ID=VNTR,Description=\"Allele comprised of VNTR repeat units\">" << "\n";

    out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName << "\n";
    cerr << "finish writing the header" << endl;
    return;
}

void writeSingleBody(VNTR * it, ofstream &out)
{
	if (it->skip) return;
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

	out << it->chr << "\t";
	out	<< to_string(it->ref_start) << "\t";
	out	<< ".\t";
	out	<< "N\t";
	out	<< "<VNTR>\t";
	out	<< ".\t";
	out	<< "PASS\t";
	out	<< "END=" + to_string(it->ref_end) + ";";
	out	<< "RU=" + motif_list + ";";
	out	<< "SVTYPE=VNTR;";
	out	<< "ALTANNO_H1=" + motif_anno_h1 + ";";

	if (GT == "1/2")
		out << "ALTANNO_H2=" + motif_anno_h2 + ";\t";

	out	<< "LEN=" + to_string(it->cur_len) + ";";
	out	<< "PASS\t";
	out	<< "GT\t";
	out << GT + "\n";
	return;
}


void VcfWriter::writeBody(vector<VNTR *> &vntrs, ofstream& out, int tid, int nproc)
{
	if (nproc > 1)
	{
		for (size_t i = tid; i < vntrs.size(); i += nproc)
			writeSingleBody(vntrs[i], out);
	}
	else
	{
		for (auto &it : vntrs) 
			writeSingleBody(it, out);
	}

	return;
}

void writeSingleNullAnno(VNTR * it, ofstream &out_nullAnno)
{
	if (!it->nullAnno) return;

	bool na;
	for (int t = 0; t < it->nreads; ++t)
	{
		na = it->nullAnnos[t];
		if (na)
		{
			out_nullAnno << it->region << "\t";
			out_nullAnno.write(it->reads[t]->seq, it->reads[t]->len);
			out_nullAnno << "\t"; 
			size_t t = 0;
			for (auto &mt : it->motifs)
			{
				if (t < it->motifs.size() - 1)
					out_nullAnno << mt.seq << ",";
				else
					out_nullAnno << mt.seq << "\n";
				t++;
			}			
		}
	}

	return;
}

void VcfWriter::writeNullAnno(vector<VNTR *> &vntrs, ofstream &out_nullAnno, int tid, int nproc)
{
	if (nproc > 1)
	{
		for (size_t i = tid; i < vntrs.size(); i += nproc)
			writeSingleNullAnno(vntrs[i], out_nullAnno);
	}
	else
	{
		for (auto &it : vntrs) 
			writeSingleNullAnno(it, out_nullAnno);
	}

	return;
}


