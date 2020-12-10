#if 0
Copyright 2018 David Curtis

This file is part of the geneVarAssoc package.

geneVarAssoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

geneVarAssoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with geneVarAssoc. If not, see <http://www.gnu.org/licenses/>.
#endif

#ifndef CONSEQUENCETYPEHPP
#define CONSEQUENCETYPEHPP

enum consequenceType {
	NULL_CONSEQUENCE=0,
	INTERGENIC,
	DOWNSTREAM,
	INTRONIC,
	THREEPRIME_UTR,
	UPSTREAM,
	SYNONYMOUS_CODING,
	FIVEPRIME_UTR,
	CODINGINDEL,
	SPLICE_SITE,
	STOP_LOST,
	NON_SYNONYMOUS_CODING,
	FRAMESHIFT_CODING,
	ESSENTIAL_SPLICE_SITE,
	STOP_GAINED,
	KOZAK,
	NCONSEQUENCETYPES
};

enum e_consequenceType {
	E_NULL_CONSEQUENCE=0,
	E_INTERGENIC_VARIANT,
	E_FEATURE_TRUNCATION,
	E_REGULATORY_REGION_VARIANT,
	E_FEATURE_ELONGATION,
	E_REGULATORY_REGION_AMPLIFICATION,
	E_REGULATORY_REGION_ABLATION,
	E_TF_BINDING_SITE_VARIANT,
	E_TFBS_AMPLIFICATION,
	E_TFBS_ABLATION,
	E_DOWNSTREAM_GENE_VARIANT,
	E_UPSTREAM_GENE_VARIANT,
	E_NON_CODING_TRANSCRIPT_VARIANT,
	E_NMD_TRANSCRIPT_VARIANT,
	E_INTRON_VARIANT,
	E_NON_CODING_TRANSCRIPT_EXON_VARIANT,
	E_3_PRIME_UTR_VARIANT,
	E_5_PRIME_UTR_VARIANT,
	E_MATURE_MIRNA_VARIANT,
	E_CODING_SEQUENCE_VARIANT,
	E_SYNONYMOUS_VARIANT,
	E_STOP_RETAINED_VARIANT,
	E_INCOMPLETE_TERMINAL_CODON_VARIANT,
	E_SPLICE_REGION_VARIANT,
	E_PROTEIN_ALTERING_VARIANT,
	E_MISSENSE_VARIANT,
	E_INFRAME_DELETION,
	E_INFRAME_INSERTION,
	E_TRANSCRIPT_AMPLIFICATION,
	E_START_LOST,
	E_STOP_LOST,
	E_FRAMESHIFT_VARIANT,
	E_STOP_GAINED,
	E_SPLICE_DONOR_VARIANT,
	E_SPLICE_ACCEPTOR_VARIANT,
	E_TRANSCRIPT_ABLATION,
	E_NCONSEQUENCETYPES
};

struct consequenceReport_t {
	int t;
	char *str;
	double weight;
};

typedef struct consequenceReport_t consequenceReport;

extern consequenceReport consequence[],e_consequence[];
const char* getConsequenceString(int t,int use_ensembl);
int getConsequenceType(const char *name,int use_ensembl);

// I’ve talked with Andy about the splicing site definition. 
// Reading literature, we decide that the donor splicing site (moving from exon to intron) we have to consider 5bp in the exon and 6 in the intron. 
// For the acceptor splicing site (moving from intron to exon) we have to consider 20bp in the intron and 3 in the exon. I’ll ask knome their definition of splicing site.
// NSB is Number of Splicesite Bases

// I will add a new type called essential splice site which disrupts either the 5' (donor) GT or 3' (acceptor) AG. I.e. 2 bases in the intron
// "Most of the splice site mutations that lead to human disease involve the invariant GT and AG dinucleotides in the 5’ and 3’ splice sites. " http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2855871/

#define NSBDONOREXON 5
#define NSBDONORINTRON 6
#define NSBACCEPTOREXON 3
#define NSBACCEPTORINTRON 20
enum spliceSiteType {DONOR,ACCEPTOR};

#endif