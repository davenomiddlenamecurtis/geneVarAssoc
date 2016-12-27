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
	E_INTERGENIC,
	E_NC_TRANSCRIPT,
	E_DOWNSTREAM500B,
	E_DOWNSTREAM5KB,
	E_INTRONIC,
	E_THREEPRIME_UTR,
	E_TRANSCRIPTIONFACTOR,
	E_REGULATORY,
	E_UPSTREAM2K,
	E_UPSTREAM5K,
	E_NMD_TRANSCRIPT,
	E_NON_CODING_TRANSCRIPT_EXON_VARIANT,
	E_SYNONYMOUS_CODING,
	E_SYNONYMOUS_VARIANT,
	E_STOP_RETAINED,
	E_CODING_UNKNOWN,
	E_FIVEPRIME_UTR,
	E_CODINGGAIN,
	E_INFRAME_INSERTION,
	E_CODINGLOSS,
	E_INFRAME_DELETION,
	E_SPLICE_DONOR,
	E_SPLICE_SITE,
	E_STOP_LOST,
	E_NON_SYNONYMOUS_CODING,
	E_MISSENSE_VARIANT,
	E_CODON_CHANGE,
	E_NON_TERMINAL_CODON,
	E_MIRNA_VARIANT,
	E_FRAMESHIFT_CODING,
	E_START_LOST,
	E_ESSENTIAL_SPLICE_SITE,
	E_TRANSCRIPT_CHANGE,
	E_STOP_GAINED,
	E_NCONSEQUENCETYPES
};

struct consequenceReport_t {
	int t;
	char *str;
	float weight;
};

typedef struct consequenceReport_t consequenceReport;

extern consequenceReport consequence[],e_consequence[];
const char* consequenceString(int t,int use_ensembl);

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