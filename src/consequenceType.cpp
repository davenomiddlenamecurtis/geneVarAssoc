#include "consequenceType.hpp"
#include <string.h>

#ifdef USEOLDWEIGHTS
consequenceReport consequence[]={
	{ NULL_CONSEQUENCE,"NULL_CONSEQUENCE",1.0 },
	{ INTERGENIC,"INTERGENIC",1.0 },
	{ DOWNSTREAM,"DOWNSTREAM",1.0 },
	{ INTRONIC,"INTRONIC",3.0 },
	{ THREEPRIME_UTR,"3PRIME_UTR",3.0 },
	{ SYNONYMOUS_CODING,"SYNONYMOUS_CODING",5.0 },
	{ UPSTREAM,"UPSTREAM",5.0 },
	{ FIVEPRIME_UTR,"5PRIME_UTR",5.0 },
	{ SPLICE_SITE,"SPLICE_SITE",5.0 },
	{ STOP_LOST,"STOP_LOST",5.0 },
	{ NON_SYNONYMOUS_CODING,"NON_SYNONYMOUS_CODING",10.0 },
	{ CODINGINDEL,"CODINGINDEL",15.0 },
	{ FRAMESHIFT_CODING,"FRAMESHIFT_CODING",20.0 },
	{ STOP_GAINED,"STOP_GAINED",20.0 },
	{ KOZAK,"KOZAK",10.0 }
};
#else
consequenceReport consequence[]={
	{ NULL_CONSEQUENCE,"NULL_CONSEQUENCE",1.0 },
	{ INTERGENIC,"INTERGENIC",1.0 },
	{ DOWNSTREAM,"DOWNSTREAM",1.0 },
	{ INTRONIC,"INTRONIC",3.0 },
	{ THREEPRIME_UTR,"THREEPRIME_UTR",5.0 },
	{ UPSTREAM,"UPSTREAM",5.0 },
	{ SYNONYMOUS_CODING,"SYNONYMOUS_CODING",3.0 },
	{ FIVEPRIME_UTR,"FIVEPRIME_UTR",5.0 },
	{ CODINGINDEL,"CODINGINDEL",15.0 },
	{ SPLICE_SITE,"SPLICE_SITE",5.0 },
	{ STOP_LOST,"STOP_LOST",5.0 },
	{ NON_SYNONYMOUS_CODING,"NON_SYNONYMOUS_CODING",10.0 },
	{ FRAMESHIFT_CODING,"FRAMESHIFT_CODING",20.0 },
	{ ESSENTIAL_SPLICE_SITE,"ESSENTIAL_SPLICE_SITE",10.0 },
	{ STOP_GAINED,"STOP_GAINED",20.0 },
	{ KOZAK,"KOZAK",10.0 }
};

#endif
consequenceReport e_consequence[]={
	{ E_NULL_CONSEQUENCE,"NULL_CONSEQUENCE",1.0 },
	{ E_INTERGENIC_VARIANT,"intergenic_variant",1.0 },
	{ E_FEATURE_TRUNCATION,"feature_truncation",3.0 },
	{ E_REGULATORY_REGION_VARIANT,"regulatory_region_variant",3.0 },
	{ E_FEATURE_ELONGATION,"feature_elongation",3.0 },
	{ E_REGULATORY_REGION_AMPLIFICATION,"regulatory_region_amplification",3.0 },
	{ E_REGULATORY_REGION_ABLATION,"regulatory_region_ablation",3.0 },
	{ E_TF_BINDING_SITE_VARIANT,"TF_binding_site_variant",3.0 },
	{ E_TFBS_AMPLIFICATION,"TFBS_amplification",3.0 },
	{ E_TFBS_ABLATION,"TFBS_ablation",3.0 },
	{ E_DOWNSTREAM_GENE_VARIANT,"downstream_gene_variant",3.0 },
	{ E_UPSTREAM_GENE_VARIANT,"upstream_gene_variant",3.0 },
	{ E_NON_CODING_TRANSCRIPT_VARIANT,"non_coding_transcript_variant",3.0 },
	{ E_NMD_TRANSCRIPT_VARIANT,"NMD_transcript_variant",3.0 },
	{ E_INTRON_VARIANT,"intron_variant",3.0 },
	{ E_NON_CODING_TRANSCRIPT_EXON_VARIANT,"non_coding_transcript_exon_variant",3.0 },
	{ E_3_PRIME_UTR_VARIANT,"3_prime_UTR_variant",5.0 },
	{ E_5_PRIME_UTR_VARIANT,"5_prime_UTR_variant",5.0 },
	{ E_MATURE_MIRNA_VARIANT,"mature_miRNA_variant",5.0 },
	{ E_CODING_SEQUENCE_VARIANT,"coding_sequence_variant",5.0 },
	{ E_SYNONYMOUS_VARIANT,"synonymous_variant",5.0 },
	{ E_STOP_RETAINED_VARIANT,"stop_retained_variant",5.0 },
	{ E_INCOMPLETE_TERMINAL_CODON_VARIANT,"incomplete_terminal_codon_variant",5.0 },
	{ E_SPLICE_REGION_VARIANT,"splice_region_variant",5.0 },
	{ E_PROTEIN_ALTERING_VARIANT,"protein_altering_variant",10.0 },
	{ E_MISSENSE_VARIANT,"missense_variant",10.0 },
	{ E_INFRAME_DELETION,"inframe_deletion",15.0 },
	{ E_INFRAME_INSERTION,"inframe_insertion",15.0 },
	{ E_TRANSCRIPT_AMPLIFICATION,"transcript_amplification",15.0 },
	{ E_START_LOST,"start_lost",15.0 },
	{ E_STOP_LOST,"stop_lost",15.0 },
	{ E_FRAMESHIFT_VARIANT,"frameshift_variant",20.0 },
	{ E_STOP_GAINED,"stop_gained",20.0 },
	{ E_SPLICE_DONOR_VARIANT,"splice_donor_variant",20.0 },
	{ E_SPLICE_ACCEPTOR_VARIANT,"splice_acceptor_variant",20.0 },
	{ E_TRANSCRIPT_ABLATION,"transcript_ablation",20.0 }
};

const char* getConsequenceString(int t,int use_ensembl)
{
	if (use_ensembl)
		for(int i=0;i<E_NCONSEQUENCETYPES;++i)
			if(t==e_consequence[i].t)
				return e_consequence[i].str;
	else
		for (int i=0;i<NCONSEQUENCETYPES;++i)
			if (t==consequence[i].t)
				return consequence[i].str;
	return "";
}

int getConsequenceType(const char *name,int use_ensembl)
{
	if (use_ensembl)
		for (int i=0;i<E_NCONSEQUENCETYPES;++i)
			if (!strncmp(e_consequence[i].str,name,strlen(e_consequence[i].str)))
				return e_consequence[i].t;
			else
				for (int i=0;i<NCONSEQUENCETYPES;++i)
					if (!strcmp(consequence[i].str,name))
						return consequence[i].t;
	return -1;
}

