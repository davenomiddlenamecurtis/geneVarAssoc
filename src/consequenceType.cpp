#include "consequenceType.hpp"

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
	{ THREEPRIME_UTR,"3PRIME_UTR",5.0 },
	{ SYNONYMOUS_CODING,"SYNONYMOUS_CODING",3.0 },
	{ UPSTREAM,"UPSTREAM",5.0 },
	{ FIVEPRIME_UTR,"5PRIME_UTR",5.0 },
	{ SPLICE_SITE,"SPLICE_SITE",5.0 },
	{ STOP_LOST,"STOP_LOST",5.0 },
	{ NON_SYNONYMOUS_CODING,"NON_SYNONYMOUS_CODING",10.0 },
	{ CODINGINDEL,"CODINGINDEL",15.0 },
	{ FRAMESHIFT_CODING,"FRAMESHIFT_CODING",20.0 },
	{ ESSENTIAL_SPLICE_SITE,"ESSENTIAL_SPLICE_SITE",10.0 },
	{ STOP_GAINED,"STOP_GAINED",20.0 },
	{ KOZAK,"KOZAK",10.0 }
};

#endif
consequenceReport e_consequence[]={
	{		E_NULL_CONSEQUENCE,"NULL_CONSEQUENCE",1.0 },
	{ E_INTERGENIC,"intergenic_variant",1.0 },
	{ E_NC_TRANSCRIPT,"nc_transcript_variant",1.0 },
	{ E_DOWNSTREAM500B,"500B_downstream_variant",1.0 },
	{ E_DOWNSTREAM5KB,"500B_downstream_variant",1.0 },
	{ E_INTRONIC,"intron_variant",3.0 },
		{ E_THREEPRIME_UTR,"3_prime_UTR_variant",5.0 },
		{ E_TRANSCRIPTIONFACTOR,"TF_binding_site_variant",5.0 },
		{ E_REGULATORY,"regulatory_region_variant",5.0 },
		{ E_UPSTREAM2K,"2KB_upstream_variant",5.0 },
		{ E_UPSTREAM5K,"5KB_upstream_variant",5.0 },
		{ E_NMD_TRANSCRIPT,"NMD_transcript_variant",3.0 },
		{ E_NON_CODING_TRANSCRIPT_EXON_VARIANT,"non_coding_transcript_exon_variant",3.0 },
		{ E_SYNONYMOUS_CODING,"synonymous_codon",3.0 },
		{ E_SYNONYMOUS_VARIANT,"synonymous_variant",3.0 },
		{ E_STOP_RETAINED,"stop_retained_variant",3.0 },
		{ E_CODING_UNKNOWN,"coding_sequence_variant",3.0 },
		{ E_FIVEPRIME_UTR,"5_prime_UTR_variant",5.0 },
		{ E_CODINGGAIN,"inframe_codon_gain",15.0 },
		{ E_INFRAME_INSERTION,"inframe_insertion",15.0 },
		{ E_CODINGLOSS,"inframe_codon_loss",15.0 },
		{ E_INFRAME_DELETION,"inframe_deletion",15.0 },
		{ E_SPLICE_DONOR, "splice_donor_variant",5.0 },
		{ E_SPLICE_SITE, "splice_region_variant",5.0 },
		{ E_STOP_LOST, "stop_lost",5.0 },
		{ E_NON_SYNONYMOUS_CODING,"non_synonymous_codon",10.0 },
		{ E_MISSENSE_VARIANT,"missense_variant",10.0 },
		{ E_CODON_CHANGE,"initiator_codon_change",10.0 },
{ E_NON_TERMINAL_CODON,"incomplete_terminal_codon_variant",10.0 },
{ E_MIRNA_VARIANT,"mature_miRNA_variant",10.0 },
{ E_FRAMESHIFT_CODING,"frameshift_variant", 20 },
{ E_START_LOST,"start_lost", 20 },
{ E_ESSENTIAL_SPLICE_SITE,"splice_acceptor_variant",10.0 },
{ E_TRANSCRIPT_CHANGE,"complex_change_in_transcript",20.0 },
		{ E_STOP_GAINED,"stop_gained",20.0 }
};

const char* consequenceString(consequenceType t,int use_ensembl)
{
	if (use_ensembl)
		for(int i=0;i<E_NCONSEQUENCETYPES;++i)
			if(t==e_consequence[i].t)
				return e_consequence[i].str;
	else
		for (int i=0;i<NCONSEQUENCETYPES;++i)
			if (t==consequence[i].t)
				return consequence[i].str;
}
