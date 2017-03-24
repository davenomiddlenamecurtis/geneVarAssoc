#if 0
Essential splice site	In the first 2 or the last 2 basepairs of an intron	splice_acceptor_variant	SO:0001574	splice-3
splice_donor_variant	SO:0001575	splice-5
Stop gained	In coding sequence, resulting in the gain of a stop codon	stop_gained	SO:0001587	nonsense
Stop lost	In coding sequence, resulting in the loss of a stop codon	stop_lost	SO:0001578	-
Complex in/del	Insertion or deletion that spans an exon/intron or coding sequence/UTR border	complex_change_in_transcript	SO:0001577	-
Frameshift coding	In coding sequence, resulting in a frameshift	frameshift_variant	SO:0001589	frameshift
Non-synonymous coding	In coding sequence and results in an amino acid change in the encoded peptide sequence	initiator_codon_change	SO:0001582	-
inframe_codon_loss	SO:0001652	-
inframe_codon_gain	SO:0001651	-
non_synonymous_codon	SO:0001583	missense
Splice site	1-3 bps into an exon or 3-8 bps into an intron	splice_region_variant	SO:0001630	-
Partial codon	Located within the final, incomplete codon of a transcript whose end coordinate is unknown	incomplete_terminal_codon_variant	SO:0001626	-
Synonymous coding	In coding sequence, not resulting in an amino acid change (silent mutation)	stop_retained_variant	SO:0001567	-
synonymous_codon	SO:0001588	cds-synon
Coding unknown	In coding sequence with indeterminate effect	coding_sequence_variant	SO:0001580	-
Within mature miRNA	Located within a microRNA	mature_miRNA_variant	SO:0001620	-
5 prime UTR	In 5 prime untranslated region	5_prime_UTR_variant	SO:0001623	untranslated_5
3 prime UTR	In 3 prime untranslated region	3_prime_UTR_variant	SO:0001624	untranslated_3
Intronic	In intron	intron_variant	SO:0001627	intron
NMD transcript	Located within a transcript predicted to undergo nonsense-mediated decay	NMD_transcript_variant	SO:0001621	-
Within non-coding gene	Located within a gene that does not code for a protein	nc_transcript_variant	SO:0001619	-
Upstream	Within 5 kb upstream of the 5 prime end of a transcript	2KB_upstream_variant	SO:0001636	near-gene-5
5KB_upstream_variant	SO:0001635	-
Downstream	Within 5 kb downstream of the 3 prime end of a transcript	500B_downstream_variant	SO:0001634	near-gene-3
5KB_downstream_variant	SO:0001633	-
Regulatory region	In regulatory region annotated by Ensembl	regulatory_region_variant	SO:0001566	-
Transcription factor binding motif	Falls in a transcription factor binding motif within an Ensembl regulatory region	TF_binding_site_variant	SO:0001782	-
Intergenic	More than 5 kb either upstream or downstream of a transcript	intergenic_variant	SO:0001628	-

Newer version, 29/12/16, from http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences

NULL_CONSEQUENCE	0	1	NULL_CONSEQUENCE{ E_NULL_CONSEQUENCE,"NULL_CONSEQUENCE",1.0 },
intergenic_variant	1	1	INTERGENIC_VARIANT{ E_INTERGENIC_VARIANT,"intergenic_variant",1.0 },
feature_truncation	2	3	FEATURE_TRUNCATION{ E_FEATURE_TRUNCATION,"feature_truncation",3.0 },
regulatory_region_variant	3	3	REGULATORY_REGION_VARIANT{ E_REGULATORY_REGION_VARIANT,"regulatory_region_variant",3.0 },
feature_elongation	4	3	FEATURE_ELONGATION{ E_FEATURE_ELONGATION,"feature_elongation",3.0 },
regulatory_region_amplification	5	3	REGULATORY_REGION_AMPLIFICATION{ E_REGULATORY_REGION_AMPLIFICATION,"regulatory_region_amplification",3.0 },
regulatory_region_ablation	6	3	REGULATORY_REGION_ABLATION{ E_REGULATORY_REGION_ABLATION,"regulatory_region_ablation",3.0 },
TF_binding_site_variant	7	3	TF_BINDING_SITE_VARIANT{ E_TF_BINDING_SITE_VARIANT,"TF_binding_site_variant",3.0 },
TFBS_amplification	8	3	TFBS_AMPLIFICATION{ E_TFBS_AMPLIFICATION,"TFBS_amplification",3.0 },
TFBS_ablation	9	3	TFBS_ABLATION{ E_TFBS_ABLATION,"TFBS_ablation",3.0 },
downstream_gene_variant	10	3	DOWNSTREAM_GENE_VARIANT{ E_DOWNSTREAM_GENE_VARIANT,"downstream_gene_variant",3.0 },
upstream_gene_variant	11	3	UPSTREAM_GENE_VARIANT{ E_UPSTREAM_GENE_VARIANT,"upstream_gene_variant",3.0 },
non_coding_transcript_variant	12	3	NON_CODING_TRANSCRIPT_VARIANT{ E_NON_CODING_TRANSCRIPT_VARIANT,"non_coding_transcript_variant",3.0 },
NMD_transcript_variant	13	3	NMD_TRANSCRIPT_VARIANT{ E_NMD_TRANSCRIPT_VARIANT,"NMD_transcript_variant",3.0 },
intron_variant	14	3	INTRON_VARIANT{ E_INTRON_VARIANT,"intron_variant",3.0 },
non_coding_transcript_exon_variant	15	3	NON_CODING_TRANSCRIPT_EXON_VARIANT{ E_NON_CODING_TRANSCRIPT_EXON_VARIANT,"non_coding_transcript_exon_variant",3.0 },
3_prime_UTR_variant	16	5	3_PRIME_UTR_VARIANT{ E_3_PRIME_UTR_VARIANT,"3_prime_UTR_variant",5.0 },
5_prime_UTR_variant	17	5	5_PRIME_UTR_VARIANT{ E_5_PRIME_UTR_VARIANT,"5_prime_UTR_variant",5.0 },
mature_miRNA_variant	18	5	MATURE_MIRNA_VARIANT{ E_MATURE_MIRNA_VARIANT,"mature_miRNA_variant",5.0 },
coding_sequence_variant	19	5	CODING_SEQUENCE_VARIANT{ E_CODING_SEQUENCE_VARIANT,"coding_sequence_variant",5.0 },
synonymous_variant	20	5	SYNONYMOUS_VARIANT{ E_SYNONYMOUS_VARIANT,"synonymous_variant",5.0 },
stop_retained_variant	21	5	STOP_RETAINED_VARIANT{ E_STOP_RETAINED_VARIANT,"stop_retained_variant",5.0 },
incomplete_terminal_codon_variant	22	5	INCOMPLETE_TERMINAL_CODON_VARIANT{ E_INCOMPLETE_TERMINAL_CODON_VARIANT,"incomplete_terminal_codon_variant",5.0 },
splice_region_variant	23	5	SPLICE_REGION_VARIANT{ E_SPLICE_REGION_VARIANT,"splice_region_variant",5.0 },
protein_altering_variant	24	10	PROTEIN_ALTERING_VARIANT{ E_PROTEIN_ALTERING_VARIANT,"protein_altering_variant",10.0 },
missense_variant	25	10	MISSENSE_VARIANT{ E_MISSENSE_VARIANT,"missense_variant",10.0 },
inframe_deletion	26	15	INFRAME_DELETION{ E_INFRAME_DELETION,"inframe_deletion",15.0 },
inframe_insertion	27	15	INFRAME_INSERTION{ E_INFRAME_INSERTION,"inframe_insertion",15.0 },
transcript_amplification	28	15	TRANSCRIPT_AMPLIFICATION{ E_TRANSCRIPT_AMPLIFICATION,"transcript_amplification",15.0 },
start_lost	29	15	START_LOST{ E_START_LOST,"start_lost",15.0 },
stop_lost	30	15	STOP_LOST{ E_STOP_LOST,"stop_lost",15.0 },
frameshift_variant	31	20	FRAMESHIFT_VARIANT{ E_FRAMESHIFT_VARIANT,"frameshift_variant",20.0 },
stop_gained	32	20	STOP_GAINED{ E_STOP_GAINED,"stop_gained",20.0 },
splice_donor_variant	33	20	SPLICE_DONOR_VARIANT{ E_SPLICE_DONOR_VARIANT,"splice_donor_variant",20.0 },
splice_acceptor_variant	34	20	SPLICE_ACCEPTOR_VARIANT{ E_SPLICE_ACCEPTOR_VARIANT,"splice_acceptor_variant",20.0 },
transcript_ablation	35	20	TRANSCRIPT_ABLATION{ E_TRANSCRIPT_ABLATION,"transcript_ablation",20.0 },

#endif