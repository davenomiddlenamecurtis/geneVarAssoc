/* scorevcf.c */
#if 0
Read in one or more vcf files containing variants for cases and controls.
Produce a sorted list of all variants found in any file.
Sort by position then ?string describing variant?
Variants noted in one file may be absent from other file(s).
SNPs may be biallelic in one file but triallelic in another.
So sort SNP just by position and will have to keep track of how many alleles found.
Also, wt base may be different in each file and alternative base[s] may be different in each file.

May be two or more variants at same position.

Master list will be ordered by position then variant, SNP first. 
Merge all SNPs into one record.
For each record keep note of each file it was found in. Need allele list for that file. May as well have line number in that file too.
Other variants will be on subsequent lines.

Have index which points to first variant at that position - other variants will be in subsequent entries.

20012 SNP A,C,G File1 A,G File2 C,G File3 A,C,G
21133 SNP A,C File1 A,C File3 A,C
21133 VAR varstring File1 
21133 VAR varother File3
24041 SNP G,T File1 G,T File2 T,G

Then read in VCF files again assigning genotypes to subjects.
If variant absent from a file then assume homozygous wt rather than "unknown".

Recognise SNP as: 
A      G,T
A      G
A      .
- single char followed by one or more comma separated chars, or . if monomorphic in this sample.
Problem with last one is that cannot tell is definitely SNP. Could be reference for something else.
I do not know what N means as a base name.

Other variants e.g.:
GTC	G,GTCTC
C	<INV>
N	.[13:123457[
Probably just strip out whitespace and store resulting string, e.g. N_.[13:123457[
Note that this will fail if some microsatellite alleles are present in some files but not others:
GTC	G,GTCTC -> GTC_G,GTCTC
GTC	G -> GTC_G
GTC	GTCTC -> GTC_GTCTC and in fact would probably be C	CTC -> C_CTC at a position two further down
But this is just too bad for now.

Dimension of array needs to be number of different variants. 
Dimension of index needs to be number of bases in chromosome.
If want to handle whole genome probably easiest to have 24 indices, one for each chromosome.
Could have separate variant arrays or one long one. Probably only 100K for exomes but v. large for genomes.



Then analyse!

OK, try again. Have a master VCF record which keeps the allele assignments for each vcf file.
#endif