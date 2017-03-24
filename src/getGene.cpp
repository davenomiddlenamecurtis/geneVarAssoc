#include "getGene.hpp"


// getgene - download 1000 genomes entry for a gene
#if 0
example:
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/ALL.chr17.phase1_integrated_calls.20101123.snps_indels_svs.genotypes.vcf.gz 17:64873348-64881099 > CACNG5genotypes.tabix
never finished this code and I do not know whether tabix works over ftp
#endif


int main(int argc,char *argv[])
{
	ggParams p;
	refseqGeneInfo g;
	p.getParams(argc,argv);
	if (p.findGene(argv[1],g) && p.getNextGene(g) && g.checkExonLengths())
		p.downloadGene(g);
	return 0;
}
