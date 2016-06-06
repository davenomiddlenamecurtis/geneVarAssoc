#include "getGene.hpp"

int main(int argc,char*argv[])
{
	refseqGeneInfo r;
	char fn[100];
	r.setListFile("c:\\sequence\\refseqgeneshgALL040212.txt");
	r.setUpstream(1000);
	sprintf(fn,"%s.downloaded.vcf",argv[1]);
	if (r.findGene(argv[1]) && r.getNextGene())
	{
		r.checkExonLengths();
	}
	return 0;
}