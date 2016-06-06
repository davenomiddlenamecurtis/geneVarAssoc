#include <stdlib.h>
#include "getGene.hpp"
#include "masterLocusFile.hpp"

#define FREQFILE "EUR.2of4freqs.vcf.gz"
#define RSGFILE "c:\\reference\\refseqgeneshgALL040212.txt"
#define REFERENCEPATH "c:\\reference"

int main(int argc,char *argv[])
{
	geneExtractor gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	masterLocusFile vf(1);
	int ns;
	float wf;
	FILE *fp;
	consequenceType e;
	
	r.setListFile(RSGFILE);
	r.setReferencePath(REFERENCEPATH);
	r.setUpstream(1000);
	r.findGene("PIGV");
	r.getNextGene();
	e=r.getEffect(27121680,"C","T");
	printf("Effect is %d\n",e);
	return 0;
}