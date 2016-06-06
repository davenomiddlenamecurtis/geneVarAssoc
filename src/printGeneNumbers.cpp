#include <stdlib.h>
#include "geneVarUtils.hpp"

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	FILE *fp,*fo[2];
	gvaParams gp;
	analysisSpecs spec;
	float mlp;
	int problemExtractingGene;
	long varNum;
	allelePair **a;
	int i,k,l; k=0;
	varNum=1L;
	fp=fopen(argv[1],"r");
	if (argc>3)
	{
		gp.firstGeneNum=atoi(argv[3]);
		gp.lastGeneNum=atoi(argv[4]);
	}
	else
	{
		gp.firstGeneNum=0;
		gp.lastGeneNum=30000;
	}
	gp.input(fp,spec);
	fclose(fp);	

	masterLocusFile vf(gp.nCc[1]);
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setUpstream(gp.downstream);
	r.setBaitMargin(gp.margin);
	r.goToStart();
	while (r.getNextGene())
	{
		++k;
		printf("%6d %s\n",k,r.getGene());
		if (k>gp.lastGeneNum)
			break;
		if (k<gp.firstGeneNum)
			continue;
	}
	return 0;
}
