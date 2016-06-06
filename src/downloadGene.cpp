#include "getGene.hpp"
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "masterLocusFile.hpp"

int main(int argc,char*argv[])
{
	geneExtractor g;
	refseqGeneInfo r;
	masterLocusFile vf(1);
	char fn[100];
	r.setListFile("c:\\sequence\\refseqgeneshgALL040212.txt");
	r.setUpstream(1000);
	sprintf(fn,"%s.downloaded.vcf",argv[1]);
	if (r.findGene(argv[1]) && r.getNextGene())
		g.downloadGene(r,fn);
	unlink("dg.db");
	unlink("dg.vdx");
	vf.openFiles("dg.db","dg,vdx");
	vf.addLocusFile(fn,VCFFILE);
	vf.readLocusFileEntries(fn);
	vf.openLocusFiles(); // not sure if I need this
	vf.getQuickFeatures(r);
	vf.getFeatures();
	vf.print(stdout);

	return 0;
}