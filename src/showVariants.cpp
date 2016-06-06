#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "getGene.hpp"
#include "masterLocusFile.hpp"

#define FREQFILE "c:\\reference\\ALL.wgs.phase1_release_v2.20101123.snps_indels_sv.sites.vcf.gz"

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
	strcpy(geneName,argv[1]);
	
	r.setListFile(RSGFILE);
	r.setReferencePath(REFERENCEPATH);
	r.setUpstream(1000);
	gcont.setVariantFileName(FREQFILE);

	unlink("sv.vcf");
	unlink("sv.db");
	unlink("sv.vdx");

    sprintf(fn,"%s.vcf",argv[1]);
	if (r.findGene(geneName) && r.getNextGene())
		gcont.extractGene(r,fn);
;
	vf.openFiles("sv.db","sv.vdx");
	vf.addLocusFile(fn,VCFFILE);
	vf.readLocusFileEntries(fn);
	vf.getEnsemblConsequences();
	vf.getQuickConsequences(r);
		vf.printFeatures(stdout,0);

	return 0;
}