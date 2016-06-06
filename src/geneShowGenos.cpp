#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	int i;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;
	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
	strcpy(geneName,argv[2]);
	
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setBaitMargin(gp.margin);

	if (!r.findGene(geneName) || !r.getNextGene())
		return 1;
	sprintf(fn,"gva.%s.db",geneName);
	sprintf(fn2,"gva.%s.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
		for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gva.%s.cont.%d.vcf",geneName,i+1);
			gcont.extractGene(r,fn);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,0);
		}
		for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gva.%s.case.%d.vcf",geneName,i+1);
			gcase.extractGene(r,fn);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,1);
		}
	if (spec.useEnsembl)
			vf.getEnsemblConsequences(spec);
	else
			vf.getQuickConsequences(r,spec);
	sprintf(fn,"gsg.%s.txt",geneName);
	vf.writeGenos(fn,gp.useFreqs,spec);
	return 0;
}