#include <stdlib.h>
#include <string.h>

#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"
#include "riskLocus.hpp"
#include <ctype.h>

extern int outputGenotypes(FILE *fo, masterLocusFile *vf, analysisSpecs spec,int * cc,allelePair *a);

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100],groupName[100];
	int i,first;
	FILE *fp,*fg,*fo;
	gvaParams gp;
	analysisSpecs spec;
	int *cc;
	allelePair *a;
	strEntry *subName;
	fp=fopen(argv[1],"r");
	if (fp == NULL)
		{ dcerror(1, "Could not open file %s", argv[1]); exit(1); }
	gp.input(fp,spec);
	fclose(fp);
	strcpy(groupName,argv[2]);
	
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setBaitMargin(gp.margin);
	fg=fopen(groupName,"r");
	fo=fopen(argv[3],"w");
	dcerror.warn();
	first=1;
	while (fgets(line,999,fg) && sscanf(line,"%s",geneName)==1)
	{
	if (!r.findGene(geneName) || !r.getNextGene())
		continue;
//	if (!r.findGene(geneName) || !r.getNextGene())
//		return 1;
	sprintf(fn,"gGG.temp.db");
	sprintf(fn2,"gGG.temp.vdx");
	unlink(fn);
	unlink(fn2);
	masterLocusFile *vf=new masterLocusFile(gp.nCc[0]+gp.nCc[1]);
	vf->openFiles(fn,fn2);
	int ff=0;
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gGG.cont.%d.vcf",i+1);
			gcont.extractGene(r,fn,0,spec.addChrInVCF[ff++]);
		}
	for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gGG.case.%d.vcf",i+1);
			gcase.extractGene(r,fn,0,spec.addChrInVCF[ff++]);
		}
	for (i=0;i<gp.nCc[0];++i)
	{
		sprintf(fn,"gGG.cont.%d.vcf",i+1);
		vf->addLocusFile(fn,VCFFILE);
		vf->readLocusFileEntries(fn,spec,0);
	}
	for (i=0;i<gp.nCc[1];++i)
	{
		sprintf(fn,"gGG.case.%d.vcf",i+1);
		vf->addLocusFile(fn,VCFFILE);
		vf->readLocusFileEntries(fn,spec,1);
	}
	vf->getQuickConsequences(r,spec,1);
	if (first == 1)
	{
		first=0;
		subName=(strEntry *)calloc(vf->getTotalSubs(),sizeof(strEntry));
		vf->outputSubNames(subName,spec);
		cc=(int *)calloc(vf->getTotalSubs(),sizeof(int));
		a=new allelePair[vf->getTotalSubs()];
		vf->outputAffectionStatus(cc,spec);

		for (i = 0; i < vf->getTotalSubs(); ++i)
			fprintf(fo,"%s ",subName[i]);
		fprintf(fo,"\n");
		for (i = 0; i < vf->getTotalSubs(); ++i)
			fprintf(fo,"%d ",cc[i]);
		fprintf(fo,"\n");
	}
	outputGenotypes(fo,vf,spec,cc,a);
	delete vf;
	}
	free(subName);
	free(cc);
	delete[] a;
	fclose(fg);
	fclose(fo);
	return 0;
}