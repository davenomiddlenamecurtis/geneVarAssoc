#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100],groupName[100];
	int i,first;
	FILE *fp,*fg;
	gvaParams gp;
	analysisSpecs spec;
	fp=fopen(argv[1],"r");
	if (fp == NULL)
		{ dcerror(1, "Could not open file %s", argv[1]); exit(1); }
	gp.input(fp,spec);
	fclose(fp);
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
	strcpy(groupName,argv[2]);
	
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setBaitMargin(gp.margin);
	sprintf(fn,"gva.%s.db",groupName);
	sprintf(fn2,"gva.%s.vdx",groupName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	first=1;
	fg=fopen(groupName,"r");
	while (fgets(line,999,fg) && sscanf(line,"%s",geneName)==1)
	{
	if (!r.findGene(geneName) || !r.getNextGene())
		return 1;
	int ff=0;
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gva.%s.cont.%d.vcf",groupName,i+1);
			gcont.extractGene(r,fn,first==00,spec.addChrInVCF[ff++]);
		}
	for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gva.%s.case.%d.vcf",groupName,i+1);
			gcase.extractGene(r,fn,first==00,spec.addChrInVCF[ff++]);
		}
	first=0;
	}
	fclose(fg);
	for (i=0;i<gp.nCc[0];++i)
	{
		sprintf(fn,"gva.%s.cont.%d.vcf",groupName,i+1);
		vf.addLocusFile(fn,VCFFILE);
		vf.readLocusFileEntries(fn,spec,0);
	}
	for (i=0;i<gp.nCc[1];++i)
	{
		sprintf(fn,"gva.%s.case.%d.vcf",groupName,i+1);
		vf.addLocusFile(fn,VCFFILE);
		vf.readLocusFileEntries(fn,spec,1);
	}
	if (spec.consequenceThreshold!=0 || spec.useConsequenceWeights!=0)
	{
		if (spec.useEnsembl)
			vf.getEnsemblConsequences(spec);
		else
		{
			fg=fopen(groupName,"r");
			while (fgets(line,999,fg) && sscanf(line,"%s",geneName)==1)
			{
				if (!r.findGene(geneName) || !r.getNextGene())
					return 1;
				vf.getQuickConsequences(r,spec,1);
			}
		}
	}
	if (spec.consequenceThreshold!=0)
		sprintf(fn,"gva.%s.ct%02d",groupName,spec.consequenceThreshold);
	else if (spec.useConsequenceWeights!=0)
		sprintf(fn,"gva.%s.ucw",groupName);
	else
		sprintf(fn,"gva.%s",groupName);

	vf.writeOldScoreAssocFiles(fn,gp.wf,gp.wFunc,gp.useFreqs,gp.nSubs,1,gp.writeComments,spec);

	sprintf(line,"scoreassoc %s.par %s.dat %s.sao",fn,fn,fn);
	if (gp.writeScoreFile==1)
		sprintf(strchr(line,'\0')," %s.sco",fn);
	system(line);

	return 0;
}