#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"
#include "vcfLocusFile.hpp"

void removeSlash(char *s)
{
	char *ptr;
	ptr=s;
	do
	{
		if (*ptr == '/')
			*ptr = '-';
	} while (*ptr++);
}

int writeShapeItInput(masterLocusFile &vf,gvaParams &gp,analysisSpecs &spec,char *root,int nSub,long *fP,long *lP)
{
	char fn[100];
	char *genoStr[] = { "0 0 0", "W T F", "1 0 0", "0 1 0", "0 0 1" };
	int s,g;
	long firstPos,lastPos,pos;
	float mapPos,lastMapPos,called;
	FILE *fp,*fm;
	allelePair *all;
	strEntry *subName;
	subName=(strEntry *)calloc(nSub,sizeof(strEntry));
	vf.outputSubNames(subName,spec);
	sprintf(fn,"%s.sample",root);
	fp=fopen(fn,"w");
	fprintf(fp,"ID_1 ID_2 missing\n0 0 0\n");
	for (s=0;s<nSub;++s)
		fprintf(fp,"%s %s 0\n",subName[s],subName[s]);
	fclose(fp);
	free(subName);
	all=(allelePair *)calloc(nSub,sizeof(allelePair));

	sprintf(fn,"%s.gen",root);
	fp=fopen(fn,"w");
	sprintf(fn,"%s.map",root);
	fm=fopen(fn,"w");
	vf.gotoFirstInRange(spec);
#define cMperMb 0.6
	lastPos=-1;
	lastMapPos=0;
	fprintf(fm,"position	COMBINED_rate.cM.Mb.	Genetic_Map.cM.\n");
	do {
		if (vf.currentNAlls()!=2)
			continue;
		if (spec.onlyUseSNPs && !vf.currentIsSNP())
			continue;
		vf.outputCurrentAlleles(all,spec);
		called=0;
		for (s=0;s<nSub;++s)
			if (all[s][0]!=0)
				++called;
		if (called/nSub<spec.proportionCalledToPass)
			continue; // omit variants which do not have many calls
		fprintf(fp,"%03d %10.10s %ld %s %s  ",
			vf.currentChr(),(*(vf.currentID())=='\0')?".":vf.currentID(),pos=vf.currentPos(),vf.currentAll(0),vf.currentAll(1));
		for (s = 0; s < nSub; ++s)
		{
			g = all[s][0] + all[s][1];
			fprintf(fp, "%s  ", genoStr[g]);
		}
		fprintf(fp,"\n");
		if (lastPos == -1)
			firstPos=lastPos=pos;
		mapPos=lastMapPos+(pos-lastPos)*cMperMb/1000000.0;
		lastMapPos=mapPos;
		lastPos=pos;
		fprintf(fm, "%010ld %f %f\n",pos,cMperMb,mapPos);
	} while (vf.gotoNextInRange(spec));
	fclose(fp);
	fclose(fm);
	*fP=firstPos;
	*lP=lastPos;
	return 1;
}

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	int i,extractedOK;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;
	long firstPos,lastPos;

	dcerror.warn();
	hereOK();
	printf("Running program\n");
	fprintf(stderr,"Output on stderr\n");

	fp=fopen(argv[1],"r");
	if (fp == NULL)
		{ dcerror(1, "Could not open file %s", argv[1]); exit(1); }
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
	removeSlash(geneName);
	sprintf(fn,"gtr.%s.db",geneName);
	sprintf(fn2,"gtr.%s.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	extractedOK=1;
	int ff=0;
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gtr.%s.cont.%d.vcf",geneName,i+1);
			if (!gcont.extractGene(r,fn,0,spec.addChrInVCF[ff++]))
				extractedOK=0;
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,0);
		}
	for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gtr.%s.case.%d.vcf",geneName,i+1);
			if (!gcase.extractGene(r,fn,0,spec.addChrInVCF[ff++]))
				extractedOK=0;
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,1);
		}
	if (spec.consequenceThreshold!=0 || spec.useConsequenceWeights!=0)
	{
		if (spec.useEnsembl)
			vf.getEnsemblConsequences(spec);
		else
			vf.getQuickConsequences(r,spec);
	}
	sprintf(fn,"gtr.%s",geneName);
	writeShapeItInput(vf,gp,spec,fn,vf.getTotalSubs(),&firstPos,&lastPos);
#ifdef MSDOS
	sprintf(fn,"shapeit.%s.bat",geneName);
#else
	sprintf(fn,"shapeit.%s.sh",geneName);
#endif
	fp=fopen(fn,"w");
hereOK();
/*
	fprintf(fp,"shapeit --effective-size 100000 --input-gen gtr.%s.gen --input-map gtr.%s.map -O gtr.%s.phased \n",
		geneName,geneName,geneName);
*/
/*
was
fprintf(fp,"shapeit --effective-size 100000 --input-gen gtr.%s --input-map gtr.%s.map -O gtr.%s.phased \n",
		geneName,geneName,geneName,geneName);
*/
	/*
	https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#ex10
	*/
	fprintf(fp,"impute2 -phase -m gtr.%s.map -g gtr.%s.gen -int %ld %ld -Ne 20000 -o gtr.%s.phased",
		geneName,geneName,firstPos,lastPos,geneName);
	fclose(fp);
#ifdef MSDOS
	sprintf(line,"%s",fn);
#else
	sprintf(line,"sh %s",fn);
#endif
	system(line);
	masterLocusFile hapsf(1); // use for haplotypes file
hereOK();
	sprintf(fn,"gtr.%s.haps.db",geneName);
	sprintf(fn2,"gtr.%s.haps.vdx",geneName);
	unlink(fn);
	unlink(fn2);
hereOK();
	hapsf.openFiles(fn,fn2);
hereOK();
	sprintf(fn,"gtr.%s.phased_haps",geneName);
hereOK();
	hapsf.addLocusFile(fn,SHAPEITHAPSFILE);
hereOK();
	hapsf.readLocusFileEntries(fn,spec,0);
hereOK();
	hapsf.getQuickConsequences(r,spec);
hereOK();
	sprintf(fn,"gtr.%s",geneName);
		
	spec.useHaplotypes=1;
hereOK();
	hapsf.writeOldScoreAssocFiles(vf,fn,gp.wf,gp.wFunc,gp.useFreqs,gp.nSubs,1,gp.writeComments,spec);
hereOK();

	sprintf(line,"scoreassoc %s.par %s.dat %s.sao",fn,fn,fn);
hereOK();
	if (gp.writeScoreFile==1)
		sprintf(strchr(line,'\0')," %s.sco",fn);
	system(line);

	return 0;
}