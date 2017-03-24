#include <stdlib.h>
#include <string.h>
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

int getHapsFromPhase(masterLocusFile &vf,gvaParams &gp,analysisSpecs &spec,char *root,int nSub)
{
	char line[1000],fn[100],*ptr;
	int s,g,l,nLoci,aa,a,chr,nLociToUse,*useLocus,ll;
	float called;
	FILE *fp,*fm;
	allelePair **all;
	strEntry *subName,*pos,*description;
	nLoci=0;
	vf.gotoFirstInRange(spec);
	chr=vf.currentChr();
	do {
		if (vf.currentNAlls()!=2)
			continue;
		if (spec.onlyUseSNPs && !vf.currentIsSNP())
			continue;
		++nLoci;
	} while (vf.gotoNextInRange(spec));

	subName=(strEntry *)calloc(nSub,sizeof(strEntry));
	vf.outputSubNames(subName,spec);
	all=(allelePair **)calloc(nLoci,sizeof(allelePair*));
	for (l=0;l<nLoci;++l)
		all[l]=(allelePair *)calloc(nSub,sizeof(allelePair));
	pos=(strEntry *)calloc(nLoci,sizeof(strEntry));
	description=(strEntry *)calloc(nLoci,sizeof(strEntry));
	useLocus=(int *)calloc(nLoci,sizeof(int));
	nLociToUse=0; 
	l=0;
	vf.gotoFirstInRange(spec);
	do {
		if (vf.currentNAlls()!=2)
			continue;
		if (spec.onlyUseSNPs && !vf.currentIsSNP())
			continue;
		sprintf(pos[l],"%ld",vf.currentPos());
		sprintf(description[l],"%ld %s %s ",vf.currentPos(),vf.currentAll(0),vf.currentAll(1));
		vf.outputCurrentAlleles(all[l],spec);
		called=0;
		for (s=0;s<nSub;++s)
			if (all[l][s][0]!=0)
				++called;
		if (called/nSub>=spec.proportionCalledToPass) // this needs called to be a float not an int
			// I think PHASE will subsequently impute genotypes if these are unknown
		{
			++nLociToUse;
			useLocus[l]=1;
		}
		++l;
	} while (vf.gotoNextInRange(spec));

	sprintf(fn,"%s.inp",root);
	fp=fopen(fn,"w");
	fprintf(fp,"%d\n%d\nP ",nSub,nLociToUse);
	for (l=0;l<nLoci;++l)
		if (useLocus[l])
			fprintf(fp,"%s ",pos[l]);
	fprintf(fp,"\n");
	for (l=0;l<nLociToUse;++l)
		fprintf(fp,"S");
	fprintf(fp,"\n");
	for (s = 0; s < nSub; ++s)
	{
		fprintf(fp,"%s\n",subName[s]);
		for (a = 0; a < 2; ++a)
		{
			for (l = 0; l < nLoci; ++l)
				if (useLocus[l])
					fprintf(fp,"%c ",(aa=all[l][s][a])==0?'?':'0'+aa);
			fprintf(fp,"\n");
		}
	}
	fclose(fp);
	if (nLociToUse==0)
		return 0;
	#ifdef MSDOS
	sprintf(fn,"%s.bat",root);
#else
	sprintf(fn,"%s.sh",root);
#endif
	fp=fopen(fn,"w");
hereOK();
#if 0
	fprintf(fp,"phase %s.inp %s.out 50 1 50",
		root,root);
#else
	fprintf(fp,"phase %s.inp %s.out",
		root,root);
#endif
	fclose(fp);
#ifdef MSDOS
	sprintf(line,"%s",fn);
#else
	sprintf(line,"sh %s",fn);
#endif
hereOK();
	system(line);
hereOK();

	sprintf(fn,"%s.out_pairs",root);
hereOK();
	fp=fopen(fn,"r");
	s=0;
	while (fgets(line, 999, fp))
	{
		while (strncmp(line,"IND:",4))
			if (!fgets(line,999,fp))
				goto done;
		if (!fgets(line,999,fp))
			goto done;
		ptr=line;
		for (a = 0; a < 2; ++a)
		{
			for (l = 0; l < nLociToUse; ++l)
				all[l][s][a]=*ptr++-'0';
			if (a==0)
				ptr+=3;
			else break;
		}
		++s;
	}
done:
	fclose(fp);
hereOK();
	sprintf(fn,"%s.phased_haps",root);
	fp=fopen(fn,"w");
	for (ll=l = 0; l < nLociToUse; ++l,++ll)
	{
		while (useLocus[ll]==0)
			++ll;
		fprintf(fp,"%03d . %s",chr,description[ll]);
		for (s=0;s<nSub;++s)
			for (a=0;a<2;++a)
				fprintf(fp,"%d ",all[l][s][a]-1);
		fprintf(fp,"\n");
	}
	fclose(fp);

	for (l=0;l<nLoci;++l)
		free(all[l]);
	free(all);
	free(subName);
	free(pos);
	free(description);
	free(useLocus);
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
	extractedOK=1;
	vf.openFiles(fn,fn2);
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
	if (!extractedOK)
		{
		sprintf(fn,"gtr.%s.sao",geneName);
		fp=fopen(fn,"w");
		fprintf(fp,"Could not extract variants for this gene\n");
		fclose(fp);
		exit(1);
		}
	if (vf.currentChr()>22)
		{
		sprintf(fn,"gtr.%s.sao",geneName);
		fp=fopen(fn,"w");
		fprintf(fp,"Cannot perform phased analysis on sex chromosomes\n");
		fclose(fp);
		exit(1);
		}
	if (spec.consequenceThreshold!=0 || spec.useConsequenceWeights!=0)
	{
		if (spec.useEnsembl)
			vf.getEnsemblConsequences(spec);
		else
			vf.getQuickConsequences(r,spec);
	}
	sprintf(fn,"gtr.%s",geneName);
	if (getHapsFromPhase(vf, gp, spec, fn, vf.getTotalSubs()) == 0)
	{
		// failed to run for some known reason - no informative variants
		sprintf(fn,"gtr.%s.sao",geneName);
		fp=fopen(fn,"w");
		fprintf(fp,"Did not run phase - no qualifying variants\n");
		fclose(fp);
		exit(1);
	}
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