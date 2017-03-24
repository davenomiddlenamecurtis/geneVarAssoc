#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

int masterLocusFile::writeGenoCounts(FILE *fo[2],char *geneName,long *varNum,analysisSpecs &spec,allelePair **a)
{
	char fn[100],buff[1000],buff2[20],comment[1000],*ptr,alleles[MAXSTR+1];;
	int totalSub,lc,s,l,ss,i,c,gc[2][3],nFiles,f;
	long recPos;
	const char *testKey;
	nFiles=spec.phenotypes?2:1;
	openLocusFiles();
	for (i=0,totalSub=0;i<nLocusFiles;++i)
		totalSub+=nSubs[i];
	lc=outputAlleles(a,spec);
	l=0;
	recPos=findFirstInRange(spec);
	if (recPos!=0L)
	{
		while (1)
		{
		testKey=index.current_key();
		if ((c=atoi(testKey))==0 || c>spec.ec)
			break;
		if (c==spec.ec && atol(testKey+3)>spec.ep)
			break;
		load(tempRecord,recPos);
		for (f=0;f<nFiles;++f)
		{
			fprintf(fo[f]," %12ld %20s ",*varNum,geneName);
			fprintf(fo[f],"%02d %010ld ",tempRecord.chr,tempRecord.pos);
		}
		for (f=0;f<nFiles;++f)
			gc[f][0]=gc[f][1]=gc[f][2]=0;
		for (s=0,i=0;i<nLocusFiles;++i)
			for (ss=0;ss<nSubs[i];++s,++ss)
				{
					if (a[l][s][0]==0)
						continue; // unknown
					else
					{
						c=(a[l][s][0]>1)+(a[l][s][1]>1);
						if (spec.phenotypes==0)
							++gc[0][c];
						else
							++gc[spec.phenotypes[s]][c];
					}
				}
		for (f=0;f<nFiles;++f)
		{
			fprintf(fo[f],"%5d %5d %5d   ",gc[f][0],gc[f][1],gc[f][2]);
			if (tempRecord.ensemblConsequence[0]!='\0')
				fprintf(fo[f],"%s",tempRecord.ensemblConsequence);
			else if (tempRecord.quickConsequence[0]!='\0')
				fprintf(fo[f],"%s",tempRecord.quickConsequence);
			fprintf(fo[f],"\n");
		}
		++l;
		++ *varNum;
		recPos=index.get_next();
		if (recPos==0L)
			break;
		}
	}
	return 1;
}

#define MAXLOCIINGENE 50000
#define MAXSUB 10000

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

#if 0
	only use cases
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
#else
	masterLocusFile vf(gp.nCc[1]);
#endif	
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setUpstream(gp.downstream);
	r.setBaitMargin(gp.margin);
	r.goToStart();
	fo[0]=fopen(argv[2],"w");
	if (spec.phenotypes)
		fo[1]=fopen(argv[5],"w");
	a=(allelePair **)calloc(MAXLOCIINGENE,sizeof(allelePair*));
	for (l=0;l<MAXLOCIINGENE;++l)
		a[l]=(allelePair *)calloc(MAXSUB,sizeof(allelePair));
	while (r.getNextGene())
	{
		printf("%s\n",r.getGene());
		++k;
		if (k>gp.lastGeneNum)
			break;
		if (k<gp.firstGeneNum)
			continue;
		unlink("gva.all.db");
		unlink("gva.all.vdx");
		vf.openFiles("gva.all.db","gva.all.vdx");
		problemExtractingGene=0;
#if 0
		only count cases so we can use same parameter files
		for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gva.all.cont.%d.vcf",i+1);
			unlink(fn);
			if (gcont.extractGene(r,fn)==0)
			{
				problemExtractingGene=1;
				break;
			}
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,0);
		}
#endif
		int ff=0;
		for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gva.all.case.%d.vcf",i+1);
			unlink(fn);
			if (gcase.extractGene(r,fn,0,spec.addChrInVCF[ff++])==0)
			{
				problemExtractingGene=1;
				break;
			}
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,1);
		}
		if (problemExtractingGene==1)
		{
			vf.closeLocusFiles();
			vf.closeFiles();
			continue;
		}

		if (spec.consequenceThreshold!=0 || spec.useConsequenceWeights!=0)
		{
			if (spec.useEnsembl)
				vf.getEnsemblConsequences(spec);
			else
				vf.getQuickConsequences(r,spec);
		}
		vf.writeGenoCounts(fo,r.getGene(),&varNum,spec,a);
		vf.closeLocusFiles();
		vf.closeFiles();
	}
	fclose(fo[0]);
	if (spec.phenotypes)
		fclose(fo[1]);
	for (l=0;l<MAXLOCIINGENE;++l)
		free(a[l]);
	free(a);

	return 0;
}
