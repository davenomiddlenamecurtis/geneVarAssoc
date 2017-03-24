#include <stdlib.h>
#include <ctype.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

// imputExomes.cpp

// use sequence data in VCF file to impute into SNP data in plink data file

enum commandFileType { BATCH=0, SCRIPT };
#ifdef MSDOS
commandFileType runScript=BATCH;
#else
commandFileType runScript=SCRIPT;
#endif

#ifdef MSDOS
#define RUNSHAPEIT "~/shapeit/shapeit"
#define RUNIMPUTE2 "\\impute2\\impute2"
#define RUNPLINK "\\plink107\\plink"
#endif

#define LONGLINELENGTH 200000
char longLine[LONGLINELENGTH+1];

#define MAXVARS 5000
class iePars {
public:
	char *runShapeIt,*runImpute2,*runPlink;
	int extraBases;  // on either side of gene
	float minSeqFreq,maxSeqFreq;
	int nFromSeq,fromSeq[MAXVARS],nVars;
	long pos[MAXVARS];
	iePars();
	int input(FILE *fp);
};

int iePars::input(FILE *fp)
{
	char line[1001];
	fgets(line,1000,fp);
	sscanf(line,"%d",&extraBases);
	fgets(line,1000,fp);
	sscanf(line,"%f %f",&minSeqFreq,&maxSeqFreq);
	return 1;
}

iePars::iePars()
{
	runShapeIt=RUNSHAPEIT;
	runImpute2=RUNIMPUTE2;
	runPlink=RUNPLINK;
	minSeqFreq=0;
	maxSeqFreq=1;
	extraBases=100000;
	nFromSeq=nVars=0;
}

#define BADPOS 1000000000L
int makeFakeMap(masterLocusFile &vSeq, gvaParams &gp, analysisSpecs &spec, refseqGeneInfo &r, char *root, iePars &pars)
{
	char fn[100],line[1001],*linePtr[2];
	FILE *inMap[2],*impMap;
	long pos[2],lastPos;
	int f,lineLength[2],lowf;
	float mapPos,lastMapPos;
	sprintf(fn,"ie.%s.SNPs.map",root);
	inMap[0]=fopen(fn,"r");
	sprintf(fn,"ie.%s.phased.haps",root);
	inMap[1]=fopen(fn,"r");
	sprintf(fn,"ie.%s.forImpute.map",root);
	impMap=fopen(fn,"w");

	pars.nVars=pars.nFromSeq=0;
	linePtr[0]=line;
	lineLength[0]=1000;
	linePtr[1]=longLine;
	lineLength[1]=LONGLINELENGTH;
	lastPos=BADPOS;
	for (f = 0; f < 2; ++f)
	{
		fgets(linePtr[f],lineLength[f], inMap[f]);
		pos[f]=BADPOS;
		sscanf(linePtr[f],f?("%*s %*s %ld"):("%*s %*s %*s %ld"),&pos[f]);
		if (pos[f]<lastPos)
			lastPos=pos[f]; // choose smallest
	}
	fprintf(impMap,"position	COMBINED_rate.cM.Mb.	Genetic_Map.cM.\n");
	lastMapPos=0;
#define cMperMb 0.6
	while (pos[0] != BADPOS || pos[1] != BADPOS)
	{
		if (pos[0] < pos[1])
			lowf = 0;
		else
			lowf = 1;
		mapPos=lastMapPos+(pos[lowf]-lastPos)*cMperMb/1000000.0;
		fprintf(impMap, "%010ld %f %f\n",pos[lowf],cMperMb,mapPos);
		lastPos=pos[lowf];
		lastMapPos=mapPos;
		if (lowf==1)
			pars.fromSeq[pars.nFromSeq++]=pars.nVars;
		pars.pos[pars.nVars++]=pos[lowf];
		pos[lowf]=BADPOS;
		if (fgets(linePtr[lowf],lineLength[lowf], inMap[lowf]))
			sscanf(linePtr[lowf],lowf?("%*s %*s %ld"):("%*s %*s %*s %ld"),&pos[lowf]);
	}
	for (f = 0; f < 2; ++f)
		fclose(inMap[f]);
	fclose(impMap);
	return 1;
}

int extractSNPs(masterLocusFile &vSeq, gvaParams &gp, analysisSpecs &spec, refseqGeneInfo &r, char *root, iePars &pars,char *plinkDataRoot,char *plinkQualityFlags)
{
	char line[1000];
	sprintf(line,"%s --bfile %s --recode --chr %s --from-bp %d --to-bp %d %s --out ie.%s.SNPs",
		pars.runPlink,
		plinkDataRoot,		
		r.getChr()+3, 
		r.getStart()-pars.extraBases<0?0:r.getStart()-pars.extraBases,
		r.getEnd()+pars.extraBases,
		plinkQualityFlags,
		root);
	system(line);
	return(1);
}

int getHapsWithShapeIt(masterLocusFile &vSeq,gvaParams &gp,analysisSpecs &spec,refseqGeneInfo &r,char *root,iePars &pars,int nSub)
{
	// for now, we will use orginal VCF file, just stripping out the variants which are not biallelic
	// this means we will repeate the tabix extraction
	// we still want the first one because we used it to annotate all the variants
	// we are going to assume that there is just one VCF file of "cases" (actually mixed cases and controls)

	FILE *fp,*iVCF,*oVCF;
	char line[1000],fn[100],vcfFn[100];
	float MAF;
	int onGenotypes;
	sprintf(fn,"iE.%s.extracted.vcf",root);
	sprintf(line,
				"tabix -h %s %s:%d-%d > %s",
				gp.ccFn[1][0],
				r.getChr()+(spec.addChrInVCF[0]?0:3), // leave the chr prefix if in tabix file
				r.getStart()-pars.extraBases<0?0:r.getStart()-pars.extraBases,
				r.getEnd()+pars.extraBases,
				fn);
	system(line);
	iVCF=fopen(fn,"r");
	sprintf(vcfFn,"iE.%s.forShapeIt.vcf",root);
	oVCF=fopen(vcfFn,"w");
	onGenotypes=0;
	while (fgets(longLine, LONGLINELENGTH, iVCF))
	{
		if (!onGenotypes)
		{
			if (!strncmp(longLine,"#CHROM",strlen("#CHROM")))
				onGenotypes=1;
			fprintf(oVCF,"%s",longLine);
		}
		else
		{
			sscanf(longLine,"%*s %*s %*s %*s %s",line);
			if (strchr(line,','))
				continue; // more than two alleles
			else 
				fprintf(oVCF,"%s",longLine);
		}
	}
	fclose(iVCF);
	fclose(oVCF);
	sprintf(fn,"runShapeIt.%s.%s",root,runScript==BATCH?"bat":"sh");
	fp=fopen(fn,"w");
	fprintf(fp,"%s --input-vcf %s -O iE.%s.phased\n",pars.runShapeIt,vcfFn,root);
	fclose(fp);
#ifdef MSDOS
	sprintf(line,"%s",fn);
#else
	sprintf(line,"sh %s",fn);
#endif
	system(fn);
	return 1;
}

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	iePars pars;
	char fn[100],fn2[100],line[1000],geneName[100],plinkDataRoot[100],plinkQualityFlags[100];
	int i,nSub;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;

	dcerror.warn();

	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
	// because I am lazy I will use the spec.exclusionStr lines to give me any parameters I need

	masterLocusFile vSeq(gp.nCc[0]+gp.nCc[1]);
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

	sscanf(spec.exclusionStr[0],"%s",plinkDataRoot);
	sscanf(spec.exclusionStr[1],"%d",&pars.extraBases);
	sscanf(spec.exclusionStr[2],"%f %f",&pars.minSeqFreq,&pars.maxSeqFreq);
	plinkQualityFlags[0]='\0';
	sscanf(spec.exclusionStr[3],"%[^\n]",plinkQualityFlags);
	sprintf(fn,"iEseq.%s.db",geneName);
	sprintf(fn2,"iEseq.%s.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vSeq.openFiles(fn,fn2);
	int ff=0;
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"iEseq.%s.cont.%d.vcf",geneName,i+1);
			sprintf(line,
				"tabix %s -h %s:%d-%d > %s",
				gp.ccFn[0][i],
				r.getChr()+(spec.addChrInVCF[ff++]?0:3), // leave the chr prefix if in tabix file
				r.getStart()-pars.extraBases<0?0:r.getStart()-pars.extraBases,
				r.getEnd()+pars.extraBases,
				fn);
			system(line);
			vSeq.readLocusFileEntries(fn,spec,0);
		}
		for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"iEseq.%s.case.%d.vcf",geneName,i+1);
			sprintf(line,
				"tabix -h %s %s:%d-%d > %s",
				gp.ccFn[1][i],
				r.getChr()+(spec.addChrInVCF[ff++]?0:3), // leave the chr prefix if in tabix file
				r.getStart()-pars.extraBases<0?0:r.getStart()-pars.extraBases,
				r.getEnd()+pars.extraBases,
				fn);
			system(line);
			vSeq.readLocusFileEntries(fn,spec,1);
		}
	vSeq.getQuickConsequences(r,spec);

	nSub=0;
	for (i=0;i<gp.nCc[0]+gp.nCc[1];++i)
		nSub+=vSeq.getNSubs(i);

	// getHapsWithShapeIt(vSeq,gp,spec,r,geneName,pars,nSub);
	extractSNPs(vSeq,gp,spec,r,geneName,pars,plinkDataRoot,plinkQualityFlags);
	makeFakeMap(vSeq,gp,spec,r,geneName,pars);
	//imputeSeq(vSeq,gp,spec,r,geneName,pars);  // prephase SNPs, then impute after
	//runScoreAssoc(vSeq,gp,spec,r,geneName,pars);

	// testImputation(geneName,nSub);
	
	return 0;
}

