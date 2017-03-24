#include <stdlib.h>
#include <ctype.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"

// testimputation.cpp

typedef float genoProbTriple[3];
typedef struct oxfordLocus {
	char SNPName[100],rsName[100];
	long pos;
	char refAll,altAll;
} ;
#define MAXLOCIINANALYSIS 1000

float getAltFreq(allelePair *all, int nSub)
{
	int nKnown,nAlt,s;
	nAlt=nKnown=0;
	for (s=0;s<nSub;++s)
	if (all[s][0] != 0)
	{
		nAlt += (all[s][0] != 1) + (all[s][1] != 1);
		++nKnown;
	}
	return nAlt*0.5/nKnown;
}

#define LONGLINELENGTH 40000
#define BADPOS 1000000000L
char longLineM[LONGLINELENGTH+1],longLineT[LONGLINELENGTH+1];
int testImputation(char *root, int nSub)
{
	FILE *fgens[2],*fref,*ftest,*fmap;
	char fn[100],*sptr[2],*ptr,line[1000];
	long pos[2],lastPos,firstPos;
	int f,lowf,nGeno,l,testLoci[MAXLOCIINANALYSIS],nLoci,nTest,w,s,t,ss;
	genoProbTriple **geno;
	oxfordLocus *loci;
	float mapPos,lastMapPos;
	geno=(genoProbTriple **)calloc(MAXLOCIINANALYSIS,sizeof(genoProbTriple *));
	loci=(oxfordLocus *)calloc(MAXLOCIINANALYSIS,sizeof(oxfordLocus));
	sprintf(fn,"%s.marker.gens",root);
	fgens[0]=fopen(fn,"r");
	sprintf(fn,"%s.test.gens",root);
	fgens[1]=fopen(fn,"r");
	sptr[0]=longLineM;
	sptr[1]=longLineT;
	nLoci=nTest=0;
	pos[0]=pos[1]=BADPOS;
	for (f = 0; f < 2; ++f)
	{
		fgets(sptr[f],LONGLINELENGTH,fgens[f]);
		sscanf(sptr[f],"%*s %*s %ld",&pos[f]);
	}
	while (pos[0] != BADPOS || pos[1] != BADPOS)
	{
		if (pos[0]<pos[1])
			lowf=0;
		else
			lowf=1;
		sscanf(sptr[lowf],"%s %s %ld %c %c ",loci[nLoci].SNPName,loci[nLoci].rsName,&loci[nLoci].pos,&loci[nLoci].refAll,&loci[nLoci].altAll);
		ptr=sptr[lowf];
		for (w = 0; w < 5; ++w)
		{
			while (!isspace(*ptr)) ++ptr;
			while (*ptr && isspace(*ptr)) ++ptr;
		}
		geno[nLoci]=(genoProbTriple*)calloc(nSub,sizeof(genoProbTriple));
		for (s = 0; s < nSub; ++s)
		{
			sscanf(ptr,"%f %f %f",&geno[nLoci][s][0],&geno[nLoci][s][1],&geno[nLoci][s][2]);
			for (w = 0; w < 3; ++w)
			{
				while (!isspace(*ptr)) ++ptr;
				while (*ptr && isspace(*ptr)) ++ptr;
			}
		}
		if (lowf==1)
			testLoci[nTest++]=nLoci;
		pos[lowf]=BADPOS;
		if (fgets(sptr[lowf],LONGLINELENGTH,fgens[lowf]))
			sscanf(sptr[lowf],"%*s %*s %ld",&pos[lowf]);
		++nLoci;
	}
	for (f=0;f<2;++f)
		fclose(fgens[f]);
	sprintf(fn,"%s.map",root);
	fmap=fopen(fn,"w");
#define cMperMb 0.6
	lastPos=firstPos=loci[0].pos;
	lastMapPos=0;
	fprintf(fmap,"position	COMBINED_rate.cM.Mb.	Genetic_Map.cM.\n");
	for (l = 0; l < nLoci; ++l)
	{
		mapPos=lastMapPos+(loci[l].pos-lastPos)*cMperMb/1000000.0;
		fprintf(fmap, "%010ld %f %f\n",loci[l].pos,cMperMb,mapPos);
		lastPos=loci[l].pos;
		lastMapPos=mapPos;
	}
	fclose(fmap);
	for (s = 0; s < nSub; ++s)
	{
		if (geno[testLoci[2]][s][1]==0)
			continue;
		sprintf(fn,"%s.testone.gens",root);
		ftest=fopen(fn,"w");
		for (l = 0; l < nLoci; ++l)
		{
			for (t=0;t<nTest;++t)
				if (testLoci[t]==l)
					break;
			if (t!=nTest)
				continue; // skip all test loci
			fprintf(ftest,"%s %s %010ld %c %c %.0f %.0f %.0f\n",
				loci[l].SNPName,loci[l].rsName,loci[l].pos,loci[l].refAll,loci[l].altAll,
				geno[l][s][0],geno[l][s][1],geno[l][s][2]);
		}
		fclose(ftest);
		sprintf(fn,"%s.refone.gens",root);
		fref=fopen(fn,"w");
		for (l = 0; l < nLoci; ++l)
		{
			fprintf(fref, "%s %s %010ld %c %c ", loci[l].SNPName, loci[l].rsName, loci[l].pos,loci[l].refAll,loci[l].altAll);
			for (ss = 0; ss < nSub; ++ss)
				if (ss != s)
					fprintf(fref, "%.0f %.0f %.0f  ",geno[l][ss][0], geno[l][ss][1], geno[l][ss][2]);
			fprintf(fref,"\n");
		}
		fclose(fref);
		sprintf(line,
"\\impute2\\impute2 \
 -m %s.map \
 -g_ref %s.refone.gens \
 -g %s.testone.gens \
 -int %ld %ld  \
 -Ne 20000 \
 -o %s.impute2.out",root,root,root,firstPos-1L,lastPos+1L,root);

		system(line);
		return 1;
	}

	for (l=0;l<nLoci;++l)
		free(geno[l]);
	free(geno);
	free(loci);
	return 1;
}

int writeImputeGenoFiles(masterLocusFile &vTest, masterLocusFile &vSNP,analysisSpecs spec,char *root,int nSub,float minTestAF, float maxTestAF, float minSNPAF, float maxSNPAF)
{
	FILE *fp;
	char line[1000];
	float MAF;
	allelePair *alls;
	int s,loc=1;
	alls=(allelePair *)calloc(nSub,sizeof(allelePair));

	sprintf(line,"%s.test.gens",root);
	fp=fopen(line,"w");
	// in fact, I will not be using ranges in spec but assume vTest and vSNP have only the correct ranges in them
	vTest.openLocusFiles();
	if (!vTest.gotoFirstInRange(spec))
		return 0;
	do {
		if (spec.onlyUseSNPs && !vTest.currentIsSNP())
			continue;
		if (vTest.currentWorstConsequenceType()<spec.consequenceThreshold)
			continue;
		vTest.outputCurrentAlleles(alls,spec);
		MAF=getAltFreq(alls,nSub);
		if (MAF>0.5)
			MAF=1.0-MAF;
		if (MAF<minTestAF || MAF >maxTestAF)
			continue;
		fprintf(fp,"SNP%05d rs%05d %010ld C T  ",loc,loc,vTest.currentPos());
		for (s=0;s<nSub;++s)
			fprintf(fp,"%s  ",
				alls[s][0]==0?"0 0 0":
				alls[s][0]==1 && alls[s][1]==1 ?"0 0 1":
				alls[s][0]!=1 && alls[s][1]!=1 ?"1 0 0":
				"0 1 0");
		fprintf(fp,"\n");
		++loc;
	} while (vTest.gotoNextInRange(spec));
	fclose(fp);

	// we have not reset loc so names will be distinct
	sprintf(line,"%s.marker.gens",root);
	fp=fopen(line,"w");
	vSNP.openLocusFiles();
	if (!vSNP.gotoFirstInRange(spec))
		return 0;
	do {
		if (spec.onlyUseSNPs && !vSNP.currentIsSNP())
			continue;
		vSNP.outputCurrentAlleles(alls,spec);
		MAF=getAltFreq(alls,nSub);
		if (MAF>0.5)
			MAF=1.0-MAF;
		if (MAF<minSNPAF || MAF >maxSNPAF)
			continue;
		// we assume SNP and test allele frequency ranges do not overlap
		fprintf(fp,"SNP%05d rs%05d %010ld C T  ",loc,loc,vSNP.currentPos());
		for (s=0;s<nSub;++s)
			fprintf(fp,"%s  ",
				alls[s][0]==0?"0 0 0":
				alls[s][0]==1 && alls[s][1]==1 ?"0 0 1":
				alls[s][0]!=1 && alls[s][1]!=1 ?"1 0 0":
				"0 1 0");
		fprintf(fp,"\n");
		++loc;
	} while (vSNP.gotoNextInRange(spec));
	fclose(fp);

	free(alls);
	return 1;
}

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	int i,margin,nSub;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;
	float minTestAF,maxTestAF,minSNPAF,maxSNPAF;

	dcerror.warn();

	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
	// because I am lazy I will use the spec.exclusionStr lines to give me any parameters I need

	masterLocusFile vTest(gp.nCc[0]+gp.nCc[1]),vSNP(gp.nCc[0]+gp.nCc[1]);
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

	sprintf(fn,"tITest.%s.db",geneName);
	sprintf(fn2,"tITest.%s.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vTest.openFiles(fn,fn2);
	int ff=0;
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"tITest.%s.cont.%d.vcf",geneName,i+1);
			gcont.extractGene(r,fn,0,spec.addChrInVCF[ff++]);
			vTest.readLocusFileEntries(fn,spec,0);
		}
		for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"tITest.%s.case.%d.vcf",geneName,i+1);
			gcase.extractGene(r,fn,0,spec.addChrInVCF[ff++]);
			vTest.readLocusFileEntries(fn,spec,1);
		}
	vTest.getQuickConsequences(r,spec);

	sscanf(spec.exclusionStr[0],"%d",&margin);
	sprintf(fn,"tISNP.%s.db",geneName);
	sprintf(fn2,"tISNP.%s.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vSNP.openFiles(fn,fn2);
	int ff=0;
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"tISNP.%s.cont.%d.vcf",geneName,i+1);
			sprintf(line,
				"tabix %s -h %s:%d-%d > %s",
				gp.ccFn[0][i],
				r.getChr()+(spec.addChrInVCF[ff++]?0:3), // leave the chr prefix if in tabix file
				r.getStart()-margin<0?0:r.getStart()-margin,
				r.getEnd()+margin,
				fn);
			system(line);
			vSNP.readLocusFileEntries(fn,spec,0);
		}
		for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"tISNP.%s.case.%d.vcf",geneName,i+1);
			sprintf(line,
				"tabix -h %s %s:%d-%d > %s",
				gp.ccFn[1][i],
				r.getChr()+(spec.addChrInVCF[ff++]?0:3), // leave the chr prefix if in tabix file
				r.getStart()-margin<0?0:r.getStart()-margin,
				r.getEnd()+margin,
				fn);
			system(line);
			vSNP.readLocusFileEntries(fn,spec,1);
		}

	nSub=0;
	for (i=0;i<gp.nCc[0]+gp.nCc[1];++i)
		nSub+=vSNP.getNSubs(i);

	sscanf(spec.exclusionStr[1],"%f %f",&minTestAF,&maxTestAF);
	sscanf(spec.exclusionStr[2],"%f %f",&minSNPAF,&maxSNPAF);
	
	writeImputeGenoFiles(vTest,vSNP,spec,geneName,nSub,minTestAF,maxTestAF,minSNPAF,maxSNPAF);
	testImputation(geneName,nSub);
	
	return 0;
}