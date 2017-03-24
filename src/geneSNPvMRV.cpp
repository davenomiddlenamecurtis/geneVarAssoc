// geneSNPvMRV
// run ShowLD to see LD relationships between MRVs in one gene and SNPs in those nearby
#if 0
example input

ank3.SNPs.txt - output file
5 10 number of MRVs number of samples/trials
\sequence\sharedseq\gva.uk.ob.MRV.130420.par 
ank3
\sequence\sharedseq\gva.uk.ob.SNP.130420.par 
ank3
cdk1
arid5b

will produce command line:
showLD gSM.ANK3.MRVs.par gSM.ANK3.MRVs.dat gSM.ANK3.SNPs.par gSM.ANK3.SNPs.dat ANK3.SNPs.txt 5 10


#endif

#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"


	int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100],markerGeneName[100],ofn[100];
	int i,nMRV,nTrial,before,after,st,en;
	FILE *fp,*fpSpec;
	gvaParams gpMRVs,gpSNPs;
	analysisSpecs specMRVs,specSNPs;

	dcerror.warn();

	fp=fopen(argv[1],"r");
	fgets(line,999,fp);
	sscanf(line,"%s",ofn);
	fgets(line,999,fp);
	sscanf(line,"%d %d",&nMRV,&nTrial);
	fgets(line,999,fp);
	sscanf(line,"%s",fn);
	fpSpec=fopen(fn,"r");
	gpMRVs.input(fpSpec,specMRVs);
	fclose(fpSpec);

	masterLocusFile vf(gpMRVs.nCc[0]+gpMRVs.nCc[1]);

	r.setListFile(gpMRVs.geneListFn);
	r.setBaitsFile(gpMRVs.baitFn);
	if (gpMRVs.referencePath[0]!='\0')
		r.setReferencePath(gpMRVs.referencePath);
	r.setUpstream(gpMRVs.upstream);
	r.setDownstream(gpMRVs.downstream);
	r.setBaitMargin(gpMRVs.margin);

	fgets(line,999,fp);
	sscanf(line,"%s",geneName);

	if (!r.findGene(geneName) || !r.getNextGene())
		return 1;

	sprintf(fn,"gSM.%s.MRVs.db",geneName);
	sprintf(fn2,"gSM.%s.MRVs.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	for (i=0;i<gpMRVs.nCc[0];++i)
		{
			gcont.setVariantFileName(gpMRVs.ccFn[0][i]);
			sprintf(fn,"gSM.%s.cont.%d.vcf",geneName,i+1);
			gcont.extractGene(r,fn);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,specMRVs,0);
		}
	for (i=0;i<gpMRVs.nCc[1];++i)
		{
			gcase.setVariantFileName(gpMRVs.ccFn[1][i]);
			sprintf(fn,"gSM.%s.case.%d.vcf",geneName,i+1);
			gcase.extractGene(r,fn);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,specMRVs,1);
		}
	vf.getQuickConsequences(r,specMRVs);
	sprintf(fn,"gva.%s.MRVs",geneName);
	vf.writeOldScoreAssocFiles(fn,gpMRVs.wf,gpMRVs.wFunc,gpMRVs.useFreqs,gpMRVs.nSubs,1,gpMRVs.writeComments,specMRVs);

	fgets(line,999,fp);
	sscanf(line,"%s",fn);
	fpSpec=fopen(fn,"r");
	gpSNPs.input(fpSpec,specSNPs);
	fclose(fpSpec);

	masterLocusFile vfSNPs(gpSNPs.nCc[0]+gpSNPs.nCc[1]);
	sprintf(fn,"gSM.%s.SNPs.db",geneName);
	sprintf(fn2,"gSM.%s.SNPs.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vfSNPs.openFiles(fn,fn2);

	for (i=0;i<gpSNPs.nCc[0];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.cont.%d.vcf",geneName,i+1);
		unlink(fn);
	}
	for (i=0;i<gpSNPs.nCc[1];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.case.%d.vcf",geneName,i+1);
		unlink(fn);
	}

	fgets(line,999,fp);
	sscanf(line,"%d %d",&before,&after);
	st=r.getStart()-before;
	en=r.getStart()+after;

		for (i=0;i<gpSNPs.nCc[0];++i)
		{
			// gcont.setVariantFileName(gpSNPs.ccFn[0][i]);
			sprintf(fn,"gSM.%s.SNPs.cont.%d.vcf",geneName,i+1);
			sprintf(line,"tabix %s -h %s:%d-%d > %s",gpSNPs.ccFn[0][i],r.getChr()+3,st,en,fn);
			system(line);
			// gcont.extractGene(r,fn,(first==1)?0:1); // 1 means append to old
		}
		for (i=0;i<gpSNPs.nCc[1];++i)
		{
			// gcase.setVariantFileName(gpSNPs.ccFn[1][i]);
			sprintf(fn,"gSM.%s.SNPs.case.%d.vcf",geneName,i+1);
			sprintf(line,"tabix %s -h %s:%d-%d > %s",gpSNPs.ccFn[1][i],r.getChr()+3,st,en,fn);
			system(line);
			// gcase.extractGene(r,fn,(first==1)?0:1);
		}

	for (i=0;i<gpSNPs.nCc[0];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.cont.%d.vcf",geneName,i+1);
		vfSNPs.readLocusFileEntries(fn,specSNPs,0);
	}
	for (i=0;i<gpSNPs.nCc[1];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.case.%d.vcf",geneName,i+1);
		vfSNPs.readLocusFileEntries(fn,specSNPs,1);
	}

	vfSNPs.getQuickConsequences(r,specSNPs);
	sprintf(fn,"gva.%s.SNPs",geneName);
	vfSNPs.writeOldScoreAssocFiles(fn,gpSNPs.wf,gpSNPs.wFunc,gpSNPs.useFreqs,gpSNPs.nSubs,1,gpSNPs.writeComments,specSNPs);

	sprintf(line,"showLD gva.%s.MRVs.par gva.%s.MRVs.dat gva.%s.SNPs.par gva.%s.SNPs.dat %s %d %d 1 >exclusions.txt",
		geneName,geneName,geneName,geneName,ofn,nMRV,nTrial);
	system(line);

}

#if 0
int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100],markerGeneName[100],ofn[100];
	int i,nMRV,nTrial,first;
	FILE *fp,*fpSpec;
	gvaParams gpMRVs,gpSNPs;
	analysisSpecs specMRVs,specSNPs;

	dcerror.warn();

	fp=fopen(argv[1],"r");
	fgets(line,999,fp);
	sscanf(line,"%s",ofn);
	fgets(line,999,fp);
	sscanf(line,"%d %d",&nMRV,&nTrial);
	fgets(line,999,fp);
	sscanf(line,"%s",fn);
	fpSpec=fopen(fn,"r");
	gpMRVs.input(fpSpec,specMRVs);
	fclose(fpSpec);

	masterLocusFile vf(gpMRVs.nCc[0]+gpMRVs.nCc[1]);

	r.setListFile(gpMRVs.geneListFn);
	r.setBaitsFile(gpMRVs.baitFn);
	if (gpMRVs.referencePath[0]!='\0')
		r.setReferencePath(gpMRVs.referencePath);
	r.setUpstream(gpMRVs.upstream);
	r.setDownstream(gpMRVs.downstream);
	r.setBaitMargin(gpMRVs.margin);

	fgets(line,999,fp);
	sscanf(line,"%s",geneName);

	if (!r.findGene(geneName) || !r.getNextGene())
		return 1;

	sprintf(fn,"gSM.%s.MRVs.db",geneName);
	sprintf(fn2,"gSM.%s.MRVs.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
		for (i=0;i<gpMRVs.nCc[0];++i)
		{
			gcont.setVariantFileName(gpMRVs.ccFn[0][i]);
			sprintf(fn,"gSM.%s.cont.%d.vcf",geneName,i+1);
			gcont.extractGene(r,fn);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,specMRVs,0);
		}
		for (i=0;i<gpMRVs.nCc[1];++i)
		{
			gcase.setVariantFileName(gpMRVs.ccFn[1][i]);
			sprintf(fn,"gSM.%s.case.%d.vcf",geneName,i+1);
			gcase.extractGene(r,fn);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,specMRVs,1);
		}
	vf.getQuickConsequences(r,specMRVs);
	sprintf(fn,"%s.MRVs",geneName);
	vf.writeOldScoreAssocFiles(fn,gpMRVs.wf,gpMRVs.wFunc,gpMRVs.useFreqs,gpMRVs.nSubs,1,gpMRVs.writeComments,specMRVs);

	fgets(line,999,fp);
	sscanf(line,"%s",fn);
	fpSpec=fopen(fn,"r");
	gpSNPs.input(fpSpec,specSNPs);
	fclose(fpSpec);

	masterLocusFile vfSNPs(gpSNPs.nCc[0]+gpSNPs.nCc[1]);
	sprintf(fn,"gSM.%s.SNPs.db",geneName);
	sprintf(fn2,"gSM.%s.SNPs.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vfSNPs.openFiles(fn,fn2);

	for (i=0;i<gpSNPs.nCc[0];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.cont.%d.vcf",geneName,i+1);
		unlink(fn);
	}
	for (i=0;i<gpSNPs.nCc[1];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.case.%d.vcf",geneName,i+1);
		unlink(fn);
	}

	first=1;
	while (fgets(line,999,fp) && sscanf(line,"%s",markerGeneName)==1)
	{
		if (!r.findGene(markerGeneName) || !r.getNextGene())
			return 1;

		for (i=0;i<gpSNPs.nCc[0];++i)
		{
			gcont.setVariantFileName(gpSNPs.ccFn[0][i]);
			sprintf(fn,"gSM.%s.SNPs.cont.%d.vcf",geneName,i+1);
			gcont.extractGene(r,fn,(first==1)?0:1); // 1 means append to old
		}
		for (i=0;i<gpSNPs.nCc[1];++i)
		{
			gcase.setVariantFileName(gpSNPs.ccFn[1][i]);
			sprintf(fn,"gSM.%s.SNPs.case.%d.vcf",geneName,i+1);
			gcase.extractGene(r,fn,(first==1)?0:1);
		}
		first=0;
	}

	for (i=0;i<gpSNPs.nCc[0];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.cont.%d.vcf",geneName,i+1);
		vfSNPs.readLocusFileEntries(fn,specSNPs,0);
	}
	for (i=0;i<gpSNPs.nCc[1];++i)
	{
		sprintf(fn,"gSM.%s.SNPs.case.%d.vcf",geneName,i+1);
		vfSNPs.readLocusFileEntries(fn,specSNPs,1);
	}

	vfSNPs.getQuickConsequences(r,specSNPs);
	sprintf(fn,"%s.SNPs",geneName);
	vfSNPs.writeOldScoreAssocFiles(fn,gpSNPs.wf,gpSNPs.wFunc,gpSNPs.useFreqs,gpSNPs.nSubs,1,gpSNPs.writeComments,specSNPs);

	sprintf(line,"showLD gva.%s.MRVs.par gva.%s.MRVs.dat gva.%s.SNPs.par gva.%s.SNPs.dat %s %d %d 1 >exclusions.txt",
		geneName,geneName,geneName,geneName,ofn,nMRV,nTrial);
	system(line);

}
#endif
