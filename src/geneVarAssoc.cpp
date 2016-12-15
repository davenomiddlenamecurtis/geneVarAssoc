#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"
#include "vcfLocusFile.hpp"
#include "geneVarParser.hpp"

// SSS uses .:.:.:.:. to mean unknown??

#define PROGRAM "geneVarAassoc"
#define GVAVERSION "1.2"

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	int i,extractedOK;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;
	weightTable wt;
	checkSystem();

#if 1
	printf("%s %s\nRunning ",PROGRAM,GVAVERSION);
	for (i=0;i<argc;++i)
		printf("%s ",argv[i]);
	printf("\n");
#endif
	wt.init("DEFAULTWEIGHTS",consequence,NCONSEQUENCETYPES);

	dcerror.warn();
	checkSystem();
#if 0
	fp = fopen(argv[1], "r");
	if (fp == NULL)
		{ dcerror(1, "Could not open file %s", argv[1]); exit(1); }
	gp.input(fp,spec);
	fclose(fp);
	strcpy(geneName,argv[2]);
#else
	if (!gp.readParms(argc,argv,spec))
		exit(1);

	strcpy(geneName,gp.geneName);
#endif
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
	
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setBaitMargin(gp.margin);

	if (!r.findGene(geneName) || !r.getNextGene())
	{
		dcerror(1,"Could not find gene: %s\n",geneName);
		return 1;
	}
	sprintf(fn,"gva.%s.db",geneName);
	sprintf(fn2,"gva.%s.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	extractedOK=1;
	int ff=0;
	checkSystem();
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gva.%s.cont.%d.vcf",geneName,i+1);
			if (gp.dontExtractGene)
				printf("Will not attempt to produce %s because --dont-extract-gene was set",fn);
			else if (!gcont.extractGene(r,fn,0,spec.addChrInVCF[ff++]))
				extractedOK=0;
			vf.addLocusFile(fn,VCFFILE);
			if (gp.useFreqs[0])
				vf.setHoldsFreqs(0,1);
			if (!vf.readLocusFileEntries(fn,spec,0))
				extractedOK=0;
		}
	for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gva.%s.case.%d.vcf",geneName,i+1);
			if (gp.dontExtractGene)
				printf("Will not attempt to produce %s because --dont-extract-gene was set",fn);
			else if (!gcase.extractGene(r,fn,0,spec.addChrInVCF[ff++]))
				extractedOK=0;
			vf.addLocusFile(fn,VCFFILE);
			if (gp.useFreqs[1])
				vf.setHoldsFreqs(gp.nCc[1],1);
			if (!vf.readLocusFileEntries(fn,spec,1))
				extractedOK=0;
		}
	checkSystem();
	if (extractedOK)
		{
			if (spec.consequenceThreshold != 0 || spec.useConsequenceWeights != 0)
			{
				printf("Annotating...\n");
				if (spec.useEnsembl)
					vf.getEnsemblConsequences(spec);
				else
					vf.getQuickConsequences(r, spec);
			}
		}
#if 0
	if (spec.consequenceThreshold!=0)
		sprintf(fn,"gva.%s.ct%02d",geneName,spec.consequenceThreshold);
	else if (spec.useConsequenceWeights!=0)
		sprintf(fn,"gva.%s.ucw",geneName);
	else
		sprintf(fn,"gva.%s",geneName);
#else
	sprintf(fn,"%s.%s",(gp.analysisName[0]?gp.analysisName:"gva"),geneName);
#endif
	checkSystem();
	if (extractedOK)
	{
		geneVarParser::thisGene=&r; // esssential for the annotation to work
		printf("Writing scoreassoc files...\n");
		vf.writeScoreAssocFiles(fn, gp.wf,  gp.useFreqs, gp.nSubs, 1, gp.writeComments, gp.writeScoreFile, spec);
#ifndef MSDOS
		sprintf(line, "bash %s.sh\n",fn);
#else
		sprintf(line, "%s.bat\n",fn);
#endif
		if (!gp.doNotRun)
		{
		printf("Running command: %s\n", line);
		system(line);
		}
		else
			printf("Files for %s analysis written OK, to run analysis enter:\n%s\n\n",fn,line);
	}
	else
	{
		if (!gp.doNotRun)
		{
			sprintf(fn2, "%s.sao", fn);
			fp = fopen(fn2, "w");
			fprintf(fp, "Failed to extract variants for this gene\n");
			fclose(fp);
		}
		else 
			printf("Failed to extract variants for %s\n",fn);
	}
	checkSystem();
	if (!gp.keepTempFiles)
	{
		vf.closeFiles();
		vf.closeLocusFiles();
		if (!gp.dontExtractGene)
		{
			for (i = 0; i < gp.nCc[0]; ++i)
			{
				sprintf(fn, "gva.%s.cont.%d.vcf", geneName, i + 1);
				unlink(fn);
			}
			for (i = 0; i < gp.nCc[1]; ++i)
			{
				sprintf(fn, "gva.%s.case.%d.vcf", geneName, i + 1);
				unlink(fn);
			}
		}
		sprintf(fn,"gva.%s.db",geneName);
		sprintf(fn2,"gva.%s.vdx",geneName);
		unlink(fn);
		unlink(fn2);
	}
	checkSystem();
	return 0;
}