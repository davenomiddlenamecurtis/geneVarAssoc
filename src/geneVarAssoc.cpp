#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"
#include "vcfLocusFile.hpp"
#include "geneVarParser.hpp"

// SSS uses .:.:.:.:. to mean unknown??

#define PROGRAM "geneVarAassoc"
#define GVAVERSION "1.4"

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100];
	int i,extractedOK;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;
	weightTable *wt;

#if 1
	printf("%s %s\nRunning ",PROGRAM,GVAVERSION);
	for (i=0;i<argc;++i)
		printf("%s ",argv[i]);
	printf("\n");
#endif
	wt=new weightTable;
	wt->init("DEFAULTWEIGHTS",consequence,NCONSEQUENCETYPES);
	weightTableList[wt->tableName]=wt;
	wt=new weightTable;
	wt->init("DEFAULTVEPWEIGHTS",e_consequence,E_NCONSEQUENCETYPES);
	weightTableList[wt->tableName]=wt;

	dcerror.warn();
#if 0
	fp = fopen(argv[1], "r");
	if (fp == NULL)
		{ dcerror(1, "Could not open file %s", argv[1]); exit(1); }
	gp.input(fp,spec);
	fclose(fp);
	strcpy(geneName,argv[2]);
#else
	strcpy(gp.testName,"gva");
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
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gva.%s.cont.%d.vcf",geneName,i+1);
			if (gp.dontExtractVariants)
				printf("Will not attempt to produce %s because --dont-extract-variants was set",fn);
			else if (!gcont.extractVariants(r,fn,0,spec.addChrInVCF[ff++]))
				extractedOK=0;
			vf.addLocusFile(fn,VCFFILE);
			if (!vf.readLocusFileEntries(fn,spec,0))
				extractedOK=0;
		}
	for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gva.%s.case.%d.vcf",geneName,i+1);
			if (gp.dontExtractVariants)
				printf("Will not attempt to produce %s because --dont-extract-variants was set",fn);
			else if (!gcase.extractVariants(r,fn,0,spec.addChrInVCF[ff++]))
				extractedOK=0;
			vf.addLocusFile(fn,VCFFILE);
			if (!vf.readLocusFileEntries(fn,spec,1))
				extractedOK=0;
		}
	if (extractedOK)
		{
		if (spec.useEnsembl)
			spec.willNeedEnsemblConsequence=1; // I cannot think why useEnsembl would be set unless needed to get these
		if (spec.willNeedInbuiltConsequence)
		{
			if (spec.weightExpression[0]=='\0' && ! spec.willNeedEnsemblConsequence)
				strcpy(spec.weightExpression,"ANNOT(\"INBUILT\")GETWEIGHT(\"DEFAULTWEIGHTS\")");
			// default behaviour is to use these weights unless told not to
			printf("Annotating using inbuilt routines...\n");
			vf.getQuickConsequences(r,spec);
		}
		if (spec.willNeedEnsemblConsequence)
		{
			if (spec.weightExpression[0]=='\0')
				strcpy(spec.weightExpression,"ANNOT(\"VEP\")GETWEIGHT(\"DEFAULTWEIGHTS\")");
			printf("Annotating using VEP...\n");
			vf.getEnsemblConsequences(spec);
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
	sprintf(fn,"%s.%s",gp.testName,geneName);
#endif
	if (extractedOK)
	{
		geneVarParser::thisGene=&r; // esssential for the annotation to work
		printf("Writing scoreassoc files...\n");
		vf.writeScoreAssocFiles(fn, spec.wf,  gp.useFreqs, gp.nSubs, 1, gp.writeComments, gp.writeScoreFile, spec);
#ifndef MSDOS
		sprintf(line, "bash %s.sh\n",fn);
#else
		sprintf(line, "%s.bat\n",fn);
#endif
		if (!gp.doNotRun)
		{
		printf("Running command: %s\n", line);
		checkSystem();
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
	if (!gp.keepTempFiles)
	{
		vf.closeFiles();
		vf.closeLocusFiles();
		if (!gp.dontExtractVariants)
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
	return 0;
}