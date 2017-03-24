#include <stdlib.h>
#ifndef MSDOS 
#include <unistd.h>
#endif
#include "geneVarUtils.hpp"
#include "vcfLocusFile.hpp"
#include "geneVarParser.hpp"

#define PROGRAM "intVarAassoc"
#define IVAVERSION "1.4"

int main(int argc,char *argv[])
{
	char fn[100],fn2[100],line[1000],intStr[100];
	int i,first,extractedOK;
	FILE *fp,*fi;
	gvaParams gp;
	analysisSpecs spec;
	weightTable *wt;

	printf("%s %s\nRunning ",PROGRAM,IVAVERSION);
	for (i=0;i<argc;++i)
		printf("%s ",argv[i]);
	printf("\n");

	wt=new weightTable;
	wt->init("DEFAULTWEIGHTS",consequence,NCONSEQUENCETYPES);
	weightTableList[wt->tableName]=wt;
	wt=new weightTable;
	wt->init("DEFAULTVEPWEIGHTS",e_consequence,E_NCONSEQUENCETYPES);
	weightTableList[wt->tableName]=wt;

	dcerror.warn();

	strcpy(gp.testName,"iva");
	if (!gp.readParms(argc,argv,spec))
		exit(1);
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);

	if (gp.intervalListFn[0]=='\0')
	{
		dcerror(1,"Must set --interval-list-file\n");
		exit(1);
	}
	if (spec.willNeedInbuiltConsequence)
	{
		dcerror(1,"Cannot use inbuilt annotation routines with intVarAssoc\n");
		exit(1);
	}
	fi=fopen(gp.intervalListFn,"r");
	sprintf(fn,"gva.%s.db",gp.testName);
	sprintf(fn2,"gva.%s.vdx",gp.testName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	first=1;
	if (gp.dontExtractVariants)
		printf("Will not attempt to produce %s because --dont-extract-variants was set",fn);
	else
		while (fgets(line,999,fi) && sscanf(line,"%s",intStr)==1)
	{
		printf("%s\n",intStr);
		int ff=0;
		for (i=0;i<gp.nCc[0];++i)
		{
			sprintf(fn,"%s.cont.%d.vcf",gp.testName,i+1);
			sprintf(line,"tabix %s %s %s%s %s %s",gp.ccFn[0][i],first==1?"-h":"",spec.addChrInVCF[ff++]?"chr":"",intStr,first?">":">>",fn);
			checkSystem();
			system(line);
		}
		for (i=0;i<gp.nCc[1];++i)
		{
			sprintf(fn,"%s.case.%d.vcf",gp.testName,i+1);
			sprintf(line,"tabix %s %s %s%s %s %s",gp.ccFn[1][i],first==1?"-h":"",spec.addChrInVCF[ff++]?"chr":"",intStr,first?">":">>",fn);
			checkSystem();
			system(line);
		}
		first=0;
	}
	for (i=0;i<gp.nCc[0];++i)
		{
			sprintf(fn,"%s.cont.%d.vcf",gp.testName,i+1);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,0);
		}
	for (i=0;i<gp.nCc[1];++i)
		{
			sprintf(fn,"%s.case.%d.vcf",gp.testName,i+1);
			vf.addLocusFile(fn,VCFFILE);
			vf.readLocusFileEntries(fn,spec,1);
		}
	if (spec.willNeedEnsemblConsequence)
	{
		printf("Annotating using VEP...\n");
		vf.getEnsemblConsequences(spec);
	}

	printf("Writing scoreassoc files...\n");
	vf.writeScoreAssocFiles(fn,spec.wf,gp.useFreqs,gp.nSubs,1,gp.writeComments,gp.writeScoreFile,spec);
#ifndef MSDOS
	sprintf(line,"bash %s.sh\n",fn);
#else
	sprintf(line,"%s.bat\n",fn);
#endif
	if(!gp.doNotRun)
	{
		printf("Running command: %s\n",line);
		checkSystem();
		system(line);
	}
	else
		printf("Files for %s analysis written OK, to run analysis enter:\n%s\n\n",fn,line);
	if (!gp.keepTempFiles)
	{
		vf.closeFiles();
		vf.closeLocusFiles();
		if (!gp.dontExtractVariants)
		{
			for (i = 0; i < gp.nCc[0]; ++i)
			{
				sprintf(fn,"%s.cont.%d.vcf",gp.testName,i + 1);
				unlink(fn);
			}
			for (i = 0; i < gp.nCc[1]; ++i)
			{
				sprintf(fn,"%s.case.%d.vcf",gp.testName,i + 1);
				unlink(fn);
			}
		}
		sprintf(fn,"%s.db",gp.testName);
		sprintf(fn2,"%s.vdx",gp.testName);
		unlink(fn);
		unlink(fn2);
	}

	return 0;
}