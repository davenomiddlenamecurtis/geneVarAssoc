#if 0
Copyright 2018 David Curtis

This file is part of the geneVarAssoc package.

geneVarAssoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

geneVarAssoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with geneVarAssoc.If not, see <http://www.gnu.org/licenses/>.
#endif

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
	char fn[100],fn2[100],line[1000],intStr[100],chr[10], * tbiFn, * ptr, tbiFnBuff[1000];
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
		printf("Will not attempt to produce %s because --dont-extract-variants was set\n",fn);
	else
		while (fgets(line,999,fi) && sscanf(line,"%s",intStr)==1)
	{
		printf("%s\n",intStr);
		int ff=0;
		for (i=0;i<gp.nCc[0];++i)
		{
			if ((ptr = strchr(gp.ccFn[0][i], '*')) == 0)
				tbiFn = gp.ccFn[0][i];
			else
			{
				sscanf(intStr, "%[^:]", chr);
				strcpy(tbiFnBuff, gp.ccFn[0][i]);
				ptr = strchr(tbiFnBuff, '*');
				*ptr = '\0';
				strcat(tbiFnBuff,chr);
				ptr = strchr(gp.ccFn[0][i], '*');
				strcat(tbiFnBuff, ptr + 1);
				tbiFn = tbiFnBuff;
			}
			sprintf(fn,"%s.cont.%d.vcf",gp.testName,i+1);
			sprintf(line,"tabix %s %s %s%s %s %s", tbiFn,first==1?"-h":"",spec.addChrInVCF[ff++]?"chr":"",intStr,first?">":">>",fn);
			checkSystem();
			printf("Executing command:\n%s\n", line);
			system(line);
		}
		for (i=0;i<gp.nCc[1];++i)
		{
			if ((ptr = strchr(gp.ccFn[1][i], '*')) == 0)
				tbiFn = gp.ccFn[1][i];
			else
			{
				sscanf(intStr, "%[^:]", chr);
				strcpy(tbiFnBuff, gp.ccFn[1][i]);
				ptr = strchr(tbiFnBuff, '*');
				*ptr = '\0';
				strcat(tbiFnBuff, chr);
				ptr = strchr(gp.ccFn[1][i], '*');
				strcat(tbiFnBuff, ptr + 1);
				tbiFn = tbiFnBuff;
			}

			sprintf(fn,"%s.case.%d.vcf",gp.testName,i+1);
			sprintf(line,"tabix %s %s %s%s %s %s",tbiFn,first==1?"-h":"",spec.addChrInVCF[ff++]?"chr":"",intStr,first?">":">>",fn);
			checkSystem();
			printf("Executing command:\n%s\n", line);
			system(line);
		}
		first=0;
	}
	if (gp.onlyExtractVariants)
	{
		printf("Exiting because --only-extract-variants was set\n");
		exit(0);
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