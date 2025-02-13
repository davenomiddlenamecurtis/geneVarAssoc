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
#define IVAVERSION "1.5"

gvaParams gp;
analysisSpecs spec;
refseqGeneInfo r; // even though there is no gene will use this to access reference files to e.g. check for CpG changes

int main(int argc,char *argv[])
{
	char fn[100],fn2[100],line[1000], chrStr[100], startStr[100], endStr[100], chr[10], * tbiFn, * ptr, tbiFnBuff[1000];
	int i,first,extractedOK,cc;
	FILE *fp,*fi;
	weightTable *wt;
	intervalList iList;

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
	masterLocusFile vf(gp.bedFileFn[0] ? 1 : gp.nCc[0] + gp.nCc[1]);
	if (gp.referencePath[0] != '\0')
		r.setReferencePath(gp.referencePath);

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
	fi = fopen(gp.intervalListFn, "r");
	if (!fi)
	{
		dcerror(1, "Could not open interval list file %s\n", gp.intervalListFn);
		exit(1);
	}
	sprintf(fn,"iva.%s.db",gp.testName);
	sprintf(fn2,"iva.%s.vdx",gp.testName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	extractedOK = 1;
	int ff = 0;
	while (fgets(line, 999, fi) && sscanf(line, " %[^:]:%[^-]-%s", chrStr, startStr,endStr) == 3)
			iList.append(chrStr, atoi(startStr), atoi(endStr));
	if (gp.bedFileFn[0])
		{
			if (gp.nCc[0] || gp.nCc[1])
			{
				dcerror(2, "Should not use --bed-file %s if also using --case-file or --cont-file.\n", gp.bedFileFn);
				return 1;
			}
			strcpy(gp.ccFn[0][gp.nCc[0]++], gp.bedFileFn); // pretend bedFile is a control file for now
		}
	extractedOK = 1;
	for (i = 0; i < gp.nCc[0]; ++i)
	{
		sprintf(fn, "%s.cont.%d.vcf", gp.testName, i + 1);
		if (gp.dontExtractVariants)
			printf("Will not attempt to produce %s because --dont-extract-variants was set\n", fn);
		else
			//if (!gcont.extractVariants(r,fn,0,spec.addChrInVCF[ff++],spec.removeVcfSpaces,spec.omitIntrons, spec.spliceRegionSize))
			//extractedOK=0;
		{
			if (gp.bedFileFn[0] == '\0')
			{
				//					if (!r.tbiExtractGene(gp.ccFn[0][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces, spec.omitIntrons, spec.spliceRegionSize))
				if (!tbiExtractIntervals(gp.ccFn[0][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces, iList))
					extractedOK = 0;
			}
			else
			{
				//					if (!r.plinkExtractGene(gp.bedFileFn,gp.famFileFn,gp.bimFileFn, fn, spec.omitIntrons, spec.spliceRegionSize))
				if (!plinkExtractIntervals(gp.bedFileFn, gp.famFileFn, gp.bimFileFn, fn, iList, gp.testName))
					extractedOK = 0;

			}
		}
		vf.addLocusFile(fn, VCFFILE);
		if (!vf.readLocusFileEntries(fn, spec, 0))
			extractedOK = 0;
	}
	for (i = 0; i < gp.nCc[1]; ++i)
	{
		sprintf(fn, "%s.case.%d.vcf", gp.testName, i + 1);
		if (gp.dontExtractVariants)
			printf("Will not attempt to produce %s because --dont-extract-variants was set\n", fn);
		else // if (!gcase.extractVariants(r,fn,0,spec.addChrInVCF[ff++], spec.removeVcfSpaces,spec.omitIntrons,spec.spliceRegionSize))
			// extractedOK=0;
//			if (!r.tbiExtractGene(gp.ccFn[1][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces, spec.omitIntrons, spec.spliceRegionSize))
		{
			if (!tbiExtractIntervals(gp.ccFn[1][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces, iList))
				extractedOK = 0;
		}
		vf.addLocusFile(fn, VCFFILE);
		if (!vf.readLocusFileEntries(fn, spec, 1))
			extractedOK = 0;
	}
	if (gp.onlyExtractVariants)
	{
		printf("Exiting because --only-extract-variants was set\n");
		exit(0);
	}
	if (spec.willNeedEnsemblConsequence)
	{
		printf("Annotating using VEP...\n");
		vf.getEnsemblConsequences(spec);
	}

	sprintf(fn, "%s", gp.testName);
	if (extractedOK)
	{
		geneVarParser::thisGene = &r; // in case needed e.g. for CpGs
		printf("Writing scoreassoc files...\n");
		vf.writeScoreAssocFiles(fn, spec.wf, gp.useFreqs, gp.nSubs, 1, gp.writeComments, gp.writeScoreFile, gp.writeRecScoreFile, spec);
#ifndef MSDOS
		sprintf(line, "bash %s.sh\n", fn);
#else
		sprintf(line, "%s.bat\n", fn);
#endif
		if (!gp.doNotRun)
		{
			printf("Running command: %s\n", line);
			checkSystem();
			system(line);
		}
		else
			printf("Files for %s analysis written OK, to run analysis enter:\n%s\n\n", fn, line);
	}
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
		sprintf(fn,"iva.%s.db",gp.testName);
		sprintf(fn2,"iva.%s.vdx",gp.testName);
		unlink(fn);
		unlink(fn2);
	}

	return 0;
}