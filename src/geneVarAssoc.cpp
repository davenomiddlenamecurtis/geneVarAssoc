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

// SSS uses .:.:.:.:. to mean unknown??

#define PROGRAM "geneVarAassoc"
#define GVAVERSION "7.1"

refseqGeneInfo r; // moved all of these off stack
analysisSpecs spec;
gvaParams gp;

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	char fn[100], fn2[100], saoFn[100], line[1000],geneName[100];
	int i,extractedOK;
	FILE *fp;
	weightTable *wt;
	intervalList iList;

	printf("%s %s\nRunning ",PROGRAM,GVAVERSION);
	for (i=0;i<argc;++i)
		printf("%s ",argv[i]);
	printf("\n");
	printf("BUFFSIZE is %d\n", BUFFSIZE);
	printf("MAXSCOREASSOCARGS is %d\n", MAXSCOREASSOCARGS);

	wt=new weightTable;
	wt->init("DEFAULTWEIGHTS",consequence,NCONSEQUENCETYPES);
	weightTableList[wt->tableName]=wt;
	wt=new weightTable;
	wt->init("DEFAULTVEPWEIGHTS",e_consequence,E_NCONSEQUENCETYPES);
	weightTableList[wt->tableName]=wt;

	dcerror.warn();
	strcpy(gp.testName,"gva");
	if (!gp.readParms(argc,argv,spec))
		exit(1);

	strcpy(geneName,gp.geneName);
	masterLocusFile vf(gp.bedFileFn[0]?1: gp.nCc[0]+gp.nCc[1]);
	
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setBaitMargin(gp.margin);
	sprintf(saoFn, "%s.%s.sao", gp.testName, geneName);
	if (!r.findGene(geneName) || !r.getNextGene())
	{
		if (!gp.doNotRun)
		{
			fp = fopen(saoFn, "w");
			fprintf(fp, "No valid lines for this gene\n");
			fclose(fp);
		}
		else
			printf("No valid lines for %s\n", fn);
		dcerror(1, "Could not find gene: %s\n", geneName);
		return 1;
	}
	if (!gp.dontExtractVariants)
	{
		r.getGeneIntervals(iList, spec.omitIntrons, spec.spliceRegionSize);
		if (iList.nInts == 0)
		{
			if (!gp.doNotRun)
			{
				fp = fopen(saoFn, "w");
				fprintf(fp, "No valid intervals for this gene\n");
				fclose(fp);
			}
			else
				printf("No valid intervals for %s\n", fn);
			dcerror(1, "Could not find any valid intervals for gene: %s\n", geneName);
			return 1;
		}
	}
	sprintf(fn,"gva.%s.db",geneName);
	sprintf(fn2,"gva.%s.vdx",geneName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
	extractedOK=1;
	int ff=0;
	if (gp.bedFileFn[0])
	{
		if (gp.nCc[0] || gp.nCc[1])
		{
			dcerror(2,"Should not use --bed-file %s if also using --case-file or --cont-file.\n", gp.bedFileFn);
			return 1;
		}
		strcpy(gp.ccFn[0][gp.nCc[0]++], gp.bedFileFn); // pretend bedFile is a control file for now
	}
	for (i=0;i<gp.nCc[0];++i)
		{
			gcont.setVariantFileName(gp.ccFn[0][i]);
			sprintf(fn,"gva.%s.cont.%d.vcf",geneName,i+1);
			if (gp.dontExtractVariants)
				printf("Will not attempt to produce %s because --dont-extract-variants was set\n",fn);
			else
				//if (!gcont.extractVariants(r,fn,0,spec.addChrInVCF[ff++],spec.removeVcfSpaces,spec.omitIntrons, spec.spliceRegionSize))
				//extractedOK=0;
			{
				if (gp.bedFileFn[0] == '\0') 
				{
//					if (!r.tbiExtractGene(gp.ccFn[0][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces, spec.omitIntrons, spec.spliceRegionSize))
					if (!tbiExtractIntervals(gp.ccFn[0][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces,iList))
						extractedOK = 0;
				}
				else 
				{
//					if (!r.plinkExtractGene(gp.bedFileFn,gp.famFileFn,gp.bimFileFn, fn, spec.omitIntrons, spec.spliceRegionSize))
					if (!plinkExtractIntervals(gp.bedFileFn, gp.famFileFn, gp.bimFileFn, fn, iList,geneName))
						extractedOK = 0;

				}
			}
			vf.addLocusFile(fn,VCFFILE);
			if (!vf.readLocusFileEntries(fn,spec,0))
				extractedOK=0;
		}
	for (i=0;i<gp.nCc[1];++i)
		{
			gcase.setVariantFileName(gp.ccFn[1][i]);
			sprintf(fn,"gva.%s.case.%d.vcf",geneName,i+1);
			if (gp.dontExtractVariants)
				printf("Will not attempt to produce %s because --dont-extract-variants was set\n",fn);
			else // if (!gcase.extractVariants(r,fn,0,spec.addChrInVCF[ff++], spec.removeVcfSpaces,spec.omitIntrons,spec.spliceRegionSize))
				// extractedOK=0;
//			if (!r.tbiExtractGene(gp.ccFn[1][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces, spec.omitIntrons, spec.spliceRegionSize))
			{
					if (!tbiExtractIntervals(gp.ccFn[1][i], fn, 0, spec.addChrInVCF[ff++], spec.removeVcfSpaces, iList))
						extractedOK = 0;
			}
			vf.addLocusFile(fn,VCFFILE);
			if (!vf.readLocusFileEntries(fn, spec, 1))
			{
				extractedOK = 0;
				printf("No valid variants were read from file %s\n", fn);
			}
		}
	if (gp.onlyExtractVariants)
	{
		printf("Exiting because --only-extract-variants was set\n");
		exit(0);
	}
	if (extractedOK)
		{
		if (spec.willNeedInbuiltConsequence)
		{
			if (gp.referencePath[0] == '\0')
			{
				dcerror(2, "Was going to annotate using inbuilt routines but cannot do this because --reference-folder has not been set\n");
				return 1;
			}
			if (spec.weightExpressions.size() == 0 && spec.recWeightExpressions.size() == 0 && !spec.willNeedEnsemblConsequence)
			{
				std::string * defaultStr = new std::string("ANNOT(\"INBUILT\")GETWEIGHT(\"DEFAULTWEIGHTS\")");
				// default behaviour is to use these weights unless told not to
				spec.weightExpressions.push_back(*defaultStr);
			}
			printf("Annotating using inbuilt routines...\n");
			vf.getQuickConsequences(r,spec);
		}
		if (spec.willNeedEnsemblConsequence)
		{
			if (spec.weightExpressions.size() == 0 && spec.recWeightExpressions.size() == 0)
			{
				std::string* defaultStr = new std::string("ANNOT(\"VEP\")GETWEIGHT(\"DEFAULTVEPWEIGHTS\")");
				spec.weightExpressions.push_back(*defaultStr);
			}
			printf("Annotating using VEP...\n");
			vf.getEnsemblConsequences(spec);
		}
		}
	sprintf(fn,"%s.%s",gp.testName,geneName);
	if (extractedOK)
	{
		geneVarParser::thisGene=&r; // essential for the annotation to work
		printf("Writing scoreassoc files...\n");
		if (vf.writeScoreAssocFiles(fn, spec.wf, gp.useFreqs, gp.nSubs, 1, gp.writeComments, gp.writeScoreFile, gp.writeRecScoreFile, spec))
		{
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
		else
		{
			if (!gp.doNotRun)
			{
				sprintf(fn2, "%s.sao", fn);
				fp = fopen(fn2, "w");
				fprintf(fp, "No valid variants for this gene\n");
				fclose(fp);
			}
			else
				printf("No valid variants for %s\n", fn);

		}
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
