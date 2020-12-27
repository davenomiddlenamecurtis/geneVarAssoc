#include "masterLocusFile.hpp"
#include "geneVarUtils.hpp"
#ifndef MSDOS 
#include <unistd.h>
#endif
#include <ctype.h>
#include <string.h>

int masterLocusFile::writeAltSubs(char* fn, analysisSpecs& spec, char* posName, char* refAll, char* altAll)
{
	int i, totalSub, altAllNum, a, refAllNum;
	FILE* fo;
	totalSub = 0;
	for (i = 0, totalSub = 0; i < nLocusFiles; ++i)
		totalSub += nSubs[i];
	allelePair* all;
	all = (allelePair*)calloc(totalSub, sizeof(allelePair));
	strEntry* subName;
	subName = (strEntry*)calloc(totalSub, sizeof(strEntry));
	openLocusFiles();
	outputSubNames(subName, spec);
	altAllNum = -1;
	if (gotoFirstInRange(spec))
	{
		if (altAll[0] == '\0')
			outputCurrentAlleles(all, spec); // just use the first variant at this position
		else {
			do {
				if (refAll[0] && strcmp(refAll, tempRecord.alls[0]))
					continue;
				for (a = 0; a < tempRecord.nAlls; ++a)
				{
					if (!strcmp(altAll, tempRecord.alls[a]))
					{
						altAllNum = a;
						outputCurrentAlleles(all, spec);
						break;
					}
				}
				if (altAllNum != -1)
					break;
			} while (gotoNextInRange(spec));
			if (altAllNum == -1)
				dcerror(1, "Could not find variant with allele %s at position %s\n", altAll, posName);
		}
	}
	else
	{
		dcerror(1, "No valid variants found at position %s\n", posName);
		free(all);
		free(subName);
		return 0;
	}
	if (fn == 0)
		fo = stdout;
	else
		fo = fopen(fn, "w");
	if (spec.outputRef)
	{
		for (i = 0; i < totalSub; ++i)
			if (all[i][0] != 0 && (all[i][0] == 1 || all[i][1] == 1))
				fprintf(fo, "%s\t%d %d\n", subName[i], all[i][0], all[i][1]);
	}
	else
	{
		for (i = 0; i < totalSub; ++i)
			if (all[i][0] != 0 && (altAllNum == -1 && (all[i][0] != 1 || all[i][1] != 1) || (all[i][0] == altAllNum + 1 || all[i][1] == altAllNum + 1))) // all[i][0] is 1 for REF allele
				fprintf(fo, "%s\t%d %d\n", subName[i], all[i][0], all[i][1]);
	}
	if (fn!=0)
		fclose(fo);
	free(all);
	free (subName);
	return 1;
}

int main(int argc,char *argv[])
{
	char fn[100],fn2[100],line[1000],vcfFn[100],vcfFnBuff[100],chrStr[20],*ptr;
	int i,totalSub,cc,p, extractedOK;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;
	intervalList iList;
	strcpy(gp.testName, "altSubs");
	if (!gp.readParms(argc,argv,spec))
		exit(1);
	spec.useConsequenceWeights=0; // I am not going to annotate these variants
	masterLocusFile vf(gp.bedFileFn[0] ? 1 : gp.nCc[0]+gp.nCc[1]);
	if (sscanf(gp.posName, "%[^:]:%d", chrStr, &p)!=2)
		dcerror(1,"Usage: showAltSubs --arg-file something.arg --position 7:12139555 [--allele G]");
	iList.append(chrStr, p, p);

	spec.sc = spec.ec = gp.posName[0] == 'X' ? 23 : atoi(gp.posName);
	spec.sp = spec.ep = p;
	sprintf(fn,"gva.%s.db",gp.testName);
	sprintf(fn2,"gva.%s.vdx",gp.testName);
	unlink(fn);
	unlink(fn2);
	vf.openFiles(fn,fn2);
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
	int ff = 0; 
	for (cc=0;cc<2;++cc)
		for (i=0;i<gp.nCc[cc];++i)
		{
			sprintf(fn,"gva.%s.%s.%d.vcf",cc?"case":"cont",gp.testName,i+1);
			if (gp.dontExtractVariants)
				printf("Will not attempt to produce %s because --dont-extract-variants was set\n", fn);
			else
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
			vf.addLocusFile(fn,VCFFILE);
			if (!vf.readLocusFileEntries(fn, spec, 1)) // this is just an easy way to make sure that files are known about
				extractedOK = 0;
		}
	// sprintf(fn,"%s_%s.aso",argv[2],ptr+1);
//	if (gp.testName[0])
		sprintf(fn, "%s.txt", gp.testName);
//	else
//		sprintf(fn, "altSubs.%s.txt", gp.posName);
	if (extractedOK)
		vf.writeAltSubs(fn,spec, gp.posName,gp.refAll, gp.altAll);
	return 0;
}
