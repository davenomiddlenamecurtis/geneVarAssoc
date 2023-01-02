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

#include "geneVarUtils.hpp"
#include <ctype.h>
#include <assert.h>

#ifndef MAXSUB
#define MAXSUB 20000
#endif

#define isArgType(a) (a[0]=='-' && a[1]=='-')
#define FILLARG(str) (strcmp(arg,str) ? 0 : ((getNextArg(arg, argc, argv, fp,&depth, &argNum) && !isArgType(arg)) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))

int gvaParams::readParms(int argc,char *argv[],analysisSpecs &spec)
{
	FILE *fp[MAXDEPTH],*ef;
	int argNum,i,s,usePhenotypes,f,depth,nPhenotypeFile,nIDsAndPhenotypeFile,nSamplesFile;
	char arg[2000],line[2000],addChrStr[MAXVCFFILES+1],phenotypeFileName[MAXVCFFILES][100];
	depth=-1;
	argNum=1;
	FILE *phenotypeFile,*multiWeightFile,* locusPosFile;
	nPhenotypeFile=nIDsAndPhenotypeFile=nSamplesFile=0;
	spec.nScoreassocArgs = 0;
	geneListFn[0]=baitFn[0]=referencePath[0]=geneName[0]=sequencePath[0]=posName[0]= refAll[0] = altAll[0] = intervalListFn[0]='\0';
	bedFileFn[0] = bimFileFn[0] = famFileFn[0] = '\0';
	// testName can be set by calling function, e.g. default "gva" for geneVarAssoc
	spec.outputRef =  0;
	strcpy(spec.vepCommand,"perl variant_effect_predictor.pl");
//	strcpy(spec.weightExpression,"ANNOT(\"DEFAULT\")GETWEIGHT(\"DEFAULTWEIGHTS\")");
	spec.excludeExpressions.clear();
	spec.weightExpressions.clear();
	spec.recWeightExpressions.clear();
	spec.weightNames.clear();
	spec.recWeightNames.clear();
	for (i=0;i<MAXVCFFILES;++i)
		spec.addChrInVCF[i]=0;
	if (spec.phenotypes!=NULL)
	{
		free(spec.phenotypes);
		spec.phenotypes=0;
	}
	spec.subPhenos.clear();
	baitFn[0]='\0';
	upstream=1000;
	downstream=0;
	margin=0;
	spec.wf=10.0;
	writeScoreFile=0;
	*addChrStr='\0';
	spec.isQuantitative = 0;
	spec.useTrios=0;
	spec.useProbs=0;
	spec.useEnsembl=spec.willNeedEnsemblConsequence=spec.willNeedInbuiltConsequence=0;
	spec.consequenceThreshold=NULL_CONSEQUENCE;
	spec.useConsequenceWeights=0;
	spec.onlyUseSNPs=0;
	writeComments=1;
writeScoreFile = 0;
spec.doRecessiveTest = 0;
spec.recWeightThreshold = 0;
spec.LDThreshold = 1.0;
spec.unknownIfUntyped = 0;
spec.skipIfNoPass = 0;
spec.unknownIfNoPass = 0;
spec.showHapLocusNames = spec.useHaplotypes = 0;
spec.GQThreshold = 0;
spec.depthThreshold = spec.hetDevThreshold = spec.hetDevThresholdSq = spec.ABThreshold = -1;
spec.ignoreAlleles = 0;
spec.dontMergeAlleles = 0;
*referencePath = *sequencePath = *posName = *altAll = '\0';
*spec.alleleFreqStr = *spec.alleleNumberStr = *spec.alleleCountStr = '\0';
spec.nExc = 0;
dontExtractVariants = onlyExtractVariants = 0;
keepTempFiles = 0;
doNotRun = 0;
spec.debug = 0;
spec.numVcfFieldsToSkip = DEFAULTNUMVCFFIELDSTOSKIP;
spec.removeVcfSpaces = 0;
for (i = 0; i < 2; ++i)
{
	useFreqs[i] = 0; // default
	nCc[i] = 0;
	nSubs[i] = 0;
}
while (getNextArg(arg, argc, argv, fp, &depth, &argNum))
{
	if (!isArgType(arg))
	{
		dcerror(1, "Expected argument beginning -- but got this: %s\n", arg);
		return 0;
	}
	else if (!strcmp(arg, "--dottest") || !strcmp(arg, "--dolrtest")
		|| !strcmp(arg, "--dolinrtest") || !strcmp(arg, "--varfile")
		|| !strcmp(arg, "--testfile") || !strcmp(arg, "--lintestfile")
		|| !strcmp(arg, "--start-from-fitted") || !strcmp(arg, "--maxmaf")
		|| !strcmp(arg, "--maxrecloci") || !strcmp(arg, "--ldthreshold")
		|| !strcmp(arg, "--lamda") || !strcmp(arg, "--missingzero"))
	{
		strcpy(spec.scoreassocArgs[spec.nScoreassocArgs][0], arg);
		getNextArg(arg, argc, argv, fp, &depth, &argNum);
		strcpy(spec.scoreassocArgs[spec.nScoreassocArgs][1], arg);
		++spec.nScoreassocArgs;
	}
	else if (FILLARG("--arg-file"))
	{
		if (++depth >= MAXDEPTH)
		{
			dcerror(1, "Attempting to recurse too deeply into arg-files with this one: %s\n", arg);
			return 0;
		}
		else
		{
			fp[depth] = fopen(arg, "r");
			if (fp[depth] == NULL)
			{
				dcerror(1, "Could not open arg file: %s\n", arg);
				return 0;
			}
		}
	}
	else if (FILLARG("--comment-expression"))
	{
		strcpy(spec.commentExpression, arg);
	}
	else if (FILLARG("--weight-name"))
	{
		std::string* wnStr = new std::string(arg);
		spec.weightNames.push_back(*wnStr);
	}
	else if (FILLARG("--rec-weight-name"))
	{
		std::string* wnStr = new std::string(arg);
		spec.recWeightNames.push_back(*wnStr);
	}
	else if (FILLARG("--weight-name-file"))
	{
		FILE* wnf;
		char weightName[200];
		wnf = fopen(arg, "r");
		if (wnf == NULL)
			dcerror(1, "Could not open weight name file %s\n", arg);
		while (fscanf(wnf, "%s", weightName) == 1)
		{
			std::string* wnStr = new std::string(weightName);
			spec.weightNames.push_back(*wnStr);
		}
		fclose(wnf);
	}
	else if (FILLARG("--rec-weight-name-file"))
	{
		FILE* wnf;
		char weightName[200];
		wnf = fopen(arg, "r");
		if (wnf == NULL)
			dcerror(1, "Could not open weight name file %s\n", arg);
		while (fscanf(wnf, "%s", weightName) == 1)
		{
			std::string* wnStr = new std::string(weightName);
			spec.recWeightNames.push_back(*wnStr);
		}
		fclose(wnf);
	}
	else if (FILLARG("--weight-expression"))
	{
		std::string* argStr = new std::string(arg);
		spec.weightExpressions.push_back(*argStr);
		if (strstr(arg, "INBUILT"))
			spec.willNeedInbuiltConsequence = 1;
	}
	else if (FILLARG("--rec-weight-expression"))
	{
		std::string* argStr = new std::string(arg);
		spec.recWeightExpressions.push_back(*argStr);
		if (strstr(arg, "INBUILT"))
			spec.willNeedInbuiltConsequence = 1;
	}
	else if (FILLARG("--exclude-expression"))
	{
		std::string* argStr = new std::string(arg);
		spec.excludeExpressions.push_back(*argStr);
		//			if(strstr(arg,"VEP"))
		//				spec.willNeedEnsemblConsequence=1;
		if (strstr(arg, "INBUILT"))
			spec.willNeedInbuiltConsequence = 1;
	}
	else if (FILLARG("--multi-weight-file"))
	{
	char multiWeightVariantFile[MAXFILENAMELENGTH], varWeightName[MAXFILENAMELENGTH];
	multiWeightFile = fopen(arg, "r");
	if (multiWeightFile == NULL)
		dcerror(1, "Could not open multi-weight-file: %s\n", arg);
	if (!fgets(line, MAXFILENAMELENGTH - 1, multiWeightFile) || sscanf(line, "%s", multiWeightVariantFile) != 1)
		dcerror(1, "Could not read name of variant file for multiple weights from multi-weight-file: %s\n", arg);
	while (fgets(line, MAXFILENAMELENGTH - 1, multiWeightFile) && sscanf(line, "%s", varWeightName) == 1)
	{
		sprintf(line, "%c%s%cDBNSFPLOOKUP%c%s%c", '"', varWeightName, '"', '"', multiWeightVariantFile, '"');
		spec.weightExpressions.push_back(*(new std::string(line)));
		spec.weightNames.push_back(*(new std::string(varWeightName)));
	}
	fclose(multiWeightFile);
	}
	else if (FILLARG("--count-these-loci"))
	{
	char locusPos[MAXFILENAMELENGTH], alls[MAXFILENAMELENGTH];
	locusPosFile = fopen(arg, "r");
	if (locusPosFile == NULL)
		dcerror(1, "Could not open count-these-loci file: %s\n", arg);
	while (alls[0]='\0',fgets(line, MAXFILENAMELENGTH, locusPosFile) && sscanf(line, "%s %s", locusPos,alls) >= 1)
	{
		if (alls[0] == '\0')
		{
			sprintf(line, "(ATTRIB(%cPOS%c)=%c%s%c)*1", '"', '"', '"', locusPos, '"');
		}
		else 
		{
			sprintf(line, "(ATTRIB(%cPOS%c)=%c%s%c && ATTRIB(%cALLS%c)=%c%s%c)*1", '"', '"', '"', locusPos, '"', '"', '"', '"', alls, '"');
		}
		spec.weightExpressions.push_back(*(new std::string(line)));
		if (alls[0] == '\0')
		{
			sprintf(line, "POS%s", locusPos);
		}
		else
		{
			sprintf(line, "POS%s-%s", locusPos,alls);
		}
		spec.weightNames.push_back(*(new std::string(line)));
	}
	fclose(locusPosFile);
	}
	else if (FILLARG("--phenotype-file"))
	{
			if (nIDsAndPhenotypeFile || nSamplesFile)
				dcerror(1, "Can only have one of --phenotype-file, --ID-and-phenotype-file and --samples-file");
			strcpy(phenotypeFileName[nPhenotypeFile++], arg);
	}
	else if (FILLARG("--ID-and-phenotype-file"))
		{
			if (nPhenotypeFile||nSamplesFile)
				dcerror(1,"Can only have one of --phenotype-file, --ID-and-phenotype-file and --samples-file");
			strcpy(phenotypeFileName[nIDsAndPhenotypeFile++],arg);
		}
	else if (FILLARG("--samples-file"))
		{
			if (nIDsAndPhenotypeFile||nPhenotypeFile)
				dcerror(1,"Can only have one of --phenotype-file, --ID-and-phenotype-file and --samples-file");
			strcpy(phenotypeFileName[nSamplesFile++],arg);
		}
	else if (FILLARG("--trio-file"))
		{
			spec.useTrios=1;
			strcpy(spec.triosFn,arg);
			nCc[0]=0;
		}
	else if (FILLARG("--cont-freq-file"))
		{
			strcpy(ccFn[0][0],arg);
			useFreqs[0]=1;
			nCc[0]=1;
		}
	else if (FILLARG("--num-cont"))
			nSubs[0]=atoi(arg);
		else if (FILLARG("--case-freq-file"))
		{
			strcpy(ccFn[1][0],arg);
			useFreqs[1]=1;
			nCc[1]=1;
		}
		else if (FILLARG("--num-case"))
			nSubs[1] = atoi(arg);
		else if (FILLARG("--position"))
			strcpy(posName, arg);
		else if (FILLARG("--alt-all"))
			strcpy(altAll, arg);
		else if (FILLARG("--ref-all"))
			strcpy(refAll, arg);
		else if (FILLARG("--output-ref"))
			spec.outputRef = atoi(arg);
		else if (FILLARG("--cont-file"))
			strcpy(ccFn[0][nCc[0]++], arg);
		else if (FILLARG("--case-file"))
			strcpy(ccFn[1][nCc[1]++], arg);
		else if (FILLARG("--bed-file"))
			strcpy(bedFileFn, arg);
		else if (FILLARG("--fam-file"))
			strcpy(famFileFn, arg);
		else if (FILLARG("--bim-file"))
		strcpy(bimFileFn, arg);
		else if (FILLARG("--ref-gene-file"))
			strcpy(geneListFn, arg);
		else if (FILLARG("--allele-freq-str"))
			strcpy(spec.alleleFreqStr, arg);
		else if (FILLARG("--allele-number-str"))
			strcpy(spec.alleleNumberStr, arg);
		else if (FILLARG("--allele-count-str"))
			strcpy(spec.alleleCountStr, arg);
		else if (FILLARG("--bait-file"))
			strcpy(baitFn, arg);
		else if (FILLARG("--upstream"))
			upstream=atoi(arg);
		else if (FILLARG("--downstream"))
			downstream=atoi(arg);
		else if (FILLARG("--margin"))
			margin=atoi(arg);
		else if (FILLARG("--weight-factor"))
			spec.wf=atof(arg);
		else if (FILLARG("--isquantitative"))
			spec.isQuantitative = atoi(arg);
		else if (FILLARG("--use-flat-file"))
			spec.useFlatFile = atoi(arg);
		else if (FILLARG("--use-transposed-file"))
			spec.useTransposedFile = atoi(arg);
		else if (FILLARG("--multiline-vep"))
			spec.multilineVEP = atoi(arg);
		else if (FILLARG("--use-ensembl"))
			spec.useEnsembl = atoi(arg);
		else if (FILLARG("--consequence-threshold"))
			spec.consequenceThreshold = atof(arg);
		else if (FILLARG("--use-consequence-weights"))
			spec.useConsequenceWeights=atoi(arg);
		else if (FILLARG("--only-use-SNPs"))
			spec.onlyUseSNPs=atoi(arg);
		else if (FILLARG("--write-comments"))
			writeComments=atoi(arg);
		else if (FILLARG("--write-score-file"))
			writeScoreFile = atoi(arg);
		else if (FILLARG("--write-rec-score-file"))
			writeRecScoreFile = atoi(arg);
#if 0
		else if (FILLARG("--do-recessive-test"))
			spec.doRecessiveTest=atoi(arg);
		else if (FILLARG("--rec-weight-threshold"))
			spec.recWeightThreshold=atof(arg);
		else if (FILLARG("--LD-threshold"))
			spec.LDThreshold=atof(arg);
		else if (FILLARG("--use-haplotypes"))
			spec.useHaplotypes = atoi(arg);
		else if (FILLARG("--show-hap-locus-names"))
			spec.showHapLocusNames = atoi(arg);
#endif
		else if (FILLARG("--add-chr"))
			strcpy(addChrStr, arg);
		else if (FILLARG("--unknown-if-untyped"))
			spec.unknownIfUntyped = atoi(arg);
		else if (FILLARG("--wild-if-unknown"))
			spec.wildIfUnknown = atoi(arg);
		else if (FILLARG("--unknown-if-no-pass"))
			spec.unknownIfNoPass=atoi(arg);
		else if (FILLARG("--skip-if-no-pass"))
			spec.skipIfNoPass=atoi(arg);
		else if (FILLARG("--ignore-alleles")) // treat loci with different allele sets as the same locus
			spec.ignoreAlleles = atoi(arg);
		else if (FILLARG("--dont-merge-alleles")) // do not merge allelic systems from different vcf lines
			spec.dontMergeAlleles = atoi(arg);
		else if (FILLARG("--use-probs"))
			spec.useProbs=atoi(arg);
		else if (FILLARG("--dont-extract-variants"))
			dontExtractVariants = atoi(arg);
		else if (FILLARG("--only-extract-variants"))
			onlyExtractVariants = atoi(arg);
		else if (FILLARG("--keep-temp-files"))
			keepTempFiles=atoi(arg);
		else if(FILLARG("--do-not-run"))
			doNotRun=atoi(arg);
		else if (FILLARG("--debug"))
			spec.debug = atoi(arg);
		else if (FILLARG("--num-fields-to-skip"))
			spec.numVcfFieldsToSkip = atoi(arg);
		else if (FILLARG("--remove-vcf-spaces"))
			spec.removeVcfSpaces = atoi(arg);
		else if (FILLARG("--merge-alt-alleles"))
			spec.mergeAltAlleles = atoi(arg);
		else if (FILLARG("--omit-introns"))
			spec.omitIntrons = atoi(arg);
		else if (FILLARG("--splice-region-size"))
			spec.spliceRegionSize = atoi(arg);
		else if(FILLARG("--reference-path"))
			strcpy(referencePath,arg);
		else if(FILLARG("--vep"))
			strcpy(spec.vepCommand,arg);
		else if (FILLARG("--sequence-path"))
			strcpy(sequencePath, arg);
		else if (FILLARG("--test-name"))
			strcpy(testName, arg);
		else if (FILLARG("--interval-list-file"))
			strcpy(intervalListFn, arg);
		else if (FILLARG("--GQ-threshold"))
			spec.GQThreshold = atof(arg);
		else if (FILLARG("--depth-threshold"))
			spec.depthThreshold = atof(arg);
		else if (FILLARG("--AB-threshold"))
			spec.ABThreshold = atof(arg);
		else if(FILLARG("--hetdev-threshold"))
		{
			spec.hetDevThreshold=atof(arg);
			spec.hetDevThresholdSq=spec.hetDevThreshold*spec.hetDevThreshold;
		}
		else if (FILLARG("--gene"))
			strcpy(geneName, arg);
		else if (FILLARG("--clear-cont"))
		{
			if (atoi(arg) == 1)
				nCc[0] = 0;
		}
		else if (FILLARG("--clear-case"))
		{
			if (atoi(arg) == 1)
				nCc[1] = 0;
		}
		else if (FILLARG("--exclusion-list"))
		{
			if (!strcmp(arg, "-"))
			{
				if (depth == -1)
				{
					dcerror(1,"Cannot specify exclusion list on command line with no file as argument.\n");
					return 0;
				}
				ef=fp[depth];
				--depth; // because file will get closed
			}
			else
			{
				ef=fopen(arg,"r");
				if (ef == NULL)
				{
					dcerror(1,"Could not open exclusion list file: %s\n",arg);
					return 0;
				}
			}
			while (fgets(line, 1999, ef))
				if (sscanf(line, "%[^\n]", spec.exclusionStr[spec.nExc]) == 1)
					++spec.nExc;
			fclose(ef);
		}
		else
		{
			dcerror(1,"Unrecognised argument: %s\n",arg);
			return 0;
		}
		;
	}
	if (bedFileFn[0] && (!famFileFn[0] || !bimFileFn[0]))
	{
		dcerror(5, "If use --bed-file must also use --fam-file and --bim-file\n");
		return 0;
	}
	if (spec.isQuantitative && !nIDsAndPhenotypeFile)
	{
		dcerror(4, "Can only use --is-quantitative 1 if also use --ID-and-phenotype-file\n");
		return 0;
	}
	if (nPhenotypeFile || nIDsAndPhenotypeFile || nSamplesFile)
	{
		spec.phenotypes = (float*)malloc(sizeof(float)*MAXSUB);
		s=0;
		if (nPhenotypeFile)
			for (i=0;i<nPhenotypeFile;++i)
		{
		phenotypeFile = fopen(phenotypeFileName[i], "r");
		if (phenotypeFile == NULL)
				dcerror(1, "Could not open phenotype file: %s\n",phenotypeFileName[i]);
			nCc[0] = 0;
			for (; fgets(line, 1999, phenotypeFile) && sscanf(line, "%f", &spec.phenotypes[s]) == 1; ++s)
				;
		fclose(phenotypeFile);
		}
		else if (nIDsAndPhenotypeFile)
			for (i=0;i<nIDsAndPhenotypeFile;++i)
		{
			char ID[100];
			float phen;
			phenotypeFile = fopen(phenotypeFileName[i], "r");
			if (phenotypeFile == NULL)
				dcerror(1,"Could not open ID and phenotype file: %s\n",phenotypeFileName[i]);
			nCc[0] = 0;
			for (; phen=MISSINGPHENOTYPE,fgets(line,1999,phenotypeFile) && sscanf(line,"%s %f",ID,&phen) >=1; ++s) // allow that phen may be NA
				spec.subPhenos.insert(TStrFloatPair(ID,phen));
			fclose(phenotypeFile);
			}
		else
			for (i=0;i<nSamplesFile;++i)
			{
			int col,c;
			char *ptr;
			assert((phenotypeFile = fopen(phenotypeFileName[i], "r"))!=0);
			fgets(line, 1999, phenotypeFile);
			for (col = 0, ptr = line; toupper(ptr[0]) != 'P'&&toupper(ptr[0]) != 'H'; ++col)
			{
				while (!isspace(*ptr)) ++ptr;
				while (isspace(*ptr)) ++ptr;
			}
			fgets(line, 1999, phenotypeFile); //junk
			for (;fgets(line, 1999, phenotypeFile) && !isspace(line[0]);++s) // breaks if somebody wants to start their sample file with spaces
			{
				for (c = 0,ptr=line; c < col; ++c)
				{
				while (!isspace(*ptr)) ++ptr;
				while (isspace(*ptr)) ++ptr;
				}
				spec.phenotypes[s]=*ptr-'0';
			}
			fclose(phenotypeFile);
			} // I have made phenotypes a float rather than int but hopefully will all still work
	}
	int len;
	if (*addChrStr != '\0')
	{
		if ((len=strlen(addChrStr))==1)
			for (i=0;i<MAXVCFFILES;++i)
				spec.addChrInVCF[i]=addChrStr[0]-'0';
		else
			for (i=0;i<len;++i)
				spec.addChrInVCF[i]=addChrStr[i]-'0';
	}
	if (*sequencePath)
		for (i=0;i<2;++i)
			for (f = 0; f < nCc[i]; ++f)
			{
#ifdef MSDOS
				sprintf(line,"%s\\%s",sequencePath,ccFn[i][f]);
#else
				sprintf(line,"%s/%s",sequencePath,ccFn[i][f]);
#endif
				strcpy(ccFn[i][f],line);
			}
	if ((spec.consequenceThreshold || spec.useConsequenceWeights) && spec.weightExpressions.size() == 0 && spec.recWeightExpressions.size() == 0)
	{
		if (spec.useEnsembl)
			spec.willNeedEnsemblConsequence=1;
		else
			spec.willNeedInbuiltConsequence=1;
	}
	if (spec.useEnsembl)
		spec.willNeedEnsemblConsequence = 1; // I cannot think why useEnsembl would be set unless needed to get these
	return 1;
}

int gvaParams::getNextArg(char *nextArg, int argc,char *argv[], FILE *fp[MAXDEPTH],int *depth, int *argNum)
{
	char *ptr;
	int ch;
	*nextArg='\0';
	while (*depth>-1)
	{
		do {
			ch=fgetc(fp[*depth]);
		} while (ch!=EOF && isspace(ch));
		if (ch=='\'')
		{
			ptr=nextArg;
			while ((ch=fgetc(fp[*depth]))!='\'' && ch!=EOF)
				*ptr++=ch;
			*ptr='\0';
			return 1;
		}
		else if(ch!=EOF)
		{
			nextArg[0]=ch;
			ptr=nextArg+1;
			while((ch=fgetc(fp[*depth]))!=EOF && !isspace(ch))
				*ptr++=ch;
			*ptr='\0';
			return 1;
		}
		else
		{
			fclose(fp[*depth]);
			--*depth;
		}
	}
	if (*argNum < argc)
	{
		if (argv[*argNum][0]=='\'' && (ptr=strchr(argv[*argNum]+1,'\''))!=0)
		{
			*ptr='\0';
			strcpy(nextArg,argv[*argNum]+1);
		}
		else
			strcpy(nextArg,argv[*argNum]);
		++ *argNum;
		return 1;
	}
	else
		return 0;
}

int plinkExtractIntervals(char* bedFilename, char* famFilename, char* bimFilename, char* outFn,intervalList& iList, char *geneName)
{
	int i, startPos, endPos, foundOne, systemStatus;
	char buff[1000], * bedFn, * ptr, bedFnBuff[1000], * bimFn, bimFnBuff[1000],rfFnBuff[1000];
	long lineStart;
	FILE* rf;
	if (geneName == 0)
		geneName = "NOGENE";
	if ((ptr = strchr(bedFilename, '*')) == 0)
		bedFn = bedFilename;
	else
	{
		strcpy(bedFnBuff, bedFilename);
		ptr = strchr(bedFnBuff, '*');
		*ptr = '\0';
		strcat(bedFnBuff, iList.ints[0].chr);
		ptr = strchr(bedFilename, '*');
		strcat(bedFnBuff, ptr + 1);
		bedFn = bedFnBuff;
	}
	if ((ptr = strchr(bimFilename, '*')) == 0)
		bimFn = bimFilename;
	else
	{
		strcpy(bimFnBuff, bimFilename);
		ptr = strchr(bimFnBuff, '*');
		*ptr = '\0';
		strcat(bimFnBuff, iList.ints[0].chr);
		ptr = strchr(bimFilename, '*');
		strcat(bimFnBuff, ptr + 1);
		bimFn = bimFnBuff;
	}
	sprintf(rfFnBuff, "range.temp.%s.txt", geneName); // allow analyses to run simultaneously
	rf = fopen(rfFnBuff, "w");
	if (rf == 0)
	{
		dcerror(5, "Could not open file: %s for writing\n", rfFnBuff);
		return 0;
	}
	for (i=0;i<iList.nInts;++i)
		fprintf(rf, "%s %d %d %s\n", iList.ints[i].chr, iList.ints[i].st, iList.ints[i].en,geneName);
	fclose(rf);
	strcpy(buff, outFn);
	if ((ptr = strstr(buff, ".vcf")) != 0)
		*ptr = '\0'; // because plink appends .vcf to outfile name
	sprintf(refseqGeneInfo::geneLine, "plink --bed %s --fam %s --bim %s --extract range %s --set-hh-missing --recode vcf-iid --out %s",
		bedFn, famFilename, bimFn, rfFnBuff,buff);
	// added --set-hh-missing because was getting Warning: 40548 het. haploid genotypes present
	printf("Running command: %s\n", refseqGeneInfo::geneLine);
	systemStatus = system(refseqGeneInfo::geneLine);
	return 1;
}

int tbiExtractIntervals(char* tbiFilename, char* outFn, int appendToOld, int addChrInVCF, int removeSpaces, intervalList& iList)
{
	int i, systemStatus;
	char buff[1000], * tbiFn, * ptr, tbiFnBuff[1000];
	long lineStart;
	if ((ptr = strchr(tbiFilename, '*')) == 0)
		tbiFn = tbiFilename;
	else
	{
		strcpy(tbiFnBuff, tbiFilename);
		ptr = strchr(tbiFnBuff, '*');
		*ptr = '\0';
		strcat(tbiFnBuff, iList.ints[0].chr);
		ptr = strchr(tbiFilename, '*');
		strcat(tbiFnBuff, ptr + 1);
		tbiFn = tbiFnBuff;
	}

	sprintf(refseqGeneInfo::geneLine, "tabix %s%s ", tbiFn, appendToOld ? "" : " -h");
	for (i = 0; i < iList.nInts; ++i)
	{
		sprintf(strchr(refseqGeneInfo::geneLine, '\0'), "%s%s:%d-%d ",
			addChrInVCF ? "chr" : "", iList.ints[i].chr, iList.ints[i].st, iList.ints[i].en);
		if (strlen(refseqGeneInfo::geneLine) > 3900) // line getting long, extract baits so far then go back for more - maximum is 4000?
			{
				sprintf(strchr(refseqGeneInfo::geneLine, '\0'), "%s %s %s", removeSpaces ? "| sed s/' '/'_'/g " : "", appendToOld ? ">>" : ">", outFn);
				printf("Running command: %s\n", refseqGeneInfo::geneLine);
				checkSystem();
				systemStatus = system(refseqGeneInfo::geneLine);
				// printf("system returned %d\n",systemStatus);
				appendToOld = 1;
				sprintf(refseqGeneInfo::geneLine, "tabix %s ", tbiFn);
			}
		}
	sprintf(buff, "tabix %s ", tbiFn);
	if (strlen(refseqGeneInfo::geneLine) > strlen(buff))
	{
		sprintf(strchr(refseqGeneInfo::geneLine, '\0'), "%s %s %s", removeSpaces ? "| sed s/' '/'_'/g " : "", appendToOld ? ">>" : ">", outFn);
		printf("Running command: %s\n", refseqGeneInfo::geneLine);
		checkSystem();
		systemStatus = system(refseqGeneInfo::geneLine);
	}
		// printf("system returned %d\n",systemStatus);
	return 1;
}


