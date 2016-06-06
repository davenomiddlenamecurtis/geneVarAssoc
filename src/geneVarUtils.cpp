#include "geneVarUtils.hpp"
#include <ctype.h>
#include <assert.h>

#define MAXSUB 10000

#define isArgType(a) (a[0]=='-' && a[1]=='-')
#define FILLARG(str) (strcmp(arg,str) ? 0 : ((getNextArg(arg, argc, argv, fp,&depth, &argNum) && !isArgType(arg)) ? 1 : (dcerror(1,"No value provided for argument: %s\n",str), 0)))

int gvaParams::readParms(int argc,char *argv[],analysisSpecs &spec)
{
	FILE *fp[MAXDEPTH],*ef;
	int argNum,i,s,usePhenotypes,f,depth;
	char arg[2000],line[2000],addChrStr[MAXVCFFILES+1],phenotypeFileName[200],samplesFileName[200];
	depth=-1;
	argNum=1;
	FILE *phenotypeFile;
	*phenotypeFileName=*samplesFileName='\0';
	geneListFn[0]=baitFn[0]=ccFn[2][MAXVCFPERCC][0]=referencePath[0]=geneName[0]=sequencePath[0]=posName[0]=intervalListFn[0]=testName[0]='\0';
	for (i=0;i<MAXVCFFILES;++i)
		spec.addChrInVCF[i]=0;
	if (spec.phenotypes!=NULL)
	{
		free(spec.phenotypes);
		spec.phenotypes=0;
	}
	baitFn[0]='\0';
	upstream=1000;
	downstream=0;
	margin=0;
	wf=10.0;
	wFunc=0;
	writeScoreFile=0;
	*addChrStr='\0';
	spec.useTrios=0;
	spec.useProbs=0;
	spec.useEnsembl=0;
	spec.consequenceThreshold=NULL_CONSEQUENCE;
	spec.useConsequenceWeights=0;
	spec.onlyUseSNPs=0;
	writeComments=0;
	writeScoreFile=0;
	spec.doRecessiveTest=0;
	spec.weightThreshold=0;
	spec.LDThreshold=1.0;
	spec.unknownIfUntyped=0;
	spec.skipIfNoPass=0;
	spec.unknownIfNoPass=1;
	spec.useHaplotypes=0;
	spec.GQThreshold=0;
	*referencePath=*sequencePath=*posName='\0';
	*spec.alleleFreqStr=*spec.alleleNumberStr=*spec.alleleCountStr='\0';
	spec.nExc=0;
	dontExtractGene=0;
	keepTempFiles=0;
	doNotRun=0;
	for (i = 0; i < 2; ++i)
	{
		useFreqs[i] = 0; // default
		nCc[i]=0;
		nSubs[i]=0;
	}
	while (getNextArg(arg, argc, argv, fp,&depth, &argNum))
	{
		if (!isArgType(arg))
		{
			dcerror(1,"Expected argument beginning -- but got this: %s\n",arg);
			return 0;
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
		else if (FILLARG("--phenotype-file"))
		{
			strcpy(phenotypeFileName,arg);
		}
		else if (FILLARG("--samples-file"))
		{
			strcpy(samplesFileName,arg);
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
		else if (FILLARG("--position"))
			strcpy(posName,arg);
		else if (FILLARG("--num-case"))
			nSubs[1]=atoi(arg);
		else if (FILLARG("--cont-file"))
			strcpy(ccFn[0][nCc[0]++], arg);
		else if (FILLARG("--case-file"))
			strcpy(ccFn[1][nCc[1]++], arg);
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
			wf=atof(arg);
		else if (FILLARG("--use-ensembl"))
			spec.useEnsembl=atoi(arg);
		else if (FILLARG("--consequence-threshold"))
			spec.consequenceThreshold=(consequenceType)atoi(arg);
		else if (FILLARG("--use-consequence-weights"))
			spec.useConsequenceWeights=atoi(arg);
		else if (FILLARG("--only-use-SNPs"))
			spec.onlyUseSNPs=atoi(arg);
		else if (FILLARG("--write-comments"))
			writeComments=atoi(arg);
		else if (FILLARG("--write-score-file"))
			writeScoreFile=atoi(arg);
		else if (FILLARG("--do-recessive-test"))
			spec.doRecessiveTest=atoi(arg);
		else if (FILLARG("--weight-threshold"))
			spec.weightThreshold=atof(arg);
		else if (FILLARG("--LD-threshold"))
			spec.LDThreshold=atof(arg);
		else if (FILLARG("--add-chr"))
			strcpy(addChrStr, arg);
		else if (FILLARG("--unknown-if-untyped"))
			spec.unknownIfUntyped=atoi(arg);
		else if (FILLARG("--unknown-if-no-pass"))
			spec.unknownIfNoPass=atoi(arg);
		else if (FILLARG("--skip-if-no-pass"))
			spec.skipIfNoPass=atoi(arg);
		else if (FILLARG("--ignore-alleles")) // treat loci with different allele sets as the same locus
			spec.ignoreAlleles=atoi(arg);
		else if (FILLARG("--use-haplotypes"))
			spec.useHaplotypes=atoi(arg);
		else if (FILLARG("--use-probs"))
			spec.useProbs=atoi(arg);
		else if (FILLARG("--dont-extract-gene"))
			dontExtractGene=atoi(arg);
		else if (FILLARG("--keep-temp-files"))
			keepTempFiles=atoi(arg);
		else if (FILLARG("--do-not-run"))
			doNotRun=atoi(arg);
		else if (FILLARG("--reference-path"))
			strcpy(referencePath, arg);
		else if (FILLARG("--sequence-path"))
			strcpy(sequencePath, arg);
		else if (FILLARG("--test-name"))
			strcpy(testName, arg);
		else if (FILLARG("--interval-list-file"))
			strcpy(intervalListFn, arg);
		else if (FILLARG("--GQ-threshold"))
			spec.GQThreshold=atof(arg);
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
				--depth; // because file will get cclosed
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
			for (;fgets(line,1999,ef) && sscanf(line,"%[^\n]",spec.exclusionStr[spec.nExc])==1;++spec.nExc)
				;
			fclose(ef);
		}
		else
		{
			dcerror(1,"Unrecognised argument: %s\n",arg);
			return 0;
		}
		;
	}
	if (phenotypeFileName[0] || samplesFileName[0])
	{
		spec.phenotypes = (int*)malloc(sizeof(int)*MAXSUB);
		if (phenotypeFileName[0])
		{
		assert((phenotypeFile = fopen(phenotypeFileName, "r"))!=0);
			if (phenotypeFile == NULL)
				dcerror(1, "Could not open phenotype file: %s\n", phenotypeFileName);
			nCc[0] = 0;
			for (s = 0; fgets(line, 1999, phenotypeFile) && sscanf(line, "%d", &spec.phenotypes[s]) == 1; ++s)
				;
		}
		else
		{
			int col,c;
			char *ptr;
			assert((phenotypeFile = fopen(samplesFileName, "r"))!=0);
			fgets(line, 1999, phenotypeFile);
			for (col = 0, ptr = line; toupper(ptr[0]) != 'P'&&toupper(ptr[0]) != 'H'; ++col)
			{
				while (!isspace(*ptr)) ++ptr;
				while (isspace(*ptr)) ++ptr;
			}
			fgets(line, 1999, phenotypeFile); //junk
			for (s=0;fgets(line, 1999, phenotypeFile) && !isspace(line[0]);++s) // breaks if somebody wants to start their sample file with spaces
			{
				for (c = 0,ptr=line; c < col; ++c)
				{
				while (!isspace(*ptr)) ++ptr;
				while (isspace(*ptr)) ++ptr;
				}
				spec.phenotypes[s]=*ptr-'0';
			}
		}
		fclose(phenotypeFile);
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
	return 1;
}

int gvaParams::getNextArg(char *nextArg, int argc,char *argv[], FILE *fp[MAXDEPTH],int *depth, int *argNum)
{
	*nextArg='\0';
	while (*depth>-1)
	{
		if (fscanf(fp[*depth],"%s ",nextArg)==1)
			return 1;
		else
		{
			fclose(fp[*depth]);
			--*depth;
		}
	}
	if (*argNum < argc)
	{
		strcpy(nextArg,argv[*argNum]);
		++ *argNum;
		return 1;
	}
	else
		return 0;
}

int gvaParams::input(FILE *fp,analysisSpecs &spec)
{
	char line[2000],rest[2000],phenotypeFn[1000];
	int i,usePhenotypes,s;
	FILE *phenotypeFile;
	for (i=0;i<MAXVCFFILES;++i)
		spec.addChrInVCF[i]=0;
	if (spec.phenotypes!=NULL)
	{
		free(spec.phenotypes);
		spec.phenotypes=0;
	}
	// ways to provide VCF files
	// filename -1 on first line
	// means filename has phenotype codes in, then there is just one file of cases and controls on next line
	// filename -2 on first line
	// means filename has trio pedigrees in, then there is just one file of cases and controls on next line
	// otherwise, control files are on first line, case files on second line
	// filename 1 nSubs means this file has allele frequencies only, derived from nSub subjects (can be for controls and/or cases)
	// filename0 filename1 filename2 ... control or case files, depending on whether on first or second line
	spec.useTrios=0;
	for (i=0;i<2;++i)
	{
		fgets(line,1999,fp);
		useFreqs[i]=0; // default
		if (sscanf(line, "%s %d", phenotypeFn, &usePhenotypes) == 2 && usePhenotypes == -2 && i == 0)
		{
			spec.useTrios=1;
			strcpy(spec.triosFn,phenotypeFn);
			nCc[0]=0;
		}
		else if (sscanf(line,"%s %d",phenotypeFn,&usePhenotypes)==2 && usePhenotypes==-1 && i==0)
		{
			spec.phenotypes=(int*)malloc(sizeof(int)*MAXSUB);
			phenotypeFile=fopen(phenotypeFn,"r");
			nCc[0]=0;
			for (s=0;fgets(line,1999,phenotypeFile) && sscanf(line,"%d",&spec.phenotypes[s])==1;++s)
				;
			fclose(phenotypeFile);
		}
		else if (sscanf(line,"%s %d %d",ccFn[i][0],&useFreqs[i],&nSubs[i])==3 && useFreqs[i]==1 && nSubs[i]>0)
			nCc[i]=1;
		else
		{
			for (nCc[i]=0,*rest='\0';sscanf(line,"%s %[^\n]",ccFn[i][nCc[i]],rest)>=1;++nCc[i])
			{
				strcpy(line,rest);
				*rest='\0';
			}
		}
	}
	fgets(line,1999,fp);
	baitFn[0]='\0';
	sscanf(line,"%s %s",geneListFn,baitFn);
	upstream=1000;
	downstream=0;
	margin=0;
	fgets(line,1999,fp);
	sscanf(line,"%d %d %d",&upstream,&downstream,&margin);
	fgets(line,1999,fp);
	wf=10.0;
	wFunc=0;
	writeScoreFile=0;
	sscanf(line,"%f %d",&wf,&wFunc);
	fgets(line,1999,fp);
	char addChrStr[MAXVCFFILES+1];
	*addChrStr='\0';
	sscanf(line,"%d %d %d %d %d %d %d %f %f %s %d %d %d %d",
		&spec.useEnsembl,
		&spec.consequenceThreshold,
		&spec.useConsequenceWeights,
		&spec.onlyUseSNPs,
		&writeComments,
		&writeScoreFile,
		&spec.doRecessiveTest,
		&spec.weightThreshold,
		&spec.LDThreshold,
		addChrStr,
		&spec.unknownIfUntyped,
		&spec.skipIfNoPass,
		&spec.unknownIfNoPass,
		&spec.useHaplotypes);
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
	*line='\0';
	fgets(line,1999,fp);
	sscanf(line,"%s",referencePath);
	spec.GQThreshold=0;
	fgets(line,1999,fp);
	sscanf(line,"%f",&spec.GQThreshold);
	for (spec.nExc=0;fgets(line,1999,fp) && sscanf(line,"%[^\n]",spec.exclusionStr[spec.nExc])==1;++spec.nExc)
		;
	return 1;
}
