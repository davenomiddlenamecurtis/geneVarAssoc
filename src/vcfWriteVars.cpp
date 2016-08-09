#include "vcfLocusFile.hpp"
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#define MAXLOCIINGENE 2000



int masterLocusFile::writeVars(char *fn,int *useFreqs,analysisSpecs &spec)
{
	char buff[1000];
	int i,l,ll,totalSub,c,lc;
	FILEPOSITION recPos;
	const char *testKey;
	FILE *fp;
	strEntry *locusName,*consequence,*alleles,*locusInfo,*subName,**call;
	assert ((locusName=(strEntry *)calloc(MAXLOCIINGENE,sizeof(strEntry)))!=0);
	assert ((consequence=(strEntry *)calloc(MAXLOCIINGENE,sizeof(strEntry)))!=0);
	assert ((alleles=(strEntry *)calloc(MAXLOCIINGENE,sizeof(strEntry)))!=0);
	assert ((locusInfo=(strEntry *)calloc(MAXLOCIINGENE,sizeof(strEntry)))!=0);
	openLocusFiles();
	for (i=0,totalSub=0;i<nLocusFiles;++i)
		totalSub+=nSubs[i];
	assert ((subName=(strEntry *)calloc(totalSub,sizeof(strEntry)))!=0);
	assert ((call=(strEntry **)calloc(totalSub,sizeof(strEntry *)))!=0);
	for (i=0;i<totalSub;++i)
		assert ((call[i]=(strEntry *)calloc(MAXLOCIINGENE,sizeof(strEntry)))!=0);

	outputSubNames(subName,spec);
	recPos=findFirstInRange(spec);
	l=0;
	if (recPos!=0L)
	{
		while (1)
		{
			testKey=index.current_key();
			if ((c=atoi(testKey))==0 || c>spec.ec)
				break;
			if (c==spec.ec && atol(testKey+3)>spec.ep)
				break;
			load(tempRecord,recPos);
			sprintf(locusName[l],"%d:%ld",tempRecord.chr,tempRecord.pos);
			strncpy(locusInfo[l],tempRecord.myLocalLocus[0]->getInfo(),MAXSTR);
			locusInfo[l][MAXSTR]='\0';
			consequence[l][0]='\0';
			if (tempRecord.ensemblConsequence[0]!='\0')
				strncpy(consequence[l],tempRecord.ensemblConsequence,MAXSTR);
			else if (tempRecord.quickConsequence[0]!='\0')
				strncpy(consequence[l],tempRecord.quickConsequence,MAXSTR);
			consequence[l][MAXSTR]='\0';
			buff[0]='\0';
			for (i=0;i<tempRecord.nAlls;++i)
				sprintf(strchr(buff,'\0'),"%s%c",tempRecord.alls[i],i==tempRecord.nAlls-1?'\0':'/');
			strncpy(alleles[l],buff,MAXSTR);
			alleles[l][MAXSTR]='\0';
			recPos=index.get_next();
			++l;
			if (recPos==0L)
				break;
		}
	}
	lc=outputCalls(call,spec);
	assert((fp=fopen(fn,"w"))!=0);
	fprintf(fp,"Position\t");
	for (ll=0;ll<l;++ll)
		fprintf(fp,"%s\t",locusName[ll]);
	fprintf(fp,"\n");
	fprintf(fp,"Alleles\t");
	for (ll=0;ll<l;++ll)
		fprintf(fp,"%s\t",alleles[ll]);
	fprintf(fp,"\n");
	fprintf(fp,"Consequence\t");
	for (ll=0;ll<l;++ll)
		fprintf(fp,"%s\t",consequence[ll]);
	fprintf(fp,"\n");
	fprintf(fp,"Info\t");
	for (ll=0;ll<l;++ll)
		fprintf(fp,"%s\t",locusInfo[ll]);
	fprintf(fp,"\n");
	for (i=0;i<2;++i)
		if (useFreqs[i])
		{
			float *freqs;
			freqs=new float[lc];
			outputEurAltFrequencies(freqs,i,spec);
			fprintf(fp,"Frequency\t");
			for (l=0;l<lc;++l)
				fprintf(fp,"%6.4f\t",freqs[l]);
			fprintf(fp,"\n");
			delete[] freqs;
		}
	for (i=0;i<totalSub;++i)
	{
		fprintf(fp,"%s\t",subName[i]);
		for (ll=0;ll<l;++ll)
			fprintf(fp,"%s\t",call[ll][i]);
		fprintf(fp,"\n");
	}


	fclose(fp);
	for (i=0;i<totalSub;++i)
		free (call[i]);
	free(call);
	free(locusName);
	free(consequence);
	free(alleles);
	free(locusInfo);
	return 1;
}

int masterLocusFile::outputCalls(strEntry **call,analysisSpecs &spec)
{
	int locusCount,subCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c,i;
	locusCount=0;
recPos=findFirstInRange(spec);
if (recPos!=0L)
{
	while (1)
	{
		testKey=index.current_key();
		if ((c=atoi(testKey))==0 || c>spec.ec)
			break;
		if (c==spec.ec && atol(testKey+3)>spec.ep)
			break;
		load(tempRecord,recPos);

		for (subCount=0,i=0;i<nLocusFiles;++i)
		{
			tempRecord.outputCalls(call[locusCount]+subCount,locusFiles[i]->fp,i,nSubs[i],spec);
			subCount+=nSubs[i];
		}
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;

	}
}
return locusCount;
}

int masterLocus::outputCalls(strEntry *call,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec)
{
	int s;
	if (!strcmp(myLocalLocus[whichFile]->filter,"UNTYPED"))
	{
			for (s=0;s<nSubs;++s)
			{
				call[s][0]='.';
				call[s][0]='\0';
			}
	}
	else
		return myLocalLocus[whichFile]->outputCalls(call,f,locusPosInFile[whichFile],nSubs,alleleMapping[whichFile],spec);
return 1;
}

int vcfLocusFile::outputSubNames(strEntry *subName, analysisSpecs &spec)
{
	int s,sk;
	char *sPtr,*ptr;
	s=0;
	fseek(fp,0L,SEEK_SET);
	do {
		if (!fgets(locusFile::buff,BUFFSIZE-1,fp))
		{
			dcerror(99,"Could not find line beginning #CHROM in VCF file");
			return 0;
		}
	} while (strncmp(buff,"#CHROM",strlen("#CHROM")));
	for (ptr=buff,sk=0;sk<nFieldsToSkip;++sk)
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
		if (*ptr=='\0')
			break; // end of line with no entries
	}
	while (*ptr)
	{
		sPtr=subName[s];
		while (!isspace(*ptr))
			*sPtr++=*ptr++;
		*sPtr='\0';
		++s;
		while (isspace(*ptr))
			++ptr;
	}
return s;
}

int masterLocusFile::outputAffectionStatus(int *ccToWrite, analysisSpecs &spec)
{
	int i,s,ss;
	for (s=0,i=0;i<nLocusFiles;++i)
		for (ss=0;ss<nSubs[i];++s,++ss)
			ccToWrite[s]=spec.phenotypes?spec.phenotypes[s]:cc[i];
return 1;
}

int masterLocusFile::outputSubNames(strEntry *subName,analysisSpecs &spec)
{
	int totalSubs,f;
	totalSubs=0;

	if (locusFiles[0]->fp==0)
		openLocusFiles();
	for (f=0;f<nLocusFiles;++f)
	{
		totalSubs+=locusFiles[f]->outputSubNames(&subName[totalSubs],spec);
	}
return 1;
}


int vcfLocalLocus::outputCalls(strEntry *call,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec)
{
	char *ptr,*strPtr,firstStr[100],secondStr[100],allStr[100];
	allelePair all;
	int s,i;
	if (fseek(f,filePos,SEEK_SET)!=0)
	{
		dcerror(99,"Failed to fseek() correctly in localLocus::outputAlleles()");
		return 0;
	}
    if (!fgets(locusFile::buff,BUFFSIZE-1,f))
	{
		dcerror(99,"Failed read locus data after fseek() in localLocus::outputAlleles()");
		return 0;
	}
	for (s=0,ptr=locusFile::buff;s<nFieldsToSkip;++s)
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	for (s=0;s<nSubs;++s)
	{
		firstStr[0]=allStr[0]=secondStr[0]='\0';
		strPtr=firstStr;
		if (*ptr!='.') for (i=0;i<GTpos;++i)
		{
			while (*ptr!=':')
			{
				if (isspace(*ptr))
				{
					dcerror(99,"Not enough genotypes for number of subjects in this line: %s",locusFile::buff);
					return 0;
				}
				*strPtr++=*ptr++;
			}
			++ptr;
			*strPtr++=':';
		}
		*strPtr='\0';
		sscanf(ptr,"%[^: \t]%s",allStr,secondStr);
		if (allStr[0]!='.')
		{
			all[0]=alleleMap[allStr[0]-'0']+1;
			if (allStr[1]=='\0')
				sprintf(allStr,"%d",all[0]);
			else
			{
				all[1]=alleleMap[allStr[2]-'0']+1;
				sprintf(allStr,"%d/%d",all[0],all[1]);
			}
		}
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
		sprintf(call[s],"%s%s%s",firstStr,allStr,secondStr);
	}
	return 1;
}
