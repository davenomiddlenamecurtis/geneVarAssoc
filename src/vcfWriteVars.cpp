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

#include "vcfLocusFile.hpp"
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#define MAXLOCIINGENE 2000

#if 0
int masterLocusFile::writeVars(char *fn,int *useFreqs,analysisSpecs &spec)
{
	char buff[1000];
	int i,l,ll,totalSub,c,lc,*subjectSelector,a;
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
	assert ((subjectSelector=(int *)calloc(totalSub,sizeof(int)))!=0);

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
#endif

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
	char *sPtr,*ptr,ch;
	long fPos;
	s=0;
	fseek(fp,0L,SEEK_SET);
	do {
		fPos = ftell(fp); // I hope this will work OK with text file
		if (!fgets(locusFile::buff,BUFFSIZE-1,fp))
		{
			dcerror(99,"Could not find line beginning #CHROM in VCF file");
			return 0;
		}
	} while (strncmp(buff,"#CHROM",strlen("#CHROM")));
	fseek(fp, fPos, SEEK_SET);

	for (ch=fgetc(fp),sk=0;sk<nFieldsToSkip;++sk)
	{
		while (!isspace(ch))
			ch = fgetc(fp);
		while (isspace(ch))
		{
			if (ch == '\n')
				break; 
			ch = fgetc(fp);
		}
		if (ch=='\n')
			break; // end of line with no entries
	}
	while (ch!='\n')
	{
		sPtr=subName[s];
		while (!isspace(ch))
		{
			*sPtr++ = ch;
			ch = fgetc(fp);
		}
		*sPtr='\0';
		++s;
		while (isspace(ch))
		{
			if (ch == '\n')
				break;
			ch = fgetc(fp);
		}
	}
return s;
}

int masterLocusFile::outputAffectionStatus(float *ccToWrite, analysisSpecs &spec)
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
