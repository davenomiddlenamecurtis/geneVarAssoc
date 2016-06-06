#include "masterLocusFile.hpp"
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#define MAXLOCIINGENE 2000

int masterLocusFile::writeGenos(char *fn,int *useFreqs,analysisSpecs &spec)
{
	char buff[1000];
	int i,l,ll,totalSub,c,lc;
	long recPos;
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
			sprintf(locusName[l],"%ld",tempRecord.pos);
			strncpy(locusInfo[l],tempRecord.myLocalLocus[0].info,MAXSTR);
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
			if (call[ll][i][0]!='1' && call[ll][i][0]!='2')
				fprintf(fp,"-1\t");
			else
				fprintf(fp,"%d\t",(call[ll][i][0]=='2')+(call[ll][i][2]=='2'));
		 // will give males with minor allele genotype 1 for X-linked marker
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
