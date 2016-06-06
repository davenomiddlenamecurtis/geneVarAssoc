#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MAXVARS 2000
#define MAXSTR 100
#define MAXSUB 5000
#define BUFFSIZE 900000
typedef char strEntry[MAXSTR+1];

char *buff;

int main(int argc,char *argv[])
{
	FILE *fi,*fo,*fv;
	char oldGeneName[50],geneName[50],line[200],posBuff[100],*ptr;
	int v,s,nVar,nSub,i,colToSkip,c,chr;
	long pos,p;
	buff=(char *)malloc(BUFFSIZE);
	strEntry *locusName,*consequence,*alleles,*locusInfo,*subName,*gene,**call;
	assert ((gene=(strEntry *)calloc(MAXVARS,sizeof(strEntry)))!=0);
	assert ((locusName=(strEntry *)calloc(MAXVARS,sizeof(strEntry)))!=0);
	assert ((consequence=(strEntry *)calloc(MAXVARS,sizeof(strEntry)))!=0);
	assert ((alleles=(strEntry *)calloc(MAXVARS,sizeof(strEntry)))!=0);
	assert ((locusInfo=(strEntry *)calloc(MAXVARS,sizeof(strEntry)))!=0);
	assert ((subName=(strEntry *)calloc(MAXSUB,sizeof(strEntry)))!=0);
	assert ((call=(strEntry **)calloc(MAXSUB,sizeof(strEntry *)))!=0);
	for (i=0;i<MAXSUB;++i)
		assert ((call[i]=(strEntry *)calloc(MAXVARS,sizeof(strEntry)))!=0);
	v=0;
	oldGeneName[0]='\0';
	assert((fi=fopen(argv[2],"r"))!=0);
	while (fgets(line,199,fi) && sscanf(line,"%s %s",geneName,&posBuff)==2)
	{
		if (strcmp(oldGeneName,geneName))
		{
			sprintf(line,"\\msvc\\vcf\\geneShowVars %s %s",argv[1],geneName);
			system(line);
			strcpy(oldGeneName,geneName);
		}
		if ((ptr=strchr(posBuff,':'))!=0)
			++ptr;
		else
			ptr=posBuff;
		sscanf(ptr,"%ld",&pos);
		strcpy(gene[v],geneName);
		sprintf(line,"gsv.%s.txt",geneName);
		assert((fv=fopen(line,"r"))!=0);
		fgets(buff,BUFFSIZE-1,fv);
		for (ptr=strchr(buff,'\t'),colToSkip=0;sscanf(ptr,"%d:%ld",&chr,&p),p!=pos;++colToSkip)
		{
			ptr=strchr(ptr+1,'\t');
			assert(ptr!=0);
		}
		sprintf(locusName[v],"%d:%ld",chr,pos);
		fgets(buff,BUFFSIZE-1,fv);
		for (ptr=strchr(buff,'\t'),c=0;c<colToSkip;++c)
			ptr=strchr(ptr+1,'\t');
		sscanf(ptr,"%s",alleles[v]);
		fgets(buff,BUFFSIZE-1,fv);
		for (ptr=strchr(buff,'\t'),c=0;c<colToSkip;++c)
			ptr=strchr(ptr+1,'\t');
		sscanf(ptr,"%s",consequence[v]);
		fgets(buff,BUFFSIZE-1,fv);
		for (ptr=strchr(buff,'\t'),c=0;c<colToSkip;++c)
			ptr=strchr(ptr+1,'\t');
		sscanf(ptr,"%s",locusInfo[v]);
		s=0;
		while(fgets(buff,BUFFSIZE-1,fv) && sscanf(buff,"%s",subName[s]))
		{
			for (ptr=strchr(buff,'\t'),c=0;c<colToSkip;++c)
				ptr=strchr(ptr+1,'\t');
			sscanf(ptr+1,"%[^:]",call[s][v]);
			++s;
		}
		fclose(fv);
		++v;
	}
	fclose(fi);
	nSub=s;
	nVar=v;
	assert((fo=fopen(argv[3],"w"))!=0);
	fprintf(fo,"Gene\t");
	for (v=0;v<nVar;++v)
		fprintf(fo,"%s\t",gene[v]);
	fprintf(fo,"\n");
	fprintf(fo,"Position\t");
	for (v=0;v<nVar;++v)
		fprintf(fo,"%s\t",locusName[v]);
	fprintf(fo,"\n");
	fprintf(fo,"Alleles\t");
	for (v=0;v<nVar;++v)
		fprintf(fo,"%s\t",alleles[v]);
	fprintf(fo,"\n");
	fprintf(fo,"Consequence\t");
	for (v=0;v<nVar;++v)
		fprintf(fo,"%s\t",consequence[v]);
	fprintf(fo,"\n");
#if 0
	fprintf(fo,"Info\t");
	for (v=0;v<nVar;++v)
		fprintf(fo,"%s\t",locusInfo[v]);
	fprintf(fo,"\n");
#endif
	for (s=0;s<nSub;++s)
	{
		fprintf(fo,"%s\t",subName[s]);
		for (v=0;v<nVar;++v)
			fprintf(fo,"'%s'\t",call[s][v]);
		fprintf(fo,"\n");
	}
	fclose(fo);
	free(buff);
	return 0;
}
