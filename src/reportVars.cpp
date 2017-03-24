// \sequence\sharedseq\pathVars\composite.genos.out stopVars10.cond stopVars10.out
#include "masterLocusFile.hpp"
#include "dcexpr.hpp"
#include <assert.h>
#include <ctype.h>
#include <string.h>

#define MAXINTVARS 20000
#define MAXSUB 6000

#define MAXCOND 20
#define MAXCONDLENGTH 200

char condStr[MAXCOND][MAXCONDLENGTH];

typedef masterLocus *masterLocusPtr; 

masterLocusPtr currentLocus;

express condition[MAXCOND];
#define EVAL_R1 \
if ((r1=b1->eval())==NULL) return NULL; 
#define EVAL_R1_R2 \
if ((r1=b1->eval())==NULL || (r2=b2->eval())==NULL) return NULL; 

dcexpr_val *consequence_op(dcvnode *b1)
{
dcexpr_val *r1;
double rv;
EVAL_R1; // but ignore it
rv=double(currentLocus->getWorstConsequenceType());
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *genoCount_op(dcvnode *b1)
{
dcexpr_val *r1;
int c;
EVAL_R1;
c=double(*r1);
double rv=currentLocus->getGenoCount(c);
delete r1;
return new dcexpr_double(rv);
}

void addExtraFuncs()
{
	add_un_op("EFFECT",consequence_op);
	add_un_op("GCOUNT",genoCount_op);
}

int initConds(FILE *fc)
{
	char buff[1000];
	int nCond;
	for (nCond = 0;fgets(buff, 999, fc) && sscanf(buff, "%[^\n]", condStr[nCond])==1;++nCond)
	{
		assert(nCond<MAXCOND);
		if (condition[nCond].parse(condStr[nCond])==0)
			return -1;
	}
	return nCond;
}

#define LINELENGTH MAXSUB*20+1000
char line[LINELENGTH+1];
int main(int argc, char *argv[])
{
	int nCond,nIntVars,genoOffset,nSub,s,v;
	allelePair **subGenos;
	FILEPOSITION currentFilePos;
	FILE *fgenos,*fcond,*frep,*fsub;
	char *ptr,*sptr,ccBuff[10];
	strEntry *subName;
	int *cc;
	masterLocusPtr *vars;
	assert((fcond=fopen(argv[2],"r"))!=0);
	addExtraFuncs();
	nCond=initConds(fcond);
	fclose(fcond);
	if (nCond==-1) 
		return 1;
	vars=new masterLocusPtr[MAXINTVARS];
	assert((fgenos=fopen(argv[1],"rb"))!=0);
	nIntVars = 0; 
	currentLocus=new masterLocus(1);
	fgets(line, LINELENGTH, fgenos);
	fgets(line, LINELENGTH, fgenos);
	while (currentFilePos=FTELL(fgenos),fgets(line, LINELENGTH, fgenos))
	{
		bool useLocus=false;
		int c;
		currentLocus->readRiskVarInfo(line,true);
		for (c=0;c<nCond;++c)
			if (double(*condition[c].eval()) != 0)
			{
			useLocus=true;
			break;
			}
		if (useLocus)
		{
			currentLocus->setLocusPosInFile(0,currentFilePos);
			vars[nIntVars++]=currentLocus;
			currentLocus=new masterLocus(1);
		}
	}
	// now read in genotypes for all interesting variants and assign them to subjects
	FSEEK(fgenos,0L,SEEK_SET);
	fgets(line, LINELENGTH, fgenos);
	fgets(line, LINELENGTH, fgenos);
	fgets(line, LINELENGTH, fgenos);
	ptr=line;
	genoOffset=0;
	while (!isspace(*ptr))
	{
		++ptr;
		++genoOffset;
	}
	while (isspace(*ptr))
	{
		++ptr;
		++genoOffset;
	}
	while (!isspace(*ptr))
	{
		++ptr;
		++genoOffset;
	}
	while (isspace(*ptr))
	{
		++ptr;
		++genoOffset;
	}
	nSub=0;
	while (!isspace(*ptr) && *ptr) // this means line has to be formatted exactly correctly
	{
		ptr+=4;
		++nSub;
	}
	assert((subName=(strEntry*)calloc(nSub,sizeof(strEntry)))!=0);
	assert((cc=(int*)calloc(nSub,sizeof(int)))!=0);
	FSEEK(fgenos,0L,SEEK_SET);
	fgets(line, LINELENGTH, fgenos);
	ptr=line;
	for (s = 0; s < nSub; ++s)
	{
		while (isspace(*ptr))
			++ptr;
		sptr = subName[s];
		while (!isspace(*ptr))
			*sptr++ = *ptr++;
		*sptr='\0';
	}
	fgets(line, LINELENGTH, fgenos);
	ptr=line;
	for (s = 0; s < nSub; ++s)
	{
		while (isspace(*ptr))
			++ptr;
		sptr = ccBuff;
		while (!isspace(*ptr))
			*sptr++ = *ptr++;
		*sptr='\0';
		cc[s]=atoi(ccBuff);
	}

	assert((subGenos=(allelePair **)calloc(nSub,sizeof(allelePair *)))!=0);
	for (s=0;s<nSub;++s)
		assert((subGenos[s]=(allelePair*)calloc(nIntVars,sizeof(allelePair)))!=0);
	for (v = 0; v < nIntVars; ++v)
	{
		currentLocus=vars[v];
		FSEEK(fgenos,currentLocus->getLocusPosInFile(0),SEEK_SET);
		fgets(line, LINELENGTH, fgenos);
		ptr=line+genoOffset;
		for (s = 0; s < nSub; ++s)
		{
			subGenos[s][v][0]=ptr[0]-'0';
			subGenos[s][v][1]=ptr[2]-'0';
			ptr+=4;
		}
	}
	assert((frep=fopen(argv[3],"w"))!=0);
	for (s = 0; s < nSub; ++s)
	{
		bool first;
		first=true;
		fprintf(frep,"%s\t%d\t",subName[s],cc[s]);
		for (v = 0; v < nIntVars; ++v)
		{
			char info[1000];
			if (subGenos[s][v][1] >1)
			{
				vars[v]->writeRiskVarInfo(info,true);
				fprintf(frep,"[%d][%d][%d/%d]%s\t",s,v,subGenos[s][v][0],subGenos[s][v][1],info);
				first=false;
			}
		}
		fprintf(frep,"\n");
	}
	fclose(frep);
	if (argc > 4)
	{
	assert((frep=fopen(argv[4],"w"))!=0);
	for (s = 0; s < nSub; ++s)
	{
		bool first;
		first=true;
		fprintf(frep,"%s\t%d\t",subName[s],cc[s]);
		for (v = 0; v < nIntVars; ++v)
		{
			char info[1000],*ptr;
			int t;
			if (subGenos[s][v][1] >1)
			{
				ptr=0;
				vars[v]->writeRiskVarInfo(info,true);
				for (t=0;t<NCONSEQUENCETYPES;++t)
					if ((ptr = strstr(info, consequence[t].str)) != 0)
					{
					// pull out just the gene name
					ptr += strlen(consequence[t].str)+1;
					*(strchr(ptr, '_')) = '\0';
					break;
					}
				fprintf(frep,"[%d/%d]%s\t",subGenos[s][v][0],subGenos[s][v][1],ptr?ptr:info);
				first=false;
			}
		}
		fprintf(frep,"\n");
	}
	fclose(frep);
	}

	free(subName);
	free(cc);
	for (s=0;s<nSub;++s)
		free(subGenos[s]);
	free(subGenos);
	return 0;
}
