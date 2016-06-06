#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dcerror.hpp"
#include "getGene.hpp"

void writeCounts(FILE *fo,char *gene,int chr,int *nCT,float p)
{
	int ct;
	fprintf(fo,"%s\t%d\t%d\t%d\t%f\t%d\t",gene,
		nCT[5]+nCT[10],
		nCT[5],
		nCT[10],
		(nCT[5]+nCT[10])?(nCT[10]-p*(nCT[5]+nCT[10]))/sqrt(p*(1-p)*(nCT[5]+nCT[10])):0.0,
		chr);
	for (ct=0;ct<NCONSEQUENCETYPES;++ct)
		fprintf(fo,"%d\t",nCT[ct]);
	fprintf(fo,"\n");
}

int main(int argc,char *argv[])
{
	FILE *fi,*fo;
	char consequenceStr[30],line[200],geneName[30],oldGeneName[30];
	int nConsequenceType[NCONSEQUENCETYPES],ct,chr,gen[3];
	float p,symCount,nonSymCount;
	fi=fopen(argv[1],"r");
	oldGeneName[0]='\0';
	symCount=nonSymCount=0;
	while (fgets(line,199,fi))
	{
		sscanf(line,"%*ld %*s %*d %*ld %d %d %d %s",&gen[0],&gen[1],&gen[2],consequenceStr);
		if (gen[1]+gen[2]==0)
			continue;
		for (ct=0;ct<NCONSEQUENCETYPES;++ct)
			if (!strcmp(consequenceStr,consequence[ct].str))
				break;
		if (ct==NCONSEQUENCETYPES)
				dcerror(0,"Could not match consequence type %s in line %s",consequenceStr,line);
		else if (consequence[ct].t==SYNONYMOUS_CODING)
			++symCount;
		else if (consequence[ct].t==NON_SYNONYMOUS_CODING)
			++nonSymCount;
	}
	fclose(fi);
	p=nonSymCount/(nonSymCount+symCount);
	fi=fopen(argv[1],"r");
	oldGeneName[0]='\0';
	fo=fopen(argv[2],"w");
	fprintf(fo,"gene\tCODING\tSYNONYMOUS_CODING\tNON_SYNONYMOUS_CODING\tZ\tchr\t",consequence[SYNONYMOUS_CODING].str,consequence[NON_SYNONYMOUS_CODING].str);
	for (ct=0;ct<NCONSEQUENCETYPES;++ct)
		fprintf(fo,"%s\t",consequence[ct].str);
	fprintf(fo,"\n");
	while (fgets(line,199,fi))
	{
		sscanf(line,"%*ld %s %d %*ld %d %d %d %s",geneName,&chr,&gen[0],&gen[1],&gen[2],consequenceStr);
		if (strcmp(geneName,oldGeneName))
		{
			if (oldGeneName[0]!='\0')
				writeCounts(fo,geneName,chr,nConsequenceType,p);
			strcpy(oldGeneName,geneName);
			for (ct=0;ct<NCONSEQUENCETYPES;++ct)
				nConsequenceType[ct]=0;
		}
		else

		{
			for (ct=0;ct<NCONSEQUENCETYPES;++ct)
				if (!strcmp(consequenceStr,consequence[ct].str))
					break;
			if (ct==NCONSEQUENCETYPES)
				dcerror(0,"Could not match consequence type %s in line %s",consequenceStr,line);
			else if (gen[1]+gen[2]!=0)
					++nConsequenceType[ct];
		}
	}
	writeCounts(fo,geneName,chr,nConsequenceType,p);
	fclose(fi);
	fclose(fo);
	return 0;
}
