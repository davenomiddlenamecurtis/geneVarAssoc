#include <stdio.h>
#include <string.h>
#include "dcerror.hpp"
#include "getGene.hpp"

enum { AUT=0,CHRX,CHRY,ALLCHR };

#define NMAFTYPE 5
float fThreshold[NMAFTYPE]= { 0.001,0.01,0.1,0.2,0.6 }; 

long summGenoCount[NCONSEQUENCETYPES+1][4][3],chrTotal[3],total[4];
enum { ALLTYPES=NCONSEQUENCETYPES};
int nType[NCONSEQUENCETYPES+1][4],MAFType[NCONSEQUENCETYPES+1][NMAFTYPE];

int main(int argc,char *argv[])
{
	FILE *fi,*fo;
	int gc[3],chr,whichChr,ct,g,i;
	char consequenceStr[30],line[200];
	float fTotal[4],MAF,N;
	fi=fopen(argv[1],"r");
	fo=fopen(argv[2],"w");
	while (fgets(line,199,fi))
	{
		sscanf(line,"%*ld %*s %d %*ld %d %d %d %s",&chr,&gc[0],&gc[1],&gc[2],consequenceStr);
		if (gc[1]+gc[2]==0)
			continue;
		for (ct=0;ct<NCONSEQUENCETYPES;++ct)
		{
			if (!strcmp(consequenceStr,consequence[ct].str))
				break;
		}
		if (ct==NCONSEQUENCETYPES)
			dcerror(0,"Could not match consequence type %s in line %s",consequenceStr,line);
		whichChr=chr==23?CHRX:chr==24?CHRY:AUT;
		for (g=0;g<3;++g)
		{
			summGenoCount[ct][whichChr][g]+=gc[g];
			summGenoCount[ct][ALLCHR][g]+=gc[g];
			summGenoCount[ALLTYPES][whichChr][g]+=gc[g];
			summGenoCount[ALLTYPES][ALLCHR][g]+=gc[g];
			total[whichChr]+=gc[g];
			total[ALLCHR]+=gc[g];
		}
			++nType[ct][whichChr];
			++nType[ALLTYPES][whichChr];
			++nType[ct][ALLCHR];
			++nType[ALLTYPES][ALLCHR];
			if (whichChr==AUT && (N=(gc[2]+gc[1]+gc[0]))!=0)
			{
				MAF=(gc[2]+0.5*gc[1])/(N);
				if (MAF>0.5)
					MAF=1.0-MAF;
				for (i=0;i<NMAFTYPE;++i)
					if (MAF<fThreshold[i])
					{
						++MAFType[ct][i];
						++MAFType[ALLTYPES][i];
						break;
					}
			}

	}
	for (whichChr=0;whichChr<4;++whichChr)
		fTotal[whichChr]=total[whichChr];
	for (ct=0;ct<NCONSEQUENCETYPES+1;++ct)
	{
		if (ct==NCONSEQUENCETYPES)
			fprintf(fo,"%23s ","ALL");
		else
			fprintf(fo,"%23s ",consequence[ct].str);
		for (whichChr=0;whichChr<4;++whichChr)
		{
			fprintf(fo,"%8d ",nType[ct][whichChr]);
		}
		for (whichChr=0;whichChr<4;++whichChr)
		{
			N=summGenoCount[ct][whichChr][0]+summGenoCount[ct][whichChr][1]+summGenoCount[ct][whichChr][2];
			for (g=0;g<3;++g)
				fprintf(fo,"%8.6f ",(N==0)?0.0:summGenoCount[ct][whichChr][g]/N);
		}
		for (i=0;i<NMAFTYPE;++i)
			fprintf(fo,"%8.6f ",(nType[ct][AUT]==0)?0.0:MAFType[ct][i]/(float)nType[ct][AUT]);
		fprintf(fo,"\n");
	}
	
	fclose(fi);
	fclose(fo);
	return 0;
}
