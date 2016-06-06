#include <stdio.h>
#include <string.h>
#include "dcerror.hpp"
#include "dcexpr.hpp"
#include "consequenceType.hpp"
extern "C"
{
#include "cdflib.h"
};

#define MAXGROUP 5
#define MAXPERGROUP 5
#define MAXFILTER 20

#define EVAL_R1 \
if ((r1=b1->eval())==NULL) return NULL; 

#define EVAL_R2 \
if ((r2=b2->eval())==NULL) return NULL; 

#define EVAL_BOTH \
if ((r1=b1->eval())==NULL) return NULL; \
if ((r2=b2->eval())==NULL) { delete r1; return NULL; }
double chistat(double x,double df);

double chistat(double x,double df)
{
double p,q,d1=1.0,bound;
int status,which=1;
if (x==0.0) return 1.0;
// if (x<0.0 && x>-1.0) return 1.0; 
if (x<0.0) return 1.0; 
/* do not worry about negative lrt values */
cdfchi(&which,&p,&q,&x,&df,&status,&bound);
if (status!=0)
	dcerror(1,"cdfchi failed","");
return q;
}


float cumulBinom(int N,int k,float p)
// return cumulative probability of getting k successes or fewer with N attempts at probability p of success
{
	double logExactP,cumulP,logP,log1mP,binom;
	int i;
	if (p==0 || k==N)
		return 1.0;
	if (p==1)
		return 0.0;
	logP=log(p);
	log1mP=log(1-p);
	binom=0;
	for (i=1;i<=k;++i)
		binom+=log((double)N-i+1)-log((double)i);
	logExactP=k*logP+(N-k)*log1mP+binom;
	cumulP=exp(logExactP);
	for (i=k;i>=1;--i)
	{
		logExactP+=log1mP-logP-log((double)N-i+1)+log((double)i);
		cumulP+=exp(logExactP);
	}
	return cumulP;
}

float two_x_two_chisq(float tab[2][2])
{
float ex[2][2],col_tot[2],row_tot[2],N;
int r,c;
float dchi;
N=dchi=0;
for (r=0;r<2;++r)
  row_tot[r]=0;
for (c=0;c<2;++c)
  col_tot[c]=0;
for (r=0;r<2;++r)
  for (c=0;c<2;++c)
    {
    row_tot[r]+=tab[r][c];
    col_tot[c]+=tab[r][c];
    N+=tab[r][c];
    }
for (r=0;r<2;++r)
  if (row_tot[r]==0 || col_tot[r]==0)
    return 0;
for (r=0;r<2;++r)
  for (c=0;c<2;++c)
    {
    ex[r][c]=row_tot[r]*col_tot[c]/N;
    dchi+=(tab[r][c]-ex[r][c])*(tab[r][c]-ex[r][c])/ex[r][c];
    }
if (tab[1][1]<ex[1][1])
	dchi*=-1;
return dchi;
}

int summCount[MAXGROUP][3],summTot[MAXGROUP],summFreq[MAXGROUP],consequenceType,chr,nFilter;
char consequenceStr[50],filterStr[MAXFILTER][200],geneName[30];
char *gargv[100];

express filter[MAXFILTER];

dcexpr_val *consequenceType_op(dcvnode *b1)
{
int c;
if (consequenceStr[0]=='\0')
	c=0;
else
	for (c=0;c<NCONSEQUENCETYPES;++c)
		if (!strncmp(consequenceStr,consequence[c].str,strlen(consequence[c].str)))
			break;
double rv=consequence[c].t;
return new dcexpr_double(rv);
}

dcexpr_val *nsub_op(dcvnode *b1)
{
dcexpr_val *r1;
int f;
EVAL_R1;
f=double(*r1);
double rv=summCount[f][0]+summCount[f][1]+summCount[f][2];
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *gene_op(dcvnode *b1)
{
	return new dcexpr_string(geneName);
}

dcexpr_val *chr_op(dcvnode *b1)
{
	return new dcexpr_double(chr);
}

dcexpr_val *countAA_op(dcvnode *b1)
{
	dcexpr_val *r1;
	int f;
	EVAL_R1;
	f = double(*r1);
	double rv = summCount[f][0];
	delete r1;
	return new dcexpr_double(rv);
}

dcexpr_val *countAB_op(dcvnode *b1)
{
dcexpr_val *r1;
int f;
EVAL_R1;
f=double(*r1);
double rv=summCount[f][1];
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *countBB_op(dcvnode *b1)
{
dcexpr_val *r1;
int f;
EVAL_R1;
f=double(*r1);
double rv=summCount[f][2];
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *argv_op(dcvnode *b1)
{
dcexpr_val *r1;
int a;
EVAL_R1;
a=double(*r1);
char *s=gargv[a-1];
delete r1;
return new dcexpr_string(s);
}

dcexpr_val *chicomp_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
int f,g;
float chisq,tab[2][2];
EVAL_BOTH;
f=double(*r1);
g=double(*r2);
if (chr==24) //Y
{
tab[0][0]=summCount[f][0];
tab[0][1]=summCount[f][1];
tab[1][0]=summCount[g][0];
tab[1][1]=summCount[g][1];
}
else
{
tab[0][0]=summCount[f][0]*2+summCount[f][1];
tab[0][1]=summCount[f][2]*2+summCount[f][1];
tab[1][0]=summCount[g][0]*2+summCount[g][1];
tab[1][1]=summCount[g][2]*2+summCount[g][1];
}
double rv=two_x_two_chisq(tab);
delete r1;
delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *recCHIcomp_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
int f,g;
float tab[2][2];
double rv;
EVAL_BOTH;
f=double(*r1);
g=double(*r2);
if (chr==24 || chr==23)
{
	rv=0;
}
else
{
tab[0][0]=summCount[f][0]+summCount[f][1];
tab[0][1]=summCount[f][2];
tab[1][0]=summCount[g][0]+summCount[g][1];
tab[1][1]=summCount[g][2];
rv=two_x_two_chisq(tab);
}
delete r1;
delete r2;
return new dcexpr_double(rv);
}

dcexpr_val *HWE_op(dcvnode *b1)
{
dcexpr_val *r1;
int f;
float o,e,N,p,chisq;
EVAL_R1;
f=double(*r1);
o=summCount[f][1];
N=summCount[f][0]+summCount[f][1]+summCount[f][2];
p=(summCount[f][0]+summCount[f][1]/2)/N;
e=N*2*p*(1-p);
chisq=(o-e)*(o-e)/e;
if (o>e)
	chisq*=-1;
double rv=chisq;
delete r1;
return new dcexpr_double(rv);
}

dcexpr_val *recHWEcomp_op(dcvnode *b1,dcvnode *b2)
{
dcexpr_val *r1,*r2;
int r[2],rr,c;
float freq,ex[2][2],tab[2][2],N,totCase;
double dchi,p,rv;

EVAL_BOTH;
r[0]=double(*r1);
r[1]=double(*r2);
if (chr==24 || chr==23)
	rv=0;
else
{
	N=freq=0;
	for (rr=1;rr>=0;--rr)
	{
		N+=summCount[r[rr]][0]+summCount[r[rr]][1]+summCount[r[rr]][2];
		if (rr==1)
			totCase=N; // save counting twice
		freq+=summCount[r[rr]][1]*0.5+summCount[r[rr]][2];
	}
	freq=freq/N;
	if (freq==0 || freq==1.0)
		rv=0;
	else
		// this code is just stolen from scoreassoc
	{
	tab[0][1]=summCount[1][2];
	ex[0][1]=totCase*freq*freq;
	tab[0][0]=totCase-tab[0][1];
	ex[0][0]=totCase-ex[0][1];
	if (ex[0][1]>5)
	{
		dchi=0;
		for (c=0;c<2;++c)
			if (ex[0][c]!=0 && tab[0][c]!=0)
				dchi+=(fabs(tab[0][c]-ex[0][c])-0.5)*(fabs(tab[0][c]-ex[0][c])-0.5)/(ex[0][c]<1?1:ex[0][c]);
		p=chistat(dchi,1.0)/2;
		if (ex[0][1]>tab[0][1])
			p=1-p;
	}
	else
	{
		if (tab[0][1]==0)
			p=1.0;
		else
			p=1-cumulBinom(totCase,tab[0][1]-1,freq*freq);
		if (p<=0)
			p=pow((double)10,(double)-20);
	}
	rv=-log10(p);
	}
}
delete r1;
delete r2;
return new dcexpr_double(rv);
}

int initFilters()
{
	int i;
	add_un_op("ARGV", argv_op);
	add_un_op("NSUB",nsub_op);
	add_un_op("COUNTAA",countAA_op);
	add_un_op("COUNTAB",countAB_op);
	add_un_op("COUNTBB", countBB_op);
	add_un_op("CONSEQUENCE", consequenceType_op);
	add_un_op("GENE", gene_op);
	add_un_op("CHR", chr_op);
	add_un_op("HWE", HWE_op);
	add_bin_op_next("CHICOMP", chicomp_op);
	add_bin_op_same("RECCHICOMP",recCHIcomp_op);
	add_bin_op_same("RECHWECOMP",recHWEcomp_op);
	for (i=0;i<nFilter;++i)
	{
		if (filter[i].parse(filterStr[i])==0)
			return 0;
	}
	return 1;
}

int resync(FILE *fg[MAXGROUP][MAXPERGROUP],long lineStart[MAXGROUP][MAXPERGROUP],int nGroup,int nInGroup[MAXGROUP])
{
	// a variant listed in one vcf file is not listed in another file
	// we hope this will be a rare occurrence
	// we hope that the variant positions are listed in ascending order
	int g,c,notSynced;
	long pos[MAXGROUP][MAXPERGROUP],maxPos;
	char line[200];
	maxPos=0L;
	// reset to start of bad lines
	for (g=0;g<nGroup;++g)
		for (c=0;c<nInGroup[g];++c)
		{
			fseek(fg[g][c],lineStart[g][c],SEEK_SET);
			pos[g][c]=-1L;
		}
	notSynced=1;
	while (notSynced)
	{
		notSynced=0;
		for (g=0;g<nGroup;++g)
			for (c=0;c<nInGroup[g];++c)
			{
				lineStart[g][c]=ftell(fg[g][c]);
				 // idea is to only read lines from files which have fallen behind
				if (pos[g][c]<maxPos)
				{
					notSynced=1;
					if (!fgets(line,199,fg[g][c]))
						return 0;
					sscanf(line,"%*ld %*s %*d %ld",&pos[g][c]);
					if (pos[g][c]>maxPos)
						maxPos=pos[g][c];
				}
			}
	}

	for (g=0;g<nGroup;++g)
		for (c=0;c<nInGroup[g];++c)
			fseek(fg[g][c],lineStart[g][c],SEEK_SET);
	return 1;
}

void outputGeneCounts(FILE *fgo,char *oldGenName,int nGroup,int countByGene[MAXGROUP][3],int summCount[MAXGROUP][3])
{
	int g,i;
	fprintf(fgo,"%s\t||\t",oldGenName);
	for (g = 0; g < nGroup; ++g)
	{
		fprintf(fgo, "%d\t%d\t%d\t||",countByGene[g][0],countByGene[g][1],countByGene[g][2]);
		for (i=0;i<3;++i)
			countByGene[g][i]=summCount[g][i];
	}
	fprintf(fgo,"\n");
}

int main(int argc,char *argv[])
{
	FILE *fi,*fo,*fgo,*fg[MAXGROUP][MAXPERGROUP];
	int genotypeCount[MAXGROUP][MAXPERGROUP][3],countByGene[MAXGROUP][3],g,c,i,first,nGroup,nInGroup[MAXGROUP],oneGc[3],failed;
	long pos,oldPos;
	long lineStart[MAXGROUP][MAXPERGROUP];
	char fn[30],line[200],rest[200],oldGeneName[30];
	float rv[MAXFILTER];
	fi=fopen(argv[1],"r");
	fo=fopen(argv[2],"w");
	fgo=fopen(argv[3],"w");
	for (i = 4; i < argc; ++i)
		gargv[i-4] = argv[i];
	fgets(line,199,fi);
	sscanf(line,"%d",&nGroup);
	for (g=0;g<nGroup;++g)
	{
		fgets(line,199,fi);
		nInGroup[g]=0;
		while (*rest='\0',sscanf(line,"%s %[^\n\r]",fn,rest)>=1)
		{
			fg[g][nInGroup[g]]=fopen(fn,"rb");
			++nInGroup[g];
			strcpy(line,rest);
		}
		for (i=0;i<3;++i)
			countByGene[g][i]=0;
	}
	fgets(line,199,fi);
	sscanf(line,"%d",&nFilter);
	for (c=0;c<nFilter;++c)
		fgets(filterStr[c],199,fi);
	fclose(fi);
	initFilters();
	first=1;
	oldGeneName[0]='\0';
	while (1)
	{
synchronised:
		oldPos=-1L;
		for (g=0;g<nGroup;++g)
			for (c=0;c<nInGroup[g];++c)
				lineStart[g][c]=ftell(fg[g][c]);

		for (g=0;g<nGroup;++g)
		{
			summCount[g][0]=summCount[g][1]=summCount[g][2]=0;
			for (c=0;c<nInGroup[g];++c)
			{
				if (!fgets(line,199,fg[g][c]))
				{
					if (oldPos==-1L)
						goto done;
					else
						dcerror(1,"Could not read line from fg[%d][%d]\n",g,c);
				}
			consequenceStr[0]='\0';
			sscanf(line,"%*ld %s %d %ld %d %d %d %[^\n\r]",geneName,&chr,&pos,&oneGc[0],&oneGc[1],&oneGc[2],consequenceStr);
			if (oldPos==-1L)
				oldPos=pos;
			else
				if (pos!=oldPos)
				{
					fprintf(stderr,"Positions do not match: %ld and %ld in line from fg[%d][%d]\n",pos,oldPos,g,c);
					if (!resync(fg,lineStart,nGroup,nInGroup))
						dcerror(1,"Failed to resynchronise positions between files.\n");
					goto synchronised;
				}
			for (i=0;i<3;++i)
				summCount[g][i]+=oneGc[i];
			}
		}
		failed=0;
		for (c=0;c<nFilter;++c)
			rv[c]=0;
		for (c=0;c<nFilter;++c)
			if ((rv[c]=double(*filter[c].eval()))==0)
			{
				failed=1;
				break;
			}
		if (failed==0)
		{
		fprintf(fo," %s\t%d\t%ld\t%s\t||\t",geneName,chr,pos,consequenceStr);
		for (g=0;g<nGroup;++g)
		{
			fprintf(fo, "%d\t%d\t%d\t||\t", summCount[g][0], summCount[g][1], summCount[g][2]);
		}
		for (c=0;c<nFilter;++c)
			fprintf(fo,"%6.2f\t",rv[c]);
		fprintf(fo,"\n");
		if (first)
			{
				strcpy(oldGeneName,geneName);
				first=0;
			}
		if (strcmp(geneName,oldGeneName))
		{
			outputGeneCounts(fgo,oldGeneName,nGroup,countByGene,summCount);
			strcpy(oldGeneName,geneName);
		}
		else
			for (g=0;g<nGroup;++g)
				for (i=0;i<3;++i)
					countByGene[g][i]+=summCount[g][i];
		}
	}
done:
	outputGeneCounts(fgo,oldGeneName,nGroup,countByGene,summCount);
	fclose(fgo);
	fclose(fo);
	return 0;
}
