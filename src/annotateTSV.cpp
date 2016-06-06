#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "geneVarUtils.hpp"
#define MAXALT 10
extern "C"
{
#include "cdflib.h"
};

#define EURTOT 395

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

float logfac(int h,int l)
 {
 float lf=0;
 if (h==l) return 0.;
 for (;h>l;--h) lf+=log((float)h);
 return lf;   /* ln(h!/l!) */
 }

float chisqfish(float tab[2][2])
  {
	  float ex[2][2],cell[3][3];
  int r=2,c=2,i,j,l,fisher;
  float test,xp,chisq,logfac_N,wrongWay,A,B,C,D;
  fisher=0;
  wrongWay=0;
  for (j=0;j<=c;++j) for (i=0;i<r+1;++i) cell[i][j]=0;
  for (i=0;i<r;++i)
   {
   for (j=0;j<c;++j)
     {
		cell[i][j]=tab[i][j];
     cell[i][c]+=cell[i][j];
     cell[r][j]+=cell[i][j];
     }
   }
  for (j=0;j<c;++j)
    cell[r][c]+=cell[r][j];
  chisq=0;
  for (i=0;i<r;++i)
   {
   for (j=0;j<c;++j)
     {
     xp=(cell[r][j]*cell[i][c])/cell[r][c];
	 if (i==0&&j==0&&xp>cell[i][j])
		 wrongWay=1;
     if (xp<5) fisher=1;
     test=cell[i][j]-xp;
     chisq+=test*test/xp;
     }
   }
  if (cell[r][c]<20) fisher=1;
  if (!fisher)
	  return -log10(chistat(chisq,1))*(wrongWay?-1:1);
  else
   {
   l=1000;
   for (i=0;i<2;++i)
    for (j=0;j<2;++j)
     if (cell[i][j]<l)
      {
      A=l=cell[i][j];
      B=cell[(i+1)%2][j];
      C=cell[i][(j+1)%2];
      D=cell[(i+1)%2][(j+1)%2];
      }
   logfac_N=logfac(A+B+C+D,0);
   chisq=0;
   i=A;D-=A;B+=A;C+=A;
   for (A=0;A<=i;++A,++D,--B,--C)
    {
    xp=logfac(A+B,B)+logfac(A+C,A)+logfac(B+D,D)+logfac(C+D,C)-logfac_N;
    chisq+=exp(xp);
    }
   if (chisq>.5) chisq=1-chisq+exp(xp);
//   outres("\nFisher's exact test, p = %.3g\n",chisq);
   return -log10(chisq)*(wrongWay?-1:1);
   }
  }

void get1000GCounts(int chr,int pos,char *rsname,float *AC,float *AN,float *EUR_AF)
{
	char line[1000],chrStr[3];
	FILE *fp;
	float ac,an,eur_AF,highestFreq;
	unlink("get1000Gfreq.vcf");
	if (chr<23)
		sprintf(chrStr,"%d",chr);
	else
		sprintf(chrStr,"%c",chr==23?'X':chr==24?'Y':'?');
	sprintf(line,"tabix \\reference\\ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.gz %s:%d-%d > get1000Gfreq.vcf",chrStr,pos,pos);
	printf("%s\n",line);
	system(line);
	ac=0;
	an=2184;
	eur_AF=0;
	*AC=ac;
	*AN=an;
	*EUR_AF=eur_AF;
	highestFreq=0;
	rsname[0]='\0';
	if ((fp=fopen("get1000Gfreq.vcf","r"))!=0)
	{
		while (fgets(line,999,fp))
		if (strstr(line,"AN=") && strstr(line,"AC="))
		{
			sscanf(line,"%*s %*s %s",rsname);
			sscanf(strstr(line,"AN="),"AN=%f",&an);
			sscanf(strstr(line,"AC="),"AC=%f",&ac);
			if (strstr(line,"EUR_AF="))
				sscanf(strstr(line,"EUR_AF="),"EUR_AF=%f",&eur_AF);
			if (ac/an>highestFreq || eur_AF>highestFreq)
			{
				highestFreq=ac/an;
				if (eur_AF>highestFreq)
					highestFreq=eur_AF;
				*AC=ac;
				*AN=an;
				*EUR_AF=eur_AF;
			}
		}
		fclose(fp);
	}
}


int main(int argc,char *argv[])
{
	refseqGeneInfo r;
	char fn[100],fn2[100],line[20001],geneName[100],*fno,rsname[100];
	char dpStr[1001],chrStr[101],stStr[101],enStr[101],refStr[1001],all0Str[1001],all1Str[MAXALT][1001],geneStr[101],freqStr[101],featureStr[101],*ptr;
	int i,genoCount[3],pos,nAll,a;
	FILE *fp,*fo,*fpvo,*fvepo;
	float G1000MAF,AC,AN,tab[2][2],EUR_AF,EUR_AC,EUR_AN;
	gvaParams gp;
	analysisSpecs spec;
	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
		r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setUpstream(1000);
	r.setDownstream(1000);
	r.setBaitMargin(gp.margin);
	fp=fopen(argv[2],"r");
	fo=fopen(argv[3],"w");
	strcpy(fn2,argv[3]);
	strcpy(strchr(fn2,'.'),".pv.inp");
	fpvo=fopen(fn2,"w");
	strcpy(strchr(fn2,'.'),".vep.inp");
	fvepo=fopen(fn2,"w");
	fprintf(fo,"Chr\tPos\tName\tGene\tEffect\tFreq\tAA\tAB\tBB\tMAF\tAC\tAN\t1000GMAF\tEUR_AF\tMLP\teurMLP\tUCSC\tdbSNP\tVEP\t");
	fgets(line,1000,fp);
	fprintf(fo,"%s",line);
	while (*geneStr='\0',fgets(line,20000,fp))
	{
		int f;
		ptr=line;
		f=0;
		do {
			switch (f)
			{
			case 1: sscanf(ptr,"\t'%[^']'",dpStr); break;			 
			case 2: sscanf(ptr,"\t'%[^']'",chrStr); break;
			case 3: sscanf(ptr,"\t%s",stStr); break;
			case 4: sscanf(ptr,"\t%s",enStr); break;
			case 5: sscanf(ptr,"\t'%[^']'",refStr); break;
			case 6: nAll=sscanf(ptr,"\t'%[^',],%[^',],%[^',],%[^',],%[^',],%[^',],%[^',],%[^',],%[^',],%[^',]'",
						all0Str,
						all1Str[0],
						all1Str[1],
						all1Str[2],
						all1Str[3],
						all1Str[4],
						all1Str[5],
						all1Str[6],
						all1Str[7],
						all1Str[8],
						all1Str[9]); 
				break;
			case 8: sscanf(ptr,"\t%s",freqStr); break;
			case 12: sscanf(ptr,"\t'%[^:']'",geneStr); break;
			default: ;
			}
			++f;
		} while (ptr!=NULL && (ptr=strchr(ptr+1,'\t'))!=NULL);
		// sscanf(line,"%*s '%[^']' '%[^']' %s %s %*s '%[^,],%[^']'\t%*[^\t]\t%*s %*s '%[^:']", // the %*[^\t] argument is for a field which may be blank
		//	dpStr,chrStr,stStr,enStr,all0Str,all1Str,geneStr)>=6)
		if (strcmp(refStr,all0Str))
		{
			for (i=0;i<nAll-1;++i)
				if (!strcmp(refStr,all1Str[i]))
				{
					strcpy(all1Str[i],all0Str);
					strcpy(all0Str,refStr);
					break;
				}
			if (i==nAll-1)
				// sometimes the ref allele is not listed !!
			{
				strcpy(all1Str[nAll-1],all0Str);
				strcpy(all0Str,refStr);
				++nAll;
			}
		}
		pos=(atoi(stStr)+1>atoi(enStr))?atoi(enStr):atoi(stStr)+1;
		if ((ptr=strchr(line,'\n'))!=0)
			*ptr='\0';
		if ((ptr=strchr(line,'\r'))!=0)
			*ptr='\0';
		//	fprintf(fp,"Chr\tPos\tGene\tEffect\tFreq\tAA\tAB\tBB\tUCSC\tdbSNP\t");
		get1000GCounts(chrStr[0]=='Y'?24:chrStr[0]=='X'?23:atoi(chrStr),all0Str[0]=='-'?atoi(enStr)+1:atoi(stStr)+1,rsname,&AC,&AN,&EUR_AF);
		fprintf(fo,"%s\t%d\t%s\t%s\t",chrStr,pos,rsname,geneStr);
		dcerror.warn(); //  we do not mind if cannot find gene
		if (*geneStr && r.findGene(geneStr) && r.getNextGene() && r.getEffect(atoi(stStr)+1,all0Str,all1Str[0],r.getUpstream(),r.getDownstream(),1)!=NULL_CONSEQUENCE)
			strcpy(featureStr,r.tellEffect());
		else
			strcpy(featureStr,"UNKNOWN");
		dcerror.kill();
		fprintf(fo,"%s\t",featureStr);
		fprintf(fo,"%s\t",freqStr);
		for (f=0;f<3;++f)
			genoCount[f]=0;
		for (ptr=dpStr;*ptr;++ptr)
		{
			++genoCount[*ptr-'0'];
			if (*ptr=='+') // indel, not sure what + means but not reference
				++genoCount[1];
		}
		for (f=0;f<3;++f)
			fprintf(fo,"%d\t",genoCount[f]);
		fprintf(fo,"%.3f\t",(genoCount[1]*0.5+genoCount[2])/(genoCount[0]+genoCount[1]+genoCount[2]));
		G1000MAF=AC/AN;
		fprintf(fo,"%.0f\t%.0f\t%.4f\t%.3f\t",AC,AN,G1000MAF,EUR_AF);
		tab[0][0]=AN-AC;
		tab[0][1]=AC;
		tab[1][0]=genoCount[0]*2+genoCount[1];
		tab[1][1]=genoCount[1]+genoCount[2]*2;
		fprintf(fo,"%.2f\t",chisqfish(tab));
		EUR_AN=EURTOT*2;
		EUR_AC=floor(EUR_AN*EUR_AF+0.5);
		tab[0][0]=EUR_AN-EUR_AC;
		tab[0][1]=EUR_AC;
		tab[1][0]=genoCount[0]*2+genoCount[1];
		tab[1][1]=genoCount[1]+genoCount[2]*2;
		fprintf(fo,"%.2f\t",chisqfish(tab));
		fprintf(fo,"chr%s:%d-%s\t",chrStr,all0Str[0]=='-'?atoi(enStr)+1:atoi(stStr)+1,enStr);
		// (10[Chromosome]) AND (61842497[Base Position]) AND ("homo sapiens"[Organism])
		fprintf(fo,"http://www.ncbi.nlm.nih.gov/snp/?term=(%s[Chromosome]) AND (%d[Base Position]) AND (\"homo sapiens\"[Organism])\t",chrStr,atoi(stStr)+1);
		fprintf(fo,"%s %d %s %s/%s +\t",chrStr,all0Str[0]=='-'?atoi(enStr)+1:atoi(stStr)+1,enStr,all0Str,all1Str[0]);
		fprintf(fo,"%s\t",line);
		for (ptr=dpStr;*ptr;++ptr)
			fprintf(fo,"%c\t",*ptr);
		fprintf(fo,"\n");
		// seems like vep will want both start and end one-coded
		// "An insertion (of any size) is indicated by start coordinate = end coordinate + 1."
		for (a=0;a<nAll-1;++a)
		{
			//fprintf(fvepo,"%s\t%s\t%d\t%s/%s\t+\n",chrStr,stStr,all0Str[0]=='-'?atoi(stStr)-1:atoi(stStr)+strlen(all0Str)-1,all0Str,all1Str[i]);
			fprintf(fvepo,"%s\t%d\t%s\t%s/%s\t+\n",chrStr,all0Str[0]=='-'?atoi(enStr)+1:atoi(stStr)+1,enStr,all0Str,all1Str[a]);
			// use end+1 if deletion else convert start from 0-coded to 1-coded
		}
		if (all0Str[0]=='-')
			all0Str[0]='.';
		for (a=0;a<nAll-1;++a)
		{
			if (all1Str[a][0]=='-')
				all1Str[a][0]='.';
			fprintf(fpvo,"%s,%s,%s,%s\n",chrStr,stStr,all0Str,all1Str[a]);
			// PROVEAN does not code insertion in same way as Knome. I would have to look at the base before the start of the insertion and use this as the reference allele, then with the inserted bases appended to it for alt.
		}
	}
	fclose(fp);
	fclose(fo);
	fclose(fvepo);
	fclose(fpvo);
	return 0;
}