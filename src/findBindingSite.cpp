#include <stdlib.h>
#include "geneVarUtils.hpp"

#define MAXSITES 5
#define MAXSITELENGTH 20
#define MAXWHOLESITELENGTH 200

class bindingSite {
	char site[MAXSITES][MAXSITELENGTH],minusSite[MAXSITES][MAXSITELENGTH],
	wholeSite[MAXWHOLESITELENGTH],minusWholeSite[MAXWHOLESITELENGTH];
	int nSites,pos[MAXSITES];
public:
	void clear() { nSites=0; }
	bindingSite() { nSites=0; }
	void addSite(char *s,int p);
	void addWholeSite(char *s);
	int getSite(char *result,char *toSearch,int isMinus);
};

void bindingSite::addSite(char *s,int p)
{
	strcpy(site[nSites],s);
	getComplementarySequence(minusSite[nSites],site[nSites]);
	pos[nSites]=p;
	++nSites;
}

void bindingSite::addWholeSite(char *s)
{
	strcpy(wholeSite,s);
	getComplementarySequence(minusWholeSite,wholeSite);
}

int bindingSite::getSite(char *res,char *toSearch,int isMinus)
{
	char *p,tempRes[MAXWHOLESITELENGTH],tempSeq[MAXWHOLESITELENGTH],*test;
	int foundOne,s;
	foundOne=0;
	p=toSearch;
	*res='\0';
	if (isMinus==0)
	{
	while ((p=strstr(p,site[0]))!=0)
	{
		for (s=1;s<nSites;++s)
			if (strncmp(p+pos[s]-pos[0],site[s],strlen(site[s]))!=0)
				goto noMatch;
		foundOne=1;
		strncpy(tempRes,p-pos[0],strlen(wholeSite));
		sprintf(tempRes+strlen(wholeSite),"\t%d\t",p-toSearch-pos[0]);
		strcat(res,tempRes);
noMatch:
		++p;
	}
	}
	else
	{
	while ((p=strstr(p,minusSite[0]))!=0)
	{
		for (s=1;s<nSites;++s)
		{
			test=p-pos[s]+pos[0]+strlen(minusSite[0])-strlen(minusSite[s]);
			if (strncmp(p-pos[s]+pos[0]+strlen(minusSite[0])-strlen(minusSite[s]),minusSite[s],strlen(minusSite[s]))!=0)
				goto noMinusMatch;
		}
		foundOne=1;
		strncpy(tempSeq,p+pos[0]+strlen(site[0])-strlen(wholeSite),strlen(wholeSite));
		tempSeq[strlen(wholeSite)]='\0';
		getComplementarySequence(tempRes,tempSeq);
		sprintf(tempRes+strlen(wholeSite),"\t%d\t",p-toSearch+pos[0]+strlen(site[0])-strlen(wholeSite));
		strcat(res,tempRes);
noMinusMatch:
		++p;
	}
	}
	return foundOne;
}

void lookFormiRNA137(bindingSite &b)
{
	char temp[100],minusSequence[100];
	b.clear();
	             // 5 UAAUAUGUAUCAGCUAGCAAUAU 3 sequence in paper, unclear if this is RNA or DNA
	             // 5                AGCAAUA  3 first site 
	             // 5       GUAU              3 second site
	             // 5 ATATTGCTAGCTGATACATATTA 3 actual sequence on + strand from BLAST (CDK6 is on - strand)
	             // but I need the sequence as would be read from - strand
	             // 5 TAATATGTATCAGCTAGCAATAT 3 
	  b.addWholeSite("taatatgtatcagctagcaatat");
	                      b.addSite("agcaata",15);
	             b.addSite("gtat",6);
#if 0	
	             // 5 UAAUAUGUAUCAGCUAGCAAUAU 3 mRNA sequence from the minus strand reading 5' to 3'
	             // 3 ATTATACATAGTCGATCGTTATA 5 DNA which would code this
	             // 5 ATATTGCTAGCTGATACATATTA 3 same DNA sequence correctly oriented
	  b.addWholeSite("atattgctagctgatacatatta");
	             // 5                AGCAAUA  3 first site on mRNA
	             // 3                TCGTTAT  5 DNA which would code this
	             // 5  TATTGCT                3 same DNA sequence correctly oriented, at 1
	        b.addSite("tattgct",1);
	             // 5       GUAU              3  mRNA sequence on the minus strand reading 5' to 3'
	             // 3       CATA              5  DNA which would code this
	             // 5              ATAC       3 same DNA sequence correctly oriented, at 13
	                    b.addSite("atac",13);

	             // UAAUAUGUAUCAGCUAGCAAUAU mRNA sequence from the minus strand reading 5' to 3'
	             // ATTATACATAGTCGATCGTTATA DNA which would code this
	b.addWholeSite("attatacatagtcgatcgttata");
	             //                AGCAAUA  first site on mRNA
	             //                TCGTTAT  DNA which would code this, at 15
	                    b.addSite("tcgttat",15);
	             //       GUAU              mRNA sequence on the minus strand reading 5' to 3'
	             //       CATA              DNA which would code this, at 6
	           b.addSite("cata",6);
	             // UAAUAUGUAUCAGCUAGCAAUAU mRNA sequence on the minus strand reading 5' to 3'
	             // UAUAACGAUCGACUAUGUAUAAU mRNA sequence as it would look opposite minus strand
	             // ATATTGCTAGCTGATACATATTA DNA sequence on - strand
	             // TATAACGATCGACTATGTATAAT DNA sequence on + strand, I think
	b.addWholeSite("tataacgatcgactatgtataat");
	             //                AGCAAUA  first site on mRNA
	             //  AUAACGA                mRNA site the opposite - strand
	             //  TATTGCT                DNA sequence of - strand
	             //  ATAAGCA                DNA on + strand at position 1
	      b.addSite("ataagca",1);
	             //       GUAU              mRNA sequence on the minus strand reading 5' to 3'
	             //              UAUG       mRNA sequence as it would look opposite minus strand
	             //              ATAC       DNA sequence on - strand
	             //              TATG       DNA on + strand at position 13
	                  b.addSite("tatg",13);

	b.addWholeSite("taatatgtatcagctagcaatat");
	                    b.addSite("agcaata",15);
	           b.addSite("gtat",6);

	b.addWholeSite("attatacatagtcgatcgttata");
	                    b.addSite("tcgttat",15);
	           b.addSite("cata",6);
#endif
}

char threePrimeUTR[60000],result[1000];
int getBindingSiteCoords(char *s,refseqGeneInfo &r,char *referencePath,bindingSite &b)
{
	char fn[100],rest[1000],siteSequence[MAXWHOLESITELENGTH],tempRes[1000];
	int st,en,l,pos,seqSt;
	faSequenceFile f;
	sprintf(fn,"%s\\%s.fa",referencePath,r.getChr());
	f.init(fn);
	r.getThreePrime(&st,&en,0);
	if (en-st>50000)
		if (r.getStrand()=='+')
			st=en-50000;
		else
			en=st+50000;
	f.getSequence(threePrimeUTR,st,en-st);
//	puts(threePrimeUTR);
	strlwr(threePrimeUTR);
	b.getSite(result,threePrimeUTR,r.getStrand()=='-');
	*s='\0';
	rest[0]='\0';
	while(sscanf(result,"%s %d %[^\n]",siteSequence,&pos,rest)>=2)
	{
		sprintf(tempRes,"%s\t%s:%d-%d\t%s\n",r.getGene(),r.getChr()+3,st+pos,st+pos+strlen(siteSequence),siteSequence);
		strcpy(result,rest);
		rest[0]='\0';
		strcat(s,tempRes);
	}
}

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100],rest[1000],compSeq[1000];
	int i;
	FILE *fp;
	gvaParams gp;
	analysisSpecs spec;

	dcerror.warn();

	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
	
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setBaitMargin(gp.margin);

	bindingSite b;
	lookFormiRNA137(b);
	r.goToStart();
	while (r.getNextGene())
	{
		getBindingSiteCoords(line,r,gp.referencePath,b);
		printf(line);
	}
}