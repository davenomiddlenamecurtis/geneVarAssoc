#include <stdlib.h>
#include "geneVarUtils.hpp"

// assume that for bed file starts are 0 based and ends 1 based
// same applies to exons, etc.
// actual pos:
// 123456789
// consider position 7
// if it is a start would be called 6
// if it is an end will be called 7
// so if transcript start is 7 it will be called 6 in the code
// three bases before it are 456
// in the bed file we would output these as chr 3 6
// if exon end is 7 it will be called 7
// 2 bases after are 89
// in the bed file we output these as chr 7 9
// I hope that is clear

void refseqGeneInfo::getAllExons()
{
	int eStart,lowestPos,lowestT,lowestE,eCount,foundOne;
	eCount=0;
	eStart=firstExonStart;
	do {
		foundOne=0;
		lowestPos=999999999;
		for (int t=0;t<nTranscript;++t)
		{
			for (int e=0;e<transcript[t].exonCount;++e)
			{
				if (transcript[t].exonStarts[e]>=eStart)
				{
					foundOne=1;
					if (transcript[t].exonStarts[e]<lowestPos)
					{
						lowestPos=transcript[t].exonStarts[e];
						lowestT=t;
						lowestE=e;
					}
					break;
				}
			}
		}
		if (foundOne)
		{
			exonStarts[eCount]=transcript[lowestT].exonStarts[lowestE];
			eStart=exonEnds[eCount]=transcript[lowestT].exonEnds[lowestE];
			++eCount;
		}
	} while (foundOne);
	allExonCount=eCount;
	gotAllExons=1;
}

int refseqGeneInfo::writeBedFile(FILE *fo,writeBedFilePar &wbfp)
{
	int usedForKozak[MAXTRANSCRIPTPERGENE],nUsedForKozak,t,k,e,s;
	int usedForPromoter[MAXTRANSCRIPTPERGENE],nUsedForPromoter,p;
	if (gotAllExons==0)
		getAllExons();
	nUsedForKozak=0;
	if (strand=='+')
	{
	if (wbfp.write5UTR)
	{
		fprintf(fo,"%s\t%d\t%d\t%s_5PRIMEUTR\n",chr,transcript[0].txStart,transcript[0].cdsStart,geneName); // I think txStart and cdStart same for all transcripts
	}
	if (wbfp.write3UTR)
	{
		fprintf(fo,"%s\t%d\t%d\t%s_3PRIMEUTR\n",chr,transcript[0].cdsEnd,transcript[0].txEnd,geneName);
	}
	if (wbfp.writeKozaks)
	{
		for (t=0;t<nTranscript;++t)
		{
			for (k=0;k<nUsedForKozak;++k)
				if (transcript[t].cdsStart==usedForKozak[k])
					break;
			if (k==nUsedForKozak) // not used already
			{
				int s=transcript[t].cdsStart;
				usedForKozak[nUsedForKozak++]=s;
				// (gcc)gccRccAUGG - A is at s, 0-coded
				fprintf(fo,"%s\t%d\t%d\t%s_KOZAK_%d\n",chr,s-9,s+5,geneName,t+1);
			}
		}
	}
	nUsedForPromoter=0;
	if (wbfp.writePromoters)
	{
		for (t=0;t<nTranscript;++t)
		{
			for (p=0;p<nUsedForPromoter;++p)
				if (transcript[t].exonStarts[0]==usedForPromoter[p])
					break;
			if (p==nUsedForPromoter) // not used already
			{
				int s=transcript[t].exonStarts[0];
				usedForPromoter[nUsedForPromoter++]=s;
				// (gcc)gccRccAUGG - A is at s, 0-coded
				fprintf(fo,"%s\t%d\t%d\t%s_PROMOTER_%d\n",chr,s-upstream,s,geneName,t+1);
			}
		}
	}
	if (wbfp.writeSpliceSites)
	{
		for (e=0;e<allExonCount-1;++e)
		{
			s=exonEnds[e];
			fprintf(fo,"%s\t%d\t%d\t%s_SPLICE%d_DONOR\n",chr,s-NSBDONOREXON,s+NSBDONORINTRON,geneName,e+1);
			s=exonStarts[e+1];
			fprintf(fo,"%s\t%d\t%d\t%s_SPLICE%d_ACCEPTOR\n",chr,s-NSBACCEPTORINTRON,s+NSBACCEPTOREXON,geneName,e+2);
		}
	}
	if (wbfp.writeExons)
	{
		for (e=0;e<allExonCount;++e)
		{
			fprintf(fo,"%s\t%d\t%d\t%s_EXON%d\n",chr,exonStarts[e],exonEnds[e],geneName,e+1);
		}
	}
	}
	else // - strand
	{
	if (wbfp.write5UTR)
	{
		fprintf(fo,"%s\t%d\t%d\t%s_5PRIMEUTR\n",chr,transcript[0].cdsEnd,transcript[0].txEnd,geneName);
	}
	if (wbfp.write3UTR)
	{
		fprintf(fo,"%s\t%d\t%d\t%s_3PRIMEUTR\n",chr,transcript[0].txStart,transcript[0].cdsStart,geneName);
	}
	if (wbfp.writeKozaks)
	{
		for (t=0;t<nTranscript;++t)
		{
			for (k=0;k<nUsedForKozak;++k)
				if (transcript[t].cdsEnd==usedForKozak[k])
					break;
			if (k==nUsedForKozak) // not used already
			{
				int s=transcript[t].cdsEnd;
				usedForKozak[nUsedForKozak++]=s;
				// (gcc)gccRccAUGG - A is at s, 0-coded
				fprintf(fo,"%s\t%d\t%d\t%s_KOZAK_%d\n",chr,s-5,s+9,geneName,t+1);
			}
		}
	}
	nUsedForPromoter=0;
	if (wbfp.writePromoters)
	{
		for (t=0;t<nTranscript;++t)
		{
			for (p=0;p<nUsedForPromoter;++p)
				if (transcript[t].exonEnds[transcript[t].exonCount-1]==usedForPromoter[p])
					break;
			if (p==nUsedForPromoter) // not used already
			{
				int s=transcript[t].exonEnds[transcript[t].exonCount-1];
				usedForPromoter[nUsedForPromoter++]=s;
				// (gcc)gccRccAUGG - A is at s, 0-coded
				fprintf(fo,"%s\t%d\t%d\t%s_PROMOTER_%d\n",chr,s+1,s+1+upstream,geneName,t+1);
			}
		}
	}
	if (wbfp.writeSpliceSites)
	{
		for (e=allExonCount-2;e>=0;--e)
		{
			s=exonStarts[e+1];
			fprintf(fo,"%s\t%d\t%d\t%s_SPLICE%d_DONOR\n",chr,s-NSBDONORINTRON,s+NSBDONOREXON,geneName,allExonCount-e-1);
			s=exonEnds[e];
			fprintf(fo,"%s\t%d\t%d\t%s_SPLICE%d_ACCEPTOR\n",chr,s-NSBACCEPTOREXON,s+NSBACCEPTORINTRON,geneName,allExonCount-e);
		}
	}
	if (wbfp.writeExons)
	{
		for (e=allExonCount-1;e>=0;--e)
		{
			fprintf(fo,"%s\t%d\t%d\t%s_EXON%d\n",chr,exonStarts[e],exonEnds[e],geneName,allExonCount-e);
		}
	}
	}
	return 1;
}

int main(int argc,char *argv[])
{
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100],*fno;
	int i;
	FILE *fp,*fo,*feo;
	gvaParams gp;
	writeBedFilePar exonsPar,noExonsPar;
	analysisSpecs spec;
	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
	r.setListFile(gp.geneListFn);
	r.setBaitsFile(gp.baitFn);
	if (gp.referencePath[0]!='\0')
		r.setReferencePath(gp.referencePath);
	gp.upstream=PROMOTERLENGTH;
	gp.downstream=DOWNSTREAMLENGTH;
	r.setUpstream(gp.upstream);
	r.setDownstream(gp.downstream);
	r.setBaitMargin(gp.margin);
	fp=fopen(argv[2],"r");
	noExonsPar.writeKozaks=noExonsPar.writePromoters=noExonsPar.writeSpliceSites=noExonsPar.write5UTR=noExonsPar.write3UTR=1;
	exonsPar.writeExons=1;
	dcerror.warn();
	while (fgets(line,999,fp) && sscanf(line,"%s",geneName)==1)
	{
		if (!r.findGene(geneName) || !r.getNextGene())
		{
			FILE *fb=fopen("badgenes.txt","a");
			fprintf(fb,"%s\n",geneName);
			fclose(fb);
			continue;
		//	return 1;
		}
		sprintf(line,"%s.bed",geneName);
		fo=fopen(line,"w");
		r.writeBedFile(fo,noExonsPar);
		fclose(fo);
		sprintf(line,"%s.exons.bed",geneName);
		feo=fopen(line,"w");
		r.writeBedFile(feo,exonsPar);
		fclose(feo);
	}
	fclose(fp);
	return 0;
}