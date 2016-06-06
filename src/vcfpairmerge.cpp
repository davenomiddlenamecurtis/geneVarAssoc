#include "masterLocusFile.hpp"
#include <stdlib.h>
#ifndef MSDOS
#include <unistd.h>
#endif


#define CHUNKSIZE 100000000

int main(int argc,char *argv[])
{
char *fn1,*fn2,*fn3,line[200],chrStr[10];
int chr,first,st,readSome;
FILE *fo;
masterLocusFile vf(2);
analysisSpecs spec; // just a dummy
fn1=argv[1];
fn2=argv[2];
fn3=argv[3];
fo=fopen(fn3,"wb");
first=1;
for (chr=1;chr<=23;++chr)
{
	st=0;
	if (chr<23)
		sprintf(chrStr,"%d",chr);
	else
		sprintf(chrStr,"X");
	do {
		unlink("tvcf.db");
		unlink("tvcf.vdx");
		unlink("merge1.vcf");
		unlink("merge2.vcf");
		vf.openFiles("tvcf.db","tvcf.vdx");
		sprintf(line,"tabix -h %s %s:%d-%d > merge1.vcf",fn1,chrStr,st,st+CHUNKSIZE-1);
		system(line);
		sprintf(line,"tabix -h %s %s:%d-%d > merge2.vcf",fn2,chrStr,st,st+CHUNKSIZE-1);
		system(line);
		readSome=0;
		vf.addLocusFile("merge1.vcf",VCFFILE);
		if (vf.readLocusFileEntries("merge1.vcf",spec,1))
			readSome=1;
		vf.addLocusFile("merged2.vcf",VCFFILE);
		if (vf.readLocusFileEntries("merge2.vcf",spec,1))
			readSome=1;
		if (readSome)
		{
			vf.openLocusFiles();
			if (first)
			{
				first=0;
				vf.outputMergedVCFHeader(fo);
			}
			vf.outputMergedVcfGenotypes(fo,spec);
			vf.closeLocusFiles();
			st+=CHUNKSIZE;
		}
		vf.closeFiles();
	} while (readSome>0);
}
fclose(fo);
return 0;
}