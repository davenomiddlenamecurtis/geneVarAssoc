#include "masterLocusFile.hpp"
#include <stdlib.h>
#ifndef MSDOS
#include <unistd.h>
#define COPY "cp"
#else
#define COPY "copy"
#endif


#define CHUNKSIZE 100000000

#define unlink(x) system("echo deleted >" x); _unlink(x)

int main(int argc,char *argv[])
{
char *fn1,*fn2,*fn3,line[200],chrStr[10];
int chr,first,st,readSome,a;
FILE *fo;
masterLocusFile vf(2);
analysisSpecs spec; // just a dummy
if (argc<3 || strstr(argv[1],".vcf.gz")==0)
{
	printf("run as %s merged.vcf.gz ROOT*.vcf.gz\n",argv[0]);
	printf("if merged.vcf.gz already exists then new files will be added to it\n");
	printf("tabix.exe and bgzip.exe must be available on PATH");
	return 1;
}
fo=fopen(argv[1],"rb");
if (fo==NULL)
{
	sprintf(line,"%s %s %s",COPY,argv[2],argv[1]);
	system(line);
	sprintf(line,"tabix -p vcf %s",argv[1]);
	system(line);
	a=3;
}
else
{
	fclose(fo);
	sprintf(line,"%s.tbi",argv[1]);
	fo=fopen(line,"rb");
	if (fo==0)
	{
		sprintf(line,"tabix -p vcf %s",argv[1]);
		system(line);
	}
	else
		fclose(fo);
	a=2;
}
for (;a<argc;++a)
{
fn1=argv[1];
fn2=argv[a];
printf("Merging %s and %s...\n",fn1,fn2);
unlink("tempmerged.vcf");
fn3="tempmerged.vcf";
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
		vf.addLocusFile("merge2.vcf",VCFFILE);
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
unlink("tempmerged.vcf.gz");
system("bgzip tempmerged.vcf");
sprintf(line,"%s tempmerged.vcf.gz %s",COPY,argv[1]);
system(line);
sprintf(line,"tabix -p vcf %s",argv[1]);
system(line);
}
return 0;
}