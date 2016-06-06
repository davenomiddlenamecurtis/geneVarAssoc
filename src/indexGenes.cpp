#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef MSDOS 
#include <unistd.h>
#endif

int main(int argc,char *argv[])
{
	int i;
	FILE *fi,*fo;
	char line[1001],geneName[30],fn[30],*ptr;
	long linePos,locPos;
	unlink("refseqgeneindex.txt");
	fo=fopen("refseqgeneindex.txt","w");
	for (i=1;i<=3;++i)
	{
		sprintf(fn,"refseqgene%d.genomic.gbff",i);
		fi=fopen(fn,"rb");
		while(linePos=ftell(fi),fgets(line,1000,fi))
		{
			if (!strncmp(line,"LOCUS",5))
				locPos=linePos;
			else if ((ptr=strstr(line,"/gene="))!=0)
			{
				sscanf(line," /gene=\"%[^\"]",geneName);
				fprintf(fo,"%s\t%s\t%ld\n",geneName,fn,locPos);
			}

		}
		fclose(fi);
	}
	fclose(fo);

	return 0;
}