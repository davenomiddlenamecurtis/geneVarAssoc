#include <stdlib.h>
#include "geneVarUtils.hpp"
#include <assert.h>

#if 0
  hsa-miR-137  	   ABHD6   	  NM_020676  	  9  	  2164  	   UUAUUGCUU   	  2156  	  3 UTR  	  0.0038  	  1  
  hsa-miR-137  	   ABHD6   	  NM_020676  	  8  	  2163  	   UAUUGCUU   	  2156  	  3 UTR  	  0.0150  	  2  
  hsa-miR-137  	   ABTB1   	  NM_032548  	  7  	  1961  	   UUAUUGC   	  1955  	  3 UTR  	  0.0291  	  1  
  hsa-miR-137  	   ACSBG2   	  NM_030924  	  9  	  2442  	   UUAUUGCUU   	  2434  	  3 UTR  	  0.0018  	  1  
  hsa-miR-137  	   ACSBG2   	  NM_030924  	  7  	  2307  	   UUAUUGC   	  2301  	  3 UTR  	  0.0282  	  1  
#endif

#if 0
TargetScan 6.2:
PDLIM3	NM_001114107	PDZ and LIM domain 3	2	2	0	0	1	0	0	1	hsa-miR-137	-0.66	0.89	2007, 2009	Sites in UTR
APPL2	NM_018171	adaptor protein, phosphotyrosine interaction, PH domain and leucine zipper containing 2	1	1	0	0	1	1	0	0	hsa-miR-137	-0.63	0.71	2009	Sites in UTR
MITF	NM_000248	microphthalmia-associated transcription factor	3	2	0	1	1	0	1	0	hsa-miR-137	-0.60	0.96	2007, 2009	Sites in UTR
CEP128	NM_152446	centrosomal protein 128kDa	1	1	0	0	2	0	1	1	hsa-miR-137	-0.57	0.37	2009	Sites in UTR
RAVER2	NM_018211	ribonucleoprotein, PTB-binding 2	2	2	0	0	0	0	0	0	hsa-miR-137	-0.55	0.96	2007, 2009	Sites in UTR

#endif

char *sites[]= 
{
	// mir-137 "AGCAATAA","GCAATAA","AGCAATA",""
	"AGCTCCTA","GCTCCTA","AGCTCCT",""
};


  char *reverseStr(char *s)
  {
	  // reverse string in situ
	  char c;
	  int i,j;
	  for (i=0,j=strlen(s)-1;i<j;++i,--j)
	  {
		  c=s[i];
		  s[i]=s[j];
		  s[j]=c;
	  }
  return s;
  }

#define MAXINTPERFILE 1000
#define FIVEMARGIN 15
#define THREEMARGIN 1
int main(int argc,char *argv[])
{
	refseqTranscript t;
	char line[10000],geneName[100],transcriptName[100],target[100],compSeq[100],mir[100],*wholeSeq,*ptr,root[100],ext[100],fn[100];
	int i,start,end,st,en,len,lineCount,fileCount,found;
	FILE *fi,*fo,*fref,*fnf;
	faSequenceFile fa;

	assert((fref=fopen(argv[1],"r"))!=0);
	assert((fi=fopen(argv[2],"r"))!=0);
	strcpy(root,argv[3]);
	strcpy(ext,strchr(root,'.'));
	*strchr(root,'.')='\0';
	sprintf(fn,"%s.notfound",root);
	fnf=fopen(fn,"w");

//	for (i=0;sites[i][0]!='\0';++i)
//		RNA2DNA(sites[i]);
	lineCount=MAXINTPERFILE;
	fileCount=0;
	fo=NULL;
	while (fgets(line,999,fi))
	{
		found=0;
		if (sscanf(line,"%s %s",geneName,transcriptName)!=2 || strncmp(transcriptName,"NM",2))
			continue;
		strcat(transcriptName,"\t");
		fseek(fref,0L,SEEK_SET);
		while(fgets(line,9999,fref))
		{
			if (strstr(line,transcriptName))
			{
				found=1;
				intervalList matches; // constructor/destructor called each loop
				t.read(line);
				
				sprintf(line,"c:\\reference\\%s.fa",t.getChr());
				fa.init(line);

				for (i=0;sites[i][0]!='\0';++i)
					t.get3PrimeMatches(fa,matches,sites[i]);
				if (matches.nInts==0)
				{
					fprintf(fnf,"%s\t%s - no binding site found\t",geneName,transcriptName);
					if (t.getStrand()=='+')
						t.print3PrimeSequence(fa,fnf);
					fprintf(fnf,"\n");
				}
				matches.sort();
				matches.merge();
				for (i=0;i<matches.nInts;++i)
				{
		if (lineCount==MAXINTPERFILE)
		{
			if (fo!=NULL)
				fclose(fo);
			sprintf(fn,"%s%d%s",root,++fileCount,ext);
			assert((fo=fopen(fn,"w"))!=NULL);
			lineCount=0;
		}
		++lineCount;
		fprintf(fo,"%s\t%d\t%d\t",t.getChr(),matches.ints[i].st-(t.getStrand()=='+'?FIVEMARGIN:THREEMARGIN),matches.ints[i].en+(t.getStrand()=='+'?THREEMARGIN:FIVEMARGIN));
					fa.getSequence(line,matches.ints[i].st-(t.getStrand()=='+'?FIVEMARGIN:THREEMARGIN),matches.ints[i].en-matches.ints[i].st+1+(t.getStrand()=='+'?THREEMARGIN:FIVEMARGIN));
					fprintf(fo,"%s\t",line);
					fprintf(fo,"%s\t%s\t\n",geneName,transcriptName);
				}
			}
		}
		if (found==0)
		{
			fprintf(fnf,"%s\t%s not found in reference file\n",geneName,transcriptName);
		}
	}
	fclose(fnf);
	fclose(fi);
	fclose(fo);
	fclose(fref);
	return 0;
}

