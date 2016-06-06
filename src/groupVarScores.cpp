#include <stdlib.h>
#include "geneVarUtils.hpp"

#define MAXSUB 5000
#define MAXGENE 5000

int main(int argc,char *argv[])
{
	geneExtractor gcase,gcont;
	refseqGeneInfo r;
	char fn[100],fn2[100],line[1000],geneName[100],groupName[100];
	int i,cc[MAXSUB],s,g,nsub,ngene;
	float **score,totalScore[MAXSUB];
	FILE *fp,*fg,*fs;
	gvaParams gp;
	analysisSpecs spec;
	fp=fopen(argv[1],"r");
	gp.input(fp,spec);
	fclose(fp);
	masterLocusFile vf(gp.nCc[0]+gp.nCc[1]);
	strcpy(groupName,argv[2]);
	sprintf(fn,"%s.dat",groupName);
	fp=fopen(fn,"w");
	fg=fopen(groupName,"r");
	fprintf(fp,"ID     CC ");
	ngene=0;
	while (fgets(line,999,fg) && sscanf(line,"%s",geneName)==1)
	{
		geneName[8]='\0'; // truncate
		fprintf(fp,"%-8s ",geneName);
		++ngene;
	}
	fprintf(fp,"Total     \n");
	fseek(fg,0L,SEEK_SET);
	score=(float **)malloc(MAXSUB*sizeof(float*));
	for (s=0;s<MAXSUB;++s)
	{
		score[s]=(float*)malloc(ngene*sizeof(float));
		totalScore[s]=0;
	}
	
	for (g=0;fgets(line,999,fg) && sscanf(line,"%s",geneName)==1;++g)
	{
	if (spec.consequenceThreshold!=0)
		sprintf(fn,"gva.%s.ct%02d.sco",geneName,spec.consequenceThreshold);
	else if (spec.useConsequenceWeights!=0)
		sprintf(fn,"gva.%s.ucw.sco",geneName);
	else
		sprintf(fn,"gva.%s.sco",geneName);
	fs=fopen(fn,"r");
	if (fs==0)
	{
		sprintf(line,"geneVarAssoc %s %s",argv[1],geneName);
		system(line);
		fs=fopen(fn,"r");
	}
	s=0;
	while (fgets(line,999,fs) && sscanf(line,"%*s %d  %f",&cc[s],&score[s][g])==2)
	{
		totalScore[s]+=score[s][g];
		++s;
	}
	fclose(fs);
	}
	fclose(fg);
	nsub=s;
	for (s=0;s<nsub;++s)
	{
		fprintf(fp,"%-6d %d  ",s+1,cc[s]);
		for (g=0;g<ngene;++g)
			fprintf(fp,"%8.4f ",score[s][g]);
		fprintf(fp,"%8.4f\n",totalScore[s]);
	}
	fclose(fp);
	for (s=0;s<MAXSUB;++s)
		free(score[s]);
	free(score);
	sprintf(fn,"%s.inp",groupName);
	fp=fopen(fn,"w");
	fprintf(fp,"d %s.dat\n",groupName);
	fprintf(fp,"o %s.out\n",groupName);
	fprintf(fp,"tt\nTotal\ncc=0\ncc=1\n");
	fprintf(fp,"co\n%d\n",ngene);
	for (g=0;g<ngene;++g)
		fprintf(fp,"C%d\n",g+3);
	fprintf(fp,"0.1\n");
	fprintf(fp,"wid\nm\nCC\n%d\n",ngene);
	for (g=0;g<ngene;++g)
		fprintf(fp,"C%d\n",g+3);
	fprintf(fp,"o close\n");
	fclose(fp);
	return 0;
}