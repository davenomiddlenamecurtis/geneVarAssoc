#include <windows.h>
#include <Wininet.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#ifndef MSDOS 
#include <unistd.h>
#endif

#define REFDIR "c:\\reference\\"
#define INDEX  "c:\\reference\\refseqgenesorted.txt"
#if 1
#define ALLVARFILE "\\users\\Dave\\skydrive\\sharedseq\\gva.uk.ob.all.121127.par"
#define NSVARFILE "\\users\\Dave\\skydrive\\sharedseq\\gva.uk.ob.ct08.rare.121127.par"
#else
#define ALLVARFILE "\\sequence\\sharedseq\\gva.fb.ob.120705.par"
#define NSVARFILE "\\sequence\\sharedseq\\gva.fb.ob.ns.120705.par"
#endif

// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(brca1%5BGene%20Name%5D)%20AND%20homo%20sapiens%5BOrganism%5D
// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=12190&retmode=xml


char comment[40000];
char commentStr[40000];
char pubMedComment[40000];
char IDList[40000];
char oneComment[40000];
char allComments[40000];

char *getPubMedComments(char *geneName)
{
	char query[1000],id[20],*ptr,*sptr,*eptr,*IDptr;
	FILE *fp;
	int nComments,c;
	*pubMedComment='\0';
	nComments=0;
	unlink("ids.txt");
	sprintf(query,
	"\\wget\\bin\\wget \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=schizophrenia%%20AND%%20%s\" -O ids.txt",
	geneName);
	system(query);
	fp=fopen("ids.txt","r");
	ptr=IDList;
	while ((c=fgetc(fp))!=EOF)
		*ptr++=c;
	*ptr='\0';
	fclose(fp);

printf("%s\n",IDList);
	IDptr=IDList;
	while ((ptr=strstr(IDptr,"<Id>"))!=0)
	{
	IDptr=ptr+1;
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
		continue;
	unlink("comment.txt");
	sprintf(query,
	"\\wget\\bin\\wget \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=xml\" -O comment.txt",
	id);
	system(query);
	fp=fopen("comment.txt","r");
	ptr=oneComment;
	while ((c=fgetc(fp))!=EOF)
		*ptr++=c;
	*ptr='\0';
	fclose(fp);

printf("%s\n",oneComment);
	if ((sptr=strstr(oneComment,"Name=\"Title\""))==0 || (sptr=strchr(sptr,'>'))==0 || (eptr=strchr(sptr,'<'))==0)
		continue;
	*eptr='\0';
	eptr=strchr(pubMedComment,'\0');
	sprintf(eptr,"ITEM: %s http://www.ncbi.nlm.nih.gov/pubmed?term=%s ",sptr+1,id);
	if (++nComments>=5)
		break;
	}
	unlink("ids.txt");
	sprintf(query,
	"\\wget\\bin\\wget \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%s\" -O ids.txt",
	geneName);
	system(query);
	fp=fopen("ids.txt","r");
	ptr=IDList;
	while ((c=fgetc(fp))!=EOF)
		*ptr++=c;
	*ptr='\0';
	fclose(fp);

printf("%s\n",IDList);
	IDptr=IDList;
	while ((ptr=strstr(IDptr,"<Id>"))!=0)
	{
	IDptr=ptr+1;
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
		continue;
	unlink("comment.txt");
	sprintf(query,
	"\\wget\\bin\\wget \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=xml\" -O comment.txt",
	id);
	system(query);
	fp=fopen("comment.txt","r");
	ptr=oneComment;
	while ((c=fgetc(fp))!=EOF)
		*ptr++=c;
	*ptr='\0';
	fclose(fp);
printf("%s\n",oneComment);
	if ((sptr=strstr(oneComment,"Name=\"Title\""))==0 || (sptr=strchr(sptr,'>'))==0 || (eptr=strchr(sptr,'<'))==0)
		continue;
	*eptr='\0';
	eptr=strchr(pubMedComment,'\0');
	sprintf(eptr,"ITEM: %s http://www.ncbi.nlm.nih.gov/pubmed?term=%s ",sptr+1,id);
	if (++nComments>=5)
		break;
	}
	return pubMedComment;
}

char *getGeneComment(char *geneName)
{
	char query[1000],id[20],*ptr,*eptr;
	int c;
	FILE *fp;
	unlink("comment.txt");
	sprintf(query,
	"\\wget\\bin\\wget \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(%s%%5BGene%%20Name%%5D)%%20AND%%20homo%%20sapiens%%5BOrganism%%5D\" -O comment.txt",
	geneName);
	system(query);
	fp=fopen("comment.txt","r");
	ptr=comment;
	while ((c=fgetc(fp))!=EOF)
		*ptr++=c;
	*ptr='\0';
	fclose(fp);
printf("%s\n",comment);
	if ((ptr=strstr(comment,"<Id>"))==0)
	{
		return "";
	}
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
	unlink("comment.txt");
	sprintf(query,
	"\\wget\\bin\\wget \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s&retmode=xml\" -O comment.txt",
	id);
	system(query);
	fp=fopen("comment.txt","r");
	ptr=comment;
	while ((c=fgetc(fp))!=EOF)
		*ptr++=c;
	*ptr='\0';
	fclose(fp);
	printf("%s\n",comment);
	*commentStr='\0';
	if ((ptr=strstr(comment,"Name=\"NomenclatureName\""))!=0 && (ptr=strchr(ptr,'>'))!=0 && strchr(ptr,'<')!=0)
	{
		sscanf(ptr+1,"%[^<]",commentStr);
	}
	if ((ptr=strstr(comment,"Name=\"Summary\""))!=0 && (ptr=strchr(ptr,'>'))!=0 && (eptr=strchr(ptr,'<'))!=0)
	{
		*eptr='\0';
		sprintf(strchr(commentStr,'\0')," : %s",ptr+1);
	}
	return commentStr;
}

int getMLP(char *fn,float *mlp,float *recMlp,float *HWEMlp,float *HWEContMlp,float *HOMMlp,float *HOMContMlp)
{
	FILE *fp;
	char buff[1000];
	*mlp=*recMlp=*HWEMlp=*HWEContMlp=*HOMMlp=*HOMContMlp=0;
	fp=fopen(fn,"r");
	if (fp==0)
		return 0;
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",mlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",recMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HWEMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HWEContMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HOMMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HOMContMlp);
	fclose(fp);
	return 1;
}

int getTestRecMLP(char *fn,float *mlp,float *recMlp,float *HWEMlp,float *HWEContMlp,float *HOMMlp,float *HOMContMlp)
{
	FILE *fp;
	char buff[1000],ftrn[100];
	strcpy(ftrn,fn);
	strcat(strstr(ftrn,".sao"),".tro");
	unlink(ftrn);
	sprintf(buff,"\\msvc\\gc\\testRecHWE2 %s %s",fn,ftrn);
	system(buff);
	*recMlp=*HWEMlp=*HWEContMlp=0;
	fp=fopen(ftrn,"r");
	if (fp==0)
		return 0;
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",recMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HWEMlp);
	while (*buff=0,fgets(buff,1999,fp) && strncmp("-log(p)",buff,strlen("-log(p)"))) ;
	sscanf(buff,"-log(p) =%f",HWEContMlp);
	fclose(fp);
	return 1;
}

void getFreqs(char *fn,char *freqBuff)
{
	FILE *fp;
	char buff[1000],word[30],freq0[30],freq1[30],wt[30];
	float mlp;
	mlp=0;
	fp=fopen(fn,"r");
	while (fgets(buff,999,fp) && sscanf(buff,"%s",word)==1 && strcmp(word,"AA"))
		;
	*freqBuff='\0';
	while (fgets(buff,999,fp) && sscanf(buff,"%*s %*s %*s %*s %s %*s %*s %*s %s %*s %*s %s",freq0,freq1,wt)==3)
		if (atof(freq0)!=0 || atof(freq1)!=0)
			sprintf(strchr(freqBuff,'\0'),"%s/%s:%s ",freq0,freq1,wt);
	fclose(fp);
}

void getGenos(char *fn,char *genoBuff)
{
	FILE *fp;
	char buff[1000],word[30],genos[1000];
	int cc;
	float mlp;
	mlp=0;
	fp=fopen(fn,"r");
	while (fgets(buff,999,fp) && sscanf(buff,"%s",word)==1 && strcmp(word,"-log(p)"))
		;
	*genoBuff='\0';
	while (fgets(buff,999,fp) && sscanf(buff,"%*s %d %s",&cc,genos)==2)
		sprintf(strchr(genoBuff,'\0'),"%d:%s ",cc,genos);
	fclose(fp);
}

#if 1
#define NCOHORT 0
char **cohort=0;
#else
#define NCOHORT 3
char *cohort[NCOHORT]= { "AB","ED","UK" };
#endif

char freqBuff[20000],genoBuff[20000];

int chaseGene(FILE *fo,char *geneName)

{
	char line[200],*gcCommStr,*pmCommStr;
	float mlp1,mlpThreshold,mlp10,mlp1000,mlps[10],mlpRAW,mlpSNP;
	float mlp,recMlp,HWEMlp,HWEContMlp,HOMMlp,HOMContMlp,allMlp,nsMlp;
	int a,i;
	sprintf(line,"gva.%s.ucw.sao",geneName);
	unlink(line);
	sprintf(line,"\\msvc\\vcf\\geneVarAssoc " ALLVARFILE " %s",geneName);
	system(line);
	sprintf(line,"gva.%s.ucw.sao",geneName);
	if (!getMLP(line,&allMlp,&recMlp,&HWEMlp,&HWEContMlp,&HOMMlp,&HOMContMlp))
	{
		fprintf(fo,"%s\n",geneName);
		return 0;
	}
	getFreqs(line,freqBuff);
	sprintf(line,"\\msvc\\vcf\\geneVarAssoc " NSVARFILE " %s",geneName);
	system(line);
	sprintf(line,"gva.%s.ct08.sao",geneName);
	getMLP(line,&nsMlp,&recMlp,&HWEMlp,&HWEContMlp,&HOMMlp,&HOMContMlp);
	// getTestRecMLP(line,&nsMlp,&recMlp,&HWEMlp,&HWEContMlp,&HOMMlp,&HOMContMlp);
	getGenos(line,genoBuff);
	sprintf(line,"del *%s*.db",geneName);
	system(line);
	sprintf(line,"del *%s*.vdx",geneName);
	system(line);
	sprintf(line,"del *%s*.dat",geneName);
	system(line);
	sprintf(line,"del *%s*.par",geneName);
	system(line);
		fprintf(fo,"%s\t\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t",
			geneName,allMlp,nsMlp,recMlp,HWEMlp,HWEContMlp,HOMMlp,HOMContMlp);
#if 1
		gcCommStr=getGeneComment(geneName);
		pmCommStr=getPubMedComments(geneName);
#else
		gcCommStr="";
		pmCommStr="";
#endif
		fprintf(fo,"%s\t%s\t%s\t%s\n",gcCommStr,pmCommStr,freqBuff,genoBuff);
	return 1;
}

int main(int argc,char *argv[])
{
	char line[200],geneName[20];
	FILE *fi,*fo;

	fi=fopen(argv[1],"r");
	fo=fopen(argv[2],"w");
	fprintf(fo,"Gene\tInteresting\tAll\tNS\tRec\tHWECase\tHWECont\tHOMCase\tHOMCont\tSummary\tPubMed\tFreqs\tHOMgenos\n");
	while (fgets(line,199,fi) && sscanf(line,"%s",geneName)==1)
	{
		chaseGene(fo,geneName);
	}
	fclose(fo);
	fclose(fi);
	return 0;
}