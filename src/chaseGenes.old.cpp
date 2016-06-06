#include <windows.h>
#include <Wininet.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define REFDIR "c:\\reference\\"
#define INDEX  "c:\\reference\\refseqgeneindex.txt"

// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(brca1%5BGene%20Name%5D)%20AND%20homo%20sapiens%5BOrganism%5D
// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=12190&retmode=xml


char comment_buff[100000];
char commentStr_buff[100000];
char pubMedComment_buff[100000];
char IDList_buff[100000];
char oneComment_buff[100000];
char allComments_buff[100000];

HINTERNET hInternet;

char *getPubMedComments(char *geneName,char *pubMedComment,char *unused,char *IDList)
{
	char query[1000],id[20],*ptr,*sptr,*eptr,*IDptr,*cptr,*epmptr,*oneComment;
	int nComments;
	int rv,i;
	HINTERNET eutils;
	LPCTSTR lpszAgent,lpszServerName,queryPtr;
	DWORD dwContext,nBytesRead;
	DWORD dwRequestFlags = INTERNET_FLAG_NO_UI   // no UI please
            |INTERNET_FLAG_NO_AUTH           // don't authenticate
            |INTERNET_FLAG_PRAGMA_NOCACHE    // do not try the cache or proxy
            |INTERNET_FLAG_NO_CACHE_WRITE   // don't add this to the IE cache
			| INTERNET_FLAG_RELOAD; // desperate attempt to fix read problem
	dwContext=1;
	lpszAgent=LPCTSTR("iexplorer.exe");
	lpszServerName=LPCTSTR("http://eutils.ncbi.nlm.nih.gov");
//	hInternet=InternetOpen("iexplorer.exe",INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL, 0 );
//	hInternet=InternetOpen("iexplorer.exe",INTERNET_OPEN_TYPE_DIRECT, NULL, NULL, 0 );
	*pubMedComment='\0';
	nComments=0;
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=schizophrenia%%20AND%%20%s",
	geneName);
	queryPtr=LPCTSTR(query);
	for (i=0;i<10;++i)
	{
		eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
		if (eutils!=NULL)
			break;
	}
	if (i==11)
		return("This failed: eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);"); 
	ptr=IDList;
	*ptr='\0';
	for (i=0;i<100;++i)
	{
		Sleep(100);
		rv=InternetReadFile(eutils,ptr,3999,(LPDWORD)&nBytesRead);
		if (rv==TRUE && nBytesRead==0)
			break;
		ptr=(strchr(ptr,'\0'));
	}
	if (i==101)
		return("This failed: rv=InternetReadFile(eutils,ptr,39999,(LPDWORD)&nBytesRead);"); 
	InternetCloseHandle(eutils);
printf("%s\n",IDList);
	IDptr=IDList;
	oneComment=0;
	while ((ptr=strstr(IDptr,"<Id>"))!=0)
	{
		if (oneComment)
			free(oneComment);
		oneComment=(char*)malloc(40000);
		assert(oneComment);
		printf("%s\n",pubMedComment);
		IDptr=ptr+1;
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
		continue;
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=xml",
	id);
	for (i=0;i<10;++i)
	{
		eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
		if (eutils!=NULL)
			break;
	}
	if (i==11)
		return("This failed: eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);"); 
	cptr=oneComment;
	*cptr='\0';
	for (i=0;i<100;++i)
	{
		Sleep(100);
	rv=InternetReadFile(eutils,cptr,3999,(LPDWORD)&nBytesRead);
	if (rv==TRUE && nBytesRead==0)
		break;
	cptr=strchr(cptr,'\0');
	}
	if (i==101)
		return("This failed: InternetReadFile(eutils,oneComment,39999,(LPDWORD)&nBytesRead);"); 
	InternetCloseHandle(eutils);
	sptr=strstr(oneComment,"Name=\"Title\"");
	if (sptr==0)
		sptr=strstr(oneComment,"<title>");
	if (sptr==0 || (sptr=strchr(sptr,'>'))==0 || (eptr=strchr(sptr,'<'))==0)
		continue;
	*eptr='\0';
	epmptr=strchr(pubMedComment,'\0');
	sprintf(epmptr,"ITEM: %s http://www.ncbi.nlm.nih.gov/pubmed?term=%s ",sptr+1,id);
	if (++nComments>=5)
		break;
	}
		if (oneComment)
		{
			free(oneComment);
			oneComment=0;
		}
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%s",
	geneName);
	queryPtr=LPCTSTR(query);
	for (i=0;i<10;++i)
	{
		eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
		if (eutils!=NULL)
			break;
	}
	if (i==11)
		return("This failed: eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);"); 
	ptr=IDList;
	*ptr='\0';
	for (i=0;i<100;++i)
	{Sleep(100);
		rv=InternetReadFile(eutils,ptr,3999,(LPDWORD)&nBytesRead);
		if (rv==TRUE && nBytesRead==0)
			break;
		ptr=(strchr(ptr,'\0'));
	}
	if (i==101)
		return("This failed: rv=InternetReadFile(eutils,ptr,39999,(LPDWORD)&nBytesRead);"); 
	InternetCloseHandle(eutils);
printf("%s\n",IDList);
	IDptr=IDList;
	while ((ptr=strstr(IDptr,"<Id>"))!=0)
	{
		if (oneComment)
			free(oneComment);
		oneComment=(char*)malloc(40000);
		assert(oneComment);
		printf("%s\n",pubMedComment);
	IDptr=ptr+1;
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
		continue;
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=xml",
	id);
	queryPtr=LPCTSTR(query);
	for (i=0;i<10;++i)
	{
		eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
		if (eutils!=NULL)
			break;
	}
	if (i==11)
		return("This failed: eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);"); 
	cptr=oneComment;
	*cptr='\0';
	for (i=0;i<100;++i)
	{Sleep(100);
	rv=InternetReadFile(eutils,cptr,3999,(LPDWORD)&nBytesRead);
	if (rv==TRUE && nBytesRead==0)
		break;
	cptr=strchr(cptr,'\0');
	}
	if (i==101)
		return("This failed: InternetReadFile(eutils,oneComment,39999,(LPDWORD)&nBytesRead);"); 
	InternetCloseHandle(eutils);
//	printf("%s\n",oneComment);
	if ((sptr=strstr(oneComment,"Name=\"Title\""))==0 || (sptr=strchr(sptr,'>'))==0 || (eptr=strchr(sptr,'<'))==0)
		continue;
	*eptr='\0';
	epmptr=strchr(pubMedComment,'\0');
	sprintf(epmptr,"ITEM: %s http://www.ncbi.nlm.nih.gov/pubmed?term=%s ",sptr+1,id);
	if (++nComments>=5)
		break;
	}
		if (oneComment)
			free(oneComment);
	return pubMedComment;
}

char *getGeneComment(char *geneName,char *unused,char *commentStr)
{
	char query[1000],id[20],*ptr,*eptr,*comment;
	int rv,i;
	HINTERNET eutils;
	LPCTSTR lpszAgent,lpszServerName,queryPtr;
	DWORD dwContext,nBytesRead;
	DWORD dwRequestFlags = INTERNET_FLAG_NO_UI   // no UI please
            |INTERNET_FLAG_NO_AUTH           // don't authenticate
            |INTERNET_FLAG_PRAGMA_NOCACHE    // do not try the cache or proxy
            |INTERNET_FLAG_NO_CACHE_WRITE   // don't add this to the IE cache
			| INTERNET_FLAG_RELOAD; // desperate attempt to fix read problem
	dwContext=1;
	lpszAgent=LPCTSTR("iexplorer.exe");
	lpszServerName=LPCTSTR("http://eutils.ncbi.nlm.nih.gov");
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(%s[Gene Name] AND Homo sapiens[Organism])",
//	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(%s%%5BGene%%20Name%%5D)%%20AND%%20homo%%20sapiens%%5BOrganism%%5D",
	geneName);
	queryPtr=LPCTSTR(query);
	for (i=0;i<10;++i)
	{
		if ((eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0))!=NULL)
			break;
	}
	if (i==11)
		return("Failed: InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0)");
	comment=(char*)malloc(40000);
	assert(comment);
	ptr=comment;
	*ptr='\0';
	for (i=0;i<100;++i)
	{Sleep(100);
		rv=InternetReadFile(eutils,ptr,3999,(LPDWORD)&nBytesRead);
		if (rv==TRUE && nBytesRead==0)
			break;
		ptr=strchr(ptr,'\0');
	}
	if (i==101)
		return ("Failed: rv=InternetReadFile(eutils,ptr,39999,(LPDWORD)&nBytesRead);");
	InternetCloseHandle(eutils);
printf("%s\n",comment);
	if ((ptr=strstr(comment,"<Id>"))==0)
	{
		InternetCloseHandle(eutils);
		return "Could not find <Id>";
	}
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
	{
		InternetCloseHandle(eutils);
		return "Could not sscanf Id";
	}
			free(comment);
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s&retmode=xml",id);
	for (i=0;i<10;++i)
	{
		if ((eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0))!=NULL)
			break;
	}
	if (i==11)
		return "Failed: if ((eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0))!=NULL)";

	comment=(char*)malloc(40000);
	assert(comment);
	ptr=comment;
	*ptr='\0'; ptr[1]='\0';
	for (i=0;i<100;++i)
	{Sleep(100);
		rv=InternetReadFile(eutils,ptr,3999,(LPDWORD)&nBytesRead);
		if (rv==TRUE && nBytesRead==0)
			break;
		ptr=strchr(ptr,'\0');
	}
	if (i==101)
		return "Failed: rv=InternetReadFile(eutils,ptr,39999,(LPDWORD)&nBytesRead);";
	InternetCloseHandle(eutils);
// printf("%s\n",comment);
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
	free(comment);
	return commentStr;
}

void getMLP(char *fn,float *mlp,float *recMlp,float *HWEMlp,float *HWEContMlp,float *HOMMlp,float *HOMContMlp)

{
	FILE *fp;
	char buff[1000];
	*mlp=*recMlp=*HWEMlp=*HWEContMlp=*HOMMlp=*HOMContMlp=0;
	fp=fopen(fn,"r");
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
	return;
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
	sprintf(line,"\\msvc\\vcf\\geneVarAssoc \\sequence\\sharedseq\\gva.uk.ob.120705.par %s",geneName);
//	system(line);
	sprintf(line,"gva.%s.ucw.sao",geneName);
	getMLP(line,&allMlp,&recMlp,&HWEMlp,&HWEContMlp,&HOMMlp,&HOMContMlp);
	getFreqs(line,freqBuff);
	sprintf(line,"\\msvc\\vcf\\geneVarAssoc \\sequence\\sharedseq\\gva.uk.ob.ns.120705.par %s",geneName);
//	system(line);
	sprintf(line,"gva.%s.ct09.sao",geneName);
	getMLP(line,&nsMlp,&recMlp,&HWEMlp,&HWEContMlp,&HOMMlp,&HOMContMlp);
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
		gcCommStr=getGeneComment(geneName,comment_buff,commentStr_buff);
		pmCommStr=getPubMedComments(geneName,pubMedComment_buff,oneComment_buff,IDList_buff);
		fprintf(fo,"%s\t%s\t%s\t%s\n",gcCommStr,pmCommStr,freqBuff,genoBuff);
	return 1;
}

int main(int argc,char *argv[])
{
	char line[200],geneName[20];
	FILE *fi,*fo;
	hInternet=InternetOpen("agent",INTERNET_OPEN_TYPE_DIRECT, NULL, NULL, 0 );

	fi=fopen(argv[1],"r");
	fo=fopen(argv[2],"w");
	fprintf(fo,"Gene\tInteresting\tAll\tNS\tRec\tHWECase\tHWECont\tHOMCase\tHOMCont\tSummary\tPubMed\tFreqs\tHOMgenos\n");
	while (fgets(line,199,fi) && sscanf(line,"%s",geneName)==1)
	{
		chaseGene(fo,geneName);
	}
	fclose(fo);
	fclose(fi);
	InternetCloseHandle(hInternet);
	return 0;
}