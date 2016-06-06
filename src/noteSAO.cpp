#include <windows.h>
#include <Wininet.h>
#include <stdio.h>
#include <stdlib.h>

#define REFDIR "c:\\reference\\"
#define INDEX  "c:\\reference\\refseqgeneindex.txt"

// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(brca1%5BGene%20Name%5D)%20AND%20homo%20sapiens%5BOrganism%5D
// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=12190&retmode=xml


char comment[40000];
char IDList[40000];
char oneComment[40000];
char allComments[40000];
char pmcComments[40000];

char *getPubMedComments(char *geneName)
{
	char query[1000],id[20],*ptr,*sptr,*eptr,*IDptr;
	int nComments;
	HINTERNET hInternet,eutils;
	LPCTSTR lpszAgent,lpszServerName,queryPtr;
	DWORD dwContext,nBytesRead;
	DWORD dwRequestFlags = INTERNET_FLAG_NO_UI   // no UI please
            |INTERNET_FLAG_NO_AUTH           // don't authenticate
            |INTERNET_FLAG_PRAGMA_NOCACHE    // do not try the cache or proxy
            |INTERNET_FLAG_NO_CACHE_WRITE;   // don't add this to the IE cache
	dwContext=1;
	lpszAgent=LPCTSTR("iexplorer.exe");
	lpszServerName=LPCTSTR("http://eutils.ncbi.nlm.nih.gov");
	hInternet=InternetOpen("iexplorer.exe",INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL, 0 );
	*pmcComments='\0';
	nComments=0;
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=schizophrenia%%20AND%%20%s",
	geneName);
	queryPtr=LPCTSTR(query);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	*IDList='\0';
	InternetReadFile(eutils,IDList,39999,(LPDWORD)&nBytesRead);
	InternetCloseHandle(eutils);
printf("%s\n",IDList);
	IDptr=IDList;
	while ((ptr=strstr(IDptr,"<Id>"))!=0)
	{
	IDptr=ptr+1;
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
		continue;
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=xml",
	id);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	if (InternetReadFile(eutils,oneComment,39999,(LPDWORD)&nBytesRead)!=TRUE)
	{
		InternetCloseHandle(eutils);
		continue;
	}
	InternetCloseHandle(eutils);
printf("%s\n",oneComment);
	if ((sptr=strstr(oneComment,"Name=\"Title\""))==0 || (sptr=strchr(sptr,'>'))==0 || (eptr=strchr(sptr,'<'))==0)
		continue;
	*eptr='\0';
	eptr=strchr(pmcComments,'\0');
	sprintf(eptr,"ITEM: %s http://www.ncbi.nlm.nih.gov/pubmed?term=%s ",sptr+1,id);
	if (++nComments>=5)
		break;
	}
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=%s",
	geneName);
	queryPtr=LPCTSTR(query);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	*IDList='\0';
	InternetReadFile(eutils,IDList,39999,(LPDWORD)&nBytesRead);
	InternetCloseHandle(eutils);
printf("%s\n",IDList);
	IDptr=IDList;
	while ((ptr=strstr(IDptr,"<Id>"))!=0)
	{
	IDptr=ptr+1;
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
		continue;
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=%s&retmode=xml",
	id);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	if (InternetReadFile(eutils,oneComment,39999,(LPDWORD)&nBytesRead)!=TRUE)
	{
		InternetCloseHandle(eutils);
		continue;
	}
	InternetCloseHandle(eutils);
printf("%s\n",oneComment);
	if ((sptr=strstr(oneComment,"Name=\"Title\""))==0 || (sptr=strchr(sptr,'>'))==0 || (eptr=strchr(sptr,'<'))==0)
		continue;
	*eptr='\0';
	eptr=strchr(pmcComments,'\0');
	sprintf(eptr,"ITEM: %s http://www.ncbi.nlm.nih.gov/pubmed?term=%s ",sptr+1,id);
	if (++nComments>=5)
		break;
	}
	InternetCloseHandle(hInternet);
	return pmcComments;
}

char *getGeneComment(char *geneName)
{
	char query[1000],id[20],*ptr,*eptr;
	HINTERNET hInternet,eutils;
	LPCTSTR lpszAgent,lpszServerName,queryPtr;
	DWORD dwContext,nBytesRead;
	DWORD dwRequestFlags = INTERNET_FLAG_NO_UI   // no UI please
            |INTERNET_FLAG_NO_AUTH           // don't authenticate
            |INTERNET_FLAG_PRAGMA_NOCACHE    // do not try the cache or proxy
            |INTERNET_FLAG_NO_CACHE_WRITE;   // don't add this to the IE cache
	dwContext=1;
	lpszAgent=LPCTSTR("iexplorer.exe");
	lpszServerName=LPCTSTR("http://eutils.ncbi.nlm.nih.gov");
	hInternet=InternetOpen("iexplorer.exe",INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL, 0 );
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(%s%%5BGene%%20Name%%5D)%%20AND%%20homo%%20sapiens%%5BOrganism%%5D",
	geneName);
	queryPtr=LPCTSTR(query);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	if (InternetReadFile(eutils,comment,39999,(LPDWORD)&nBytesRead)!=TRUE)
	{
		InternetCloseHandle(hInternet);
		InternetCloseHandle(eutils);
		return "";
	}
	InternetCloseHandle(eutils);
printf("%s\n",comment);
	if ((ptr=strstr(comment,"<Id>"))==0)
	{
		InternetCloseHandle(hInternet);
		InternetCloseHandle(eutils);
		return "";
	}
	if (sscanf(ptr,"<Id>%[^<]",id)!=1)
	if ((ptr=strstr(comment,"<Id>"))==0)
	{
		InternetCloseHandle(hInternet);
		InternetCloseHandle(eutils);
		return "";
	}
	sprintf(query,
	"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s&retmode=xml",
	id);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	if (InternetReadFile(eutils,comment,39999,(LPDWORD)&nBytesRead)!=TRUE)
	if ((ptr=strstr(comment,"<Id>"))==0)
	{
		InternetCloseHandle(hInternet);
		InternetCloseHandle(eutils);
		return "";
	}
	InternetCloseHandle(eutils);
printf("%s\n",comment);
	if ((ptr=strstr(comment,"Name=\"Summary\""))==0 || (ptr=strchr(ptr,'>'))==0 || (eptr=strchr(ptr,'<'))==0)
	{
		InternetCloseHandle(hInternet);
		return "";
	}
	*eptr='\0';
	InternetCloseHandle(hInternet);
	return ptr+1;
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

char freqBuff[20000],genoBuff[20000];
int main(int argc,char *argv[])
{
	char line[200],geneName[20],*gcptr,*pmcptr;
	float mlpThreshold,mlp10,recMlp10,HWEMlp,HWEContMlp,HOMMlp,HOMContMlp;
	FILE *fo;
	int a;
	mlpThreshold=2;
	if (argc<2 || sscanf(argv[1],"gva.all.%[^.]",geneName)==1)
	{
		printf("Usage: noteSAO summary.txt *.sao\n");
		return 1;
	}
	fo=fopen(argv[1],"r");
	if (fo==0)
	{
		fo=fopen(argv[1],"w");
		fprintf(fo,"Gene\tInteresting\tMLP\trecMLP\tHWECase\tHWECont\tHOMCase\tHOMCont\tFreqs\tGenos\tSummary\tPubMed\n");
	}
	fclose(fo);
	fo=fopen(argv[1],"a");
	for (a=2;a<argc;++a)
	{
	if (sscanf(argv[a],"gva.all.%[^.]",geneName)!=1)
//	if (sscanf(argv[a],"gva.%[^.]",geneName)!=1)
		continue;
	strupr(geneName);
	getMLP(argv[a],&mlp10,&recMlp10,&HWEMlp,&HWEContMlp,&HOMMlp,&HOMContMlp);
	if (mlp10>mlpThreshold || recMlp10>mlpThreshold || HWEMlp-HWEContMlp>mlpThreshold|| HOMMlp-HOMContMlp>mlpThreshold)
	{
		getFreqs(argv[a],freqBuff);
		getGenos(argv[a],genoBuff);
		pmcptr=getPubMedComments(geneName);
		gcptr=getGeneComment(geneName);
		fprintf(fo,"%s\t\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%7.2f\t%s\t%s\t%s\t%s\n",
			geneName,mlp10,recMlp10,HWEMlp,HWEContMlp,HOMMlp,HOMContMlp,
			freqBuff,genoBuff,gcptr,pmcptr);

	}
	}
	return 0;
}