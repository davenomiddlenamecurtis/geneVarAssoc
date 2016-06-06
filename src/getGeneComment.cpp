#include <windows.h>
#include <Wininet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>

// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=(brca1%5BGene%20Name%5D)%20AND%20homo%20sapiens%5BOrganism%5D
// http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=12190&retmode=xml

char comment[40000];
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
	if ((ptr=strstr(comment,"Name=\"Summary\""))==0 || (ptr=strchr(ptr,'>'))==0 || (eptr=strchr(ptr,'<'))==0)
	{
		InternetCloseHandle(hInternet);
		return "";
	}
	*eptr='\0';
	InternetCloseHandle(hInternet);
	return ptr+1;
}

int main(int argc,char *argv[])
{
	if (InternetAttemptConnect(0)!=ERROR_SUCCESS)
		return (1);
	printf(getGeneComment(argv[1]));
	return 0;
}

#if 0
char *getGeneComment(char *geneName)
{
	FILE *fi,*fg;
	char line[1001],strippedLine[1001],fn[40],path[60],comment[40000];
	long locPos;
	fi=fopen(INDEX,"r");
	while (fgets(line,1000,fi) && strncmp(line,geneName,strlen(geneName)))
		;
	sscanf(line,"%*s %s %ld",fn,&locPos);
	fclose(fi);
	sprintf(path,"%s%s",REFDIR,fn);
	fg=fopen(path,"rb");
	fseek(fg,locPos+5,SEEK_SET);
	*comment='\0';
	do {
		if (!fgets(line,1000,fg) || !strncmp(line,"LOCUS",5))
			return comment;
	}	while (strncmp(line,"COMMENT",5));
	while (!strstr(line,"Summary:"))
	{
		if (!fgets(line,1000,fg) || !strncmp(line,"LOCUS",5)|| !strncmp(line,"PRIMARY",7))
			return comment;
	}
	do {
		sscanf(line," %[^\n^\r]",strippedLine);
		strcat(comment,strippedLine);
		strcat(comment," ");
	} while (fgets(line,1000,fg) && strncmp(line,"LOCUS",5) && strncmp(line,"PRIMARY",7));
	fclose (fg);
	return comment;
}

int main(int argc,char *argv[])
{
	printf(getGeneComment(argv[1]));
}
#endif