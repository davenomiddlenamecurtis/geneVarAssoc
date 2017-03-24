#include <windows.h>
#include <Wininet.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include "dcerror.hpp"
#include <assert.h>

/*
NB results from UCSC / DAS seem to be 0 based rather than 1-based
*/

char buffer[400000];

HINTERNET openConnection()
{
	HINTERNET hInternet;
	DWORD timeOut;
	timeOut=180000;
	hInternet=InternetOpen("iexplorer.exe",INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL, 0 );
	InternetSetOption(hInternet,INTERNET_OPTION_CONNECT_TIMEOUT,&timeOut,sizeof(timeOut));
	InternetSetOption(hInternet,INTERNET_OPTION_SEND_TIMEOUT,&timeOut,sizeof(timeOut));
	InternetSetOption(hInternet,INTERNET_OPTION_RECEIVE_TIMEOUT,&timeOut,sizeof(timeOut));
	return hInternet;
}

int rsNameFromPos(char *pos, char *output)
{
	char chrStr[100],*ptr,localPos[100],query[1000],rsName[100];
	int found,rv,totalDataLength;
	HINTERNET hInternet,eutils;
	LPCTSTR lpszAgent,lpszServerName,queryPtr;
	DWORD dwContext,nBytesRead;
	DWORD dwRequestFlags = INTERNET_FLAG_NO_UI   // no UI please
            |INTERNET_FLAG_NO_AUTH           // don't authenticate
            |INTERNET_FLAG_PRAGMA_NOCACHE    // do not try the cache or proxy
            |INTERNET_FLAG_NO_CACHE_WRITE   // don't add this to the IE cache
			|INTERNET_FLAG_HYPERLINK
			|INTERNET_FLAG_RELOAD
	;
	dwContext=1;
	lpszAgent=LPCTSTR("iexplorer.exe");
	lpszServerName=LPCTSTR("http://genome.ucsc.edu");
	hInternet=InternetOpen("iexplorer.exe",INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL, 0 );
	strcpy(localPos,pos);
	ptr=strchr(localPos,':');
	if (ptr==0)
		dcerror(1,"Did not understand format of position s\n",localPos);
	*ptr='\0';
	strcpy(chrStr,localPos);
	if (atoi(chrStr)==23)
		strcpy(chrStr,"X");
	else if (atoi(chrStr)==24)
		strcpy(chrStr,"Y");
	sprintf(query,"http://genome.ucsc.edu/cgi-bin/das/hg19/features?segment=%s:%s,%s;type=snp138",chrStr,ptr+1,ptr+1);
	// puts(query);
	queryPtr=LPCTSTR(query);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	totalDataLength=0;
	do {
		rv=InternetReadFile(eutils, buffer+totalDataLength, 399999, (LPDWORD)&nBytesRead);
		totalDataLength+=nBytesRead;
	} while (rv!=TRUE || nBytesRead!=0);
	buffer[totalDataLength]='\0';
	InternetCloseHandle(eutils);
	InternetCloseHandle(hInternet);
	found=0;
	ptr=buffer;
	// printf("%s\nLength = %d\n",buffer,strlen(buffer));
	if ((ptr = strstr(buffer, "FEATURE id=")) != 0)
	{
		sscanf(ptr+strlen("FEATURE id=")+1,"[^.]",rsName);
		sprintf(output,"%s",rsName);
		found=1;
	}
	if (found==0)
		sprintf(output,"%s\t-",pos);
	return 1;
}

int rsNameFromPos_dbSNP(char *pos, char *output)
{
	char chrStr[100],*ptr,localPos[100],query[1000];
	int found,rsnum,rv,totalDataLength;
	HINTERNET hInternet,eutils;
	LPCTSTR lpszAgent,lpszServerName,queryPtr;
	DWORD dwContext,nBytesRead;
	DWORD dwRequestFlags = INTERNET_FLAG_NO_UI   // no UI please
            |INTERNET_FLAG_NO_AUTH           // don't authenticate
            |INTERNET_FLAG_PRAGMA_NOCACHE    // do not try the cache or proxy
            |INTERNET_FLAG_NO_CACHE_WRITE   // don't add this to the IE cache
			|INTERNET_FLAG_HYPERLINK
			|INTERNET_FLAG_RELOAD
	;
	dwContext=1;
	lpszAgent=LPCTSTR("iexplorer.exe");
	lpszServerName=LPCTSTR("http://www.ncbi.nlm.nih.gov");
	hInternet=InternetOpen("iexplorer.exe",INTERNET_OPEN_TYPE_PRECONFIG, NULL, NULL, 0 );
	strcpy(localPos,pos);
	ptr=strchr(localPos,':');
	if (ptr==0)
		dcerror(1,"Did not understand format of position s\n",localPos);
	*ptr='\0';
	strcpy(chrStr,localPos);
	if (atoi(chrStr)==23)
		strcpy(chrStr,"X");
	else if (atoi(chrStr)==24)
		strcpy(chrStr,"Y");
	sprintf(query,"http://www.ncbi.nlm.nih.gov/snp?term=(%s[Chromosome]) AND (%s[Base Position]) AND (\"homo sapiens\"[Organism])",chrStr,ptr+1);
	// puts(query);
	queryPtr=LPCTSTR(query);
	eutils=InternetOpenUrl(hInternet,queryPtr,NULL,0,dwRequestFlags,0);
	totalDataLength=0;
	do {
		rv=InternetReadFile(eutils, buffer+totalDataLength, 399999, (LPDWORD)&nBytesRead);
		totalDataLength+=nBytesRead;
	} while (rv!=TRUE || nBytesRead!=0);
	buffer[totalDataLength]='\0';
	InternetCloseHandle(eutils);
	InternetCloseHandle(hInternet);
	found=0;
	ptr=buffer;
	// printf("%s\nLength = %d\n",buffer,strlen(buffer));
	while ((ptr = strstr(ptr, "rs=")) && !found)
	{
		--ptr;
		if (isalnum(*ptr))
			ptr+=2;
		else
		{
			if (sscanf(ptr + 4, "%d", &rsnum) == 1)
			{
				found = 1;
				sprintf(output,"%s\trs%d",pos,rsnum);
			}
		}
	}
	if (found==0)
		sprintf(output,"%s\t-",pos);
	return 1;
}

void startQuery(char *inputString)
{
	sprintf(inputString,"http://genome.ucsc.edu/cgi-bin/das/hg19/features?");
}

void addPosToQuery(char *inputString, char *pos)
{
	char chrStr[100],localPos[100],*ptr;
	strcpy(localPos,pos);
	ptr=strchr(localPos,':');
	if (ptr==0)
		dcerror(1,"Did not understand format of position s\n",localPos);
	*ptr='\0';
	strcpy(chrStr,localPos);
	if (atoi(chrStr)==23)
		strcpy(chrStr,"X");
	else if (atoi(chrStr)==24)
		strcpy(chrStr,"Y");
	sprintf(strchr(inputString,'\0'),"segment=%s:%s,%s;",chrStr,ptr+1,ptr+1);
}

void endQuery(char *inputString)
{
	sprintf(strchr(inputString,'\0'),"type=snp138");
}

int getHTMLResult(HINTERNET h, char *query, char *buffer)
{
	int totalDataLength,rv;
	HINTERNET link;
	LPCTSTR queryPtr;
	DWORD nBytesRead,err;
	DWORD dwRequestFlags = INTERNET_FLAG_NO_UI   // no UI please
            |INTERNET_FLAG_NO_AUTH           // don't authenticate
            |INTERNET_FLAG_PRAGMA_NOCACHE    // do not try the cache or proxy
            |INTERNET_FLAG_NO_CACHE_WRITE   // don't add this to the IE cache
			|INTERNET_FLAG_HYPERLINK
			|INTERNET_FLAG_RELOAD ;
	queryPtr=LPCTSTR(query);
	link=InternetOpenUrl(h,queryPtr,NULL,0,dwRequestFlags,0);
	if (link == 0)
	{
		LPVOID lpMsgBuf;
		err=GetLastError();
		FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | 
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        err,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0, NULL );
		printf("Failed to run query URL, error number %d:%s",(int)err,lpMsgBuf);
		return 0;
	}
	totalDataLength=0;
	do {
		rv=InternetReadFile(link, buffer+totalDataLength, 399999, (LPDWORD)&nBytesRead);
		totalDataLength+=nBytesRead;
	} while (rv!=TRUE || nBytesRead!=0);
	buffer[totalDataLength]='\0';
	InternetCloseHandle(link);
	return 1;
}

void closeConnection(HINTERNET hInternet)
{
		InternetCloseHandle(hInternet);
}

void parseResults(char *outputString,char *buffer)
{
	char *ptr,chr[10],pos[20],rsName[20];
	ptr=buffer;
	while ((ptr = strstr(ptr, "FEATURE id=")) != 0)
	{
		sscanf(ptr + strlen("FEATURE id=") + 1, "%[^.].%[^.].%[^\"]", rsName, chr, pos);
		sprintf(strchr(outputString, '\0'), "%s:%s %s\n",chr,pos,rsName);
		++ptr;
	}
}

#define BATCHSIZE 10
int main(int argc, char *argv[])
{
	int b;
	FILE *fi,*fo;
	HINTERNET hInternet;
	char inputString[1000],outputString[1000],pos[100];
	assert(argc==3);
	assert((fi=fopen(argv[1],"r"))!=0);
	hInternet=openConnection();
	*outputString='\0';
	while (1)
	{
		*inputString='\0';
		startQuery(inputString);
		for (b = 0; b < BATCHSIZE; ++b)
		{
			if (fscanf(fi,"%s",pos)!=1)
				if (b==0)
					goto done;
				else
					break;
			addPosToQuery(inputString,pos);
		}
		endQuery(inputString);
		getHTMLResult(hInternet,inputString,buffer);
		parseResults(outputString,buffer);

	}
	done:
	closeConnection(hInternet);
	fclose(fi);
	assert((fo=fopen(argv[2],"w"))!=0);
	fprintf(fo,outputString);
	fclose(fo);
}

#if 0
int main(int argc, char *argv[])
{
	char pos[100],output[100];
	FILE *fi,*fo;
	if (argc==1)
	while (scanf("%s ", pos) == 1 && rsNameFromPos(pos, output))
		printf("%s\n", output); 
	else if (argc == 3)
	{
		fi=fopen(argv[1],"r");
		fo=fopen(argv[2],"w");
		while (fscanf(fi,"%s ",pos)==1 && rsNameFromPos(pos,output))
			fprintf(fo,"%s\n",output);
	}
	else if (argc==2)
		if (rsNameFromPos(argv[1],output))
			printf("%s\n",output);
	return 0;
}
#endif
/*
http://genome.ucsc.edu/cgi-bin/das/hg19/features?segment=3:196792163,196792163;segment=3:196792663,196792663;type=snp138


*/