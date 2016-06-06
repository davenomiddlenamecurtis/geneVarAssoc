#include "hapsLocusFile.hpp"
#include <ctype.h>

#define STARTOFFILE -99
// special value because btree cannot store values of 0

int hapsLocalLocus::outputAlleles(allelePair *all, FILE *f, long filePos, int nSubs, int *alleleMap, analysisSpecs const &spec)
{
	char *ptr;
	int s;
	if (filePos==STARTOFFILE)
		filePos=0;
	if (fseek(f,filePos,SEEK_SET)!=0)
	{
		dcerror(99,"Failed to fseek() correctly in hapsLocalLocus::outputAlleles()");
		return 0;
	}
    if (!fgets(locusFile::buff,BUFFSIZE-1,f))
	{
		dcerror(99,"Failed read locus data after fseek() in hapsLocalLocus::outputAlleles()");
		return 0;
	}
	for (s=0,ptr=locusFile::buff;s<5;++s)
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	for (s = 0; s < nSubs; ++s)
	{
		if (sscanf(ptr, "%d %d", &all[s][0], &all[s][1]) != 2)
		{
			dcerror(99,"Could not read enough genotypes in hapsLocalLocus::outputAlleles() from this line:\n%s\n",locusFile::buff);
			return 0;
		}
		++all[s][0]; // change 01 to 12 coding
		++all[s][1];
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}

	return 1;
}

int hapsLocalLocus::input(FILE *f, long *locusPosInFile, analysisSpecs const &spec)
{
	char *ptr,chrStr[10],*chrStrPtr,posStr[20];
	*locusPosInFile=ftell(f);
	if (*locusPosInFile==0)
		*locusPosInFile=STARTOFFILE;
	if (!fgets(locusFile::buff,BUFFSIZE-1,f))
		return 0;
	if (!strchr(locusFile::buff, '\n')) // line was too long to fit into buff
	{
		dcerror(99,"Line in HAPS file is too long to read into BUFFSIZE - will need to recompile");
		return 0;
	}
	ptr=locusFile::buff;
	if (!scanWord(&ptr,chrStr,9))
		return 0;
	if (!scanWord(&ptr,rsName,20))
		return 0;
	if (!scanWord(&ptr,posStr,20))
		return 0;
	if (!scanWord(&ptr,ref,MAXALLLENGTH-1))
		return 0;
	if (!scanWord(&ptr,alt,MAXALLLENGTH-1))
		return 0;

	if (strncmp(chrStr,"chr",3)==0)
		chrStrPtr=chrStr+3;
	else
		chrStrPtr=chrStr;
	if (chrStrPtr[0]=='X')
		chr=23;
	else if (chrStrPtr[0]=='Y')
		chr=24;
	else chr=atoi(chrStrPtr);
	pos=atol(posStr);

	strcpy(filter,"PASS");

	return 1;
}
int hapsLocusFile::readHeaderInfo()
{
	int i,s;
	char *ptr;
	fseek(fp,0,SEEK_SET);
	if (!fgets(locusFile::buff,BUFFSIZE-1,fp))
	{
		dcerror(99,"Could not read first line of HAPS file");
		return 0;
	}
	for (ptr=locusFile::buff,s=0;s<5;++s)
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
		if (*ptr=='\0')
			break; // end of line with no entries
	}
	nSubs=0;
	for (;*ptr;++nSubs) // count subjects - each is one pair of alleles
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	fseek(fp,0,SEEK_SET); // be ready to read the first locus
	return nSubs;
}

int hapsLocusFile::outputSubNames(strEntry *subName, analysisSpecs &spec)
{
	return 0;
}

#define BCOPY(src,dest) memcpy(src,dest,sizeof(src))
#define BREAD(buff,fp) fread(buff,sizeof(buff),1,fp)
#define BWRITE(buff,fp) fwrite(buff,sizeof(buff),1,fp)

int hapsLocalLocus::read(FILE *fp)
{
	if (localLocus::read(fp)==0)
		return 0;
	BREAD(rsName,fp);
}
int hapsLocalLocus::write(FILE *fp)
{
	if (localLocus::write(fp)==0)
		return 0;
	BWRITE(rsName,fp);
}

int hapsLocalLocus::typeSpecificCopy(localLocus *s)
{
	hapsLocalLocus *src=(hapsLocalLocus *)s;
	if (localLocus::typeSpecificCopy(src)==0)
		return 0;
	BCOPY(rsName,src->rsName);
	return 1;
}
