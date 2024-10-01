#if 0
Copyright 2018 David Curtis

This file is part of the geneVarAssoc package.

geneVarAssoc is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

geneVarAssoc is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with geneVarAssoc.If not, see <http://www.gnu.org/licenses/>.
#endif

#include "vcfLocusFile.hpp"
#include <ctype.h>

#define BCOPY(src,dest) memcpy(src,dest,sizeof(src))
#define BREAD(buff,fp) fread(buff,sizeof(buff),1,fp)
#define BWRITE(buff,fp) fwrite(buff,sizeof(buff),1,fp)
// these are only for arrays, not for pointers or addresses

int vcfLocalLocus::typeSpecificCopy(localLocus *s)
{
	vcfLocalLocus *src=(vcfLocalLocus *)s;
	if (localLocus::typeSpecificCopy(src)==0)
		return 0;
	GTpos=src->GTpos;
	GPpos=src->GPpos;
	GQpos=src->GQpos;
	ADpos = src->ADpos;
	DPpos = src->DPpos;
	qual=src->qual;
	BCOPY(info,src->info);
	BCOPY(format,src->format);
	nFieldsToSkip=src->nFieldsToSkip;
	return 1;
}

int vcfLocalLocus::read(FILE *fp)
{
	if (localLocus::read(fp)==0)
		return 0;
	fread(&GTpos,sizeof(GTpos),1,fp);
	fread(&GPpos,sizeof(GTpos),1,fp);
	fread(&GQpos,sizeof(GQpos),1,fp);
	fread(&ADpos, sizeof(ADpos), 1, fp);
	fread(&DPpos, sizeof(DPpos), 1, fp);
	fread(&qual,sizeof(qual),1,fp);
	BREAD(info,fp);
	BREAD(format,fp);
	fread(&nFieldsToSkip,sizeof(nFieldsToSkip),1,fp);
	return 1;
}

int vcfLocalLocus::write(FILE *fp)
{
	if (localLocus::write(fp)==0)
		return 0;
	fwrite(&GTpos,sizeof(GTpos),1,fp);
	fwrite(&GPpos,sizeof(GTpos),1,fp);
	fwrite(&GQpos,sizeof(GQpos),1,fp);
	fwrite(&ADpos, sizeof(ADpos), 1, fp);
	fwrite(&DPpos, sizeof(DPpos), 1, fp);
	fwrite(&qual,sizeof(qual),1,fp);
	BWRITE(info,fp);
	BWRITE(format,fp);
	fwrite(&nFieldsToSkip,sizeof(nFieldsToSkip),1,fp);
	return 1;
}

void vcfLocalLocus::parseFormat()
{
	// find position of GT, GP, GQ, AD fields
	int i;
	char *ptr;
	ptr=format;
	i=0;
	while (*ptr)
	{
		if (*ptr == 'G')
		{
			if (ptr[1] == 'T')
				GTpos = i;
			else if (ptr[1] == 'P')
				GPpos = i;
			else if (ptr[1] == 'Q')
				GQpos = i;
		}
		else if (*ptr == 'A')
		{
			if (ptr[1] == 'D')
				ADpos = i;
		}
		else if (*ptr == 'D')
		{
			if (ptr[1] == 'P')
				DPpos = i;
		}
		do {
			++ptr;
			if (*ptr == ':')
			{
				++ptr;
				++i;
				break;
			}
		} while (*ptr);
	}
	return;
}

int vcfLocalLocus::outputVcfGenotypes(FILE *fo,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap)
{
	char *ptr,allStr[20],otherInfo[100];
	int s;
	if (fseek(f,filePos,SEEK_SET)!=0)
	{
		dcerror(99,"Failed to fseek() correctly in localLocus::outputtVcfGenotypes()");
		return 0;
	}
    if (!fgets(locusFile::buff,BUFFSIZE-1,f))
	{
		dcerror(99,"Failed read locus data after fseek() in localLocus::outputtVcfGenotypes()");
		return 0;
	}
	if (strchr(locusFile::buff, '\n') == 0)
	{
		dcerror(99, "This line is too long:\n\n%s\n\nAbove line in locus file was too long. Need to change BUFFSIZE in masterLocusFile.hpp and recompile.\n", locusFile::buff);
		return 0;
	}
	for (s=0,ptr=locusFile::buff;s<nFieldsToSkip;++s)
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	for (s=0;s<nSubs;++s)
	{
		if (sscanf(ptr,"%[^: \t]:%s",allStr,otherInfo)<1)
		{
			dcerror(99,"Not enough genotypes for number of subjects in this line: %s",locusFile::buff);
			return 0;
		}
		if (allStr[0]=='.')
		{
			fprintf(fo,".\t");
		}
		else
		{	
			if (allStr[1]=='\0') // only one allele
				fprintf(fo,"%d%s\t",alleleMap[allStr[0]-'0'],otherInfo);
			else
				fprintf(fo,"%d%c%d:%s\t",alleleMap[allStr[0]-'0'],allStr[1],alleleMap[allStr[2]-'0'],otherInfo);
		}
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	return 1;
}

int vcfLocalLocus::outputProbs(probTriple *prob,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec)
{
	char *ptr,allStr[20],*aptr,*ptr2;
	int s,i;
	float gq;
	if (fseek(f,filePos,SEEK_SET)!=0)
	{
		dcerror(99,"Failed to fseek() correctly in vcfLocalLocus::outputProbs()");
		return 0;
	}
    if (!fgets(locusFile::buff,BUFFSIZE-1,f))
	{
		dcerror(99,"Failed read locus data after fseek() in vcfLocalLocus::outputProbs()");
		return 0;
	}
	if(strchr(locusFile::buff,'\n')==0)
	{
		dcerror(99,"This line is too long:\n\n%s\n\nAbove line in locus file was too long. Need to change BUFFSIZE in masterLocusFile.hpp and recompile.\n",locusFile::buff);
		return 0;
	}

	for (s=0,ptr=locusFile::buff;s<nFieldsToSkip;++s)
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	for (s=0;s<nSubs;++s)
	{
		gq=1000;
		ptr2=ptr;
		if (*ptr!='.') for (i=0;i<GPpos;++i)
		{
			while (*ptr!=':')
			{
				if (isspace(*ptr))
				{
					dcerror(99,"Not enough genotypes for number of subjects in this line: %s",locusFile::buff);
					return 0;
				}
				++ptr;
			}
			++ptr;
		}
		sscanf(ptr,"%f,%f,%f",&prob[s][0],&prob[s][1],&prob[s][2]);
		if (GQpos>0 && *ptr2!='.')
		{
		for (i=0;i<GQpos;++i)
		{
			while (*ptr2!=':')
			{
				if (isspace(*ptr2))
				{
					dcerror(99,"Could not read GQ in this line: %s",locusFile::buff);
					return 0;
				}
				++ptr2;
			}
			++ptr2;
		}
		gq=atof(ptr2);
		}
		if (*ptr=='.' || gq<spec.GQThreshold)
			prob[s][0]=prob[s][1]=prob[s][2]=0;
		// there may be a problem in dealing with X-linked loci
		// for now, cannot apply AD or DP criteria to genotype probabilities
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	return 1;
}

int vcfLocalLocus::outputAlleles(allelePair *all,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec)
{
	char *ptr,allStr[20],*aptr,*ptr2;
	int s,i,hdFail,dFail,abFail;
	float gq,ad[2],hetDev,dp;
	if (fseek(f,filePos,SEEK_SET)!=0)
	{
		dcerror(99,"Failed to fseek() correctly in vcfLocalLocus::outputAlleles()");
		return 0;
	}
    if (!fgets(locusFile::buff,BUFFSIZE-1,f))
	{
		dcerror(99,"Failed read locus data after fseek() in vcfLocalLocus::outputAlleles()");
		return 0;
	}
	if (strchr(locusFile::buff,'\n')==0)
	{
		dcerror(99,"This line is too long:\n\n%s\n\nAbove line in locus file was too long. Need to change BUFFSIZE in masterLocusFile.hpp and recompile.\n",locusFile::buff);
		return 0;
	}

	for (s=0,ptr=locusFile::buff;s<spec.numVcfFieldsToSkip;++s) // user specified number of columns to skip
	{
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	for (s=0;s<nSubs;++s)
	{
		gq=1000;
		ptr2=ptr;
		if (*ptr!='.') for (i=0;i<GTpos;++i)
		{
			while (*ptr!=':')
			{
				if (isspace(*ptr))
				{
					dcerror(99,"Not enough genotypes for number of subjects in this line: %s",locusFile::buff);
					return 0;
				}
				++ptr;
			}
			++ptr;
		}
		aptr=allStr;
		while (*ptr!=':' && *ptr!=' ' && *ptr != '\t' && *ptr != '\n')
			*aptr++=*ptr++;
		*aptr='\0'; // copied genotype into allStr
		ptr=ptr2; // back at start of entry
		if (spec.GQThreshold>0 && GQpos>0 && *ptr2!='.') // I had a rogue vcf which said it had a GQ entry but didn't
		{
		for (i=0;i<GQpos;++i)
		{
			while (*ptr2!=':')
			{
				if (isspace(*ptr2))
				{
					dcerror(99,"Could not read GQ in this line: %s",locusFile::buff);
					return 0;
				}
				++ptr2;
			}
			++ptr2;
		}
		gq=atof(ptr2);
		}
		if (spec.depthZeroUnknown)
		{
			for (i = 0; i < DPpos; ++i)
			{
				while (*ptr2 != ':')
				{
					if (isspace(*ptr2))
					{
						dcerror(99, "Could not read DP in this line: %s", locusFile::buff);
						return 0;
					}
					++ptr2;
				}
				++ptr2;
			}
			dp = atoi(ptr2);
		}
		if (spec.depthZeroUnknown && dp==0)
			all[s][0] = all[s][1] = 0;
		// GATK is now coding unknowns as genotype 0/0 with DP=0
		else if ((allStr[0] == '.' || allStr[1]=='.') && spec.wildIfUnknown)
			all[s][0] = all[s][1] = 1;
		else if (allStr[0] == '.' || allStr[1] == '.' || gq<spec.GQThreshold)
			all[s][0]=all[s][1]=0;
		// somebody made a vcf file which had entries like 0/.
		else if (((spec.hetDevThresholdSq!=-1 || spec.ABThreshold!=-1) && allStr[0]!=allStr[2]) || spec.depthThreshold!=-1 )
		{
			ptr2=ptr;
			hdFail = dFail = abFail = 0;
			for(i=0;i<ADpos;++i)
			{
				while(*ptr2!=':')
				{
					if(isspace(*ptr2))
					{
						dcerror(99,"Could not read AD in this line: %s",locusFile::buff);
						return 0;
					}
					++ptr2;
				}
				++ptr2;
			}
			ad[0]=atof(ptr2);
			while(*ptr2!=',')
				++ptr2;
			ad[1]=atof(ptr2+1);
			dp=ad[0]+ad[1];
			if (allStr[0] != allStr[2])
			{
				if (spec.hetDevThresholdSq != -1)
				{
					hetDev = dp / 2 - ad[0]; // avoid calling sscanf, sqrt, etc.
					if (hetDev*hetDev > spec.hetDevThresholdSq*dp*0.25)
						hdFail = 1;
				}
				if (spec.ABThreshold != -1)
				{
					if (ad[1] / dp < spec.ABThreshold || ad[0] / dp < spec.ABThreshold)
						abFail = 1;
				}
			}
			if (dp < spec.depthThreshold)
				dFail = 1;
			if (hdFail || abFail || dFail) // heterozygous counts too far from expected
				all[s][0] = all[s][1] = 0;
			else
			{
				if (allStr[1] == '\0')
					allStr[2] = allStr[0];
				// single allele - just double it for now
				all[s][0] = alleleMap[allStr[0] - '0'] + 1;
				all[s][1] = alleleMap[allStr[2] - '0'] + 1;
			}
		}
		else
		{	
			if (allStr[1]=='\0')
				allStr[2]=allStr[0];
			// single allele - just double it for now
			all[s][0]=alleleMap[allStr[0]-'0']+1; // we are assuming single digit allele characters: 0, 1, 2, etc.
			all[s][1]=alleleMap[allStr[2]-'0']+1;		
		}
		while (!isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	return 1;
}

int vcfLocusFile::readHeaderInfo()
{
	char *ptr;
	int s,ch,sk;
	long fPos;
	fseek(fp,0,SEEK_SET);
	do {
		fPos = ftell(fp); // I hope this will work OK with text file
		if (!fgets(locusFile::buff, BUFFSIZE - 1, fp))
		{
			dcerror(99, "Could not find line beginning #CHROM in VCF file");
			return 0;
		}
	} while (strncmp(buff, "#CHROM", strlen("#CHROM")));
	if (strchr(locusFile::buff, '\n') == 0)
	{
		dcerror(99, "This line is too long:\n\n%s\n\nAbove line in locus file was too long. Need to change BUFFSIZE in masterLocusFile.hpp and recompile.\n", locusFile::buff);
		return 0;
	}
	fseek(fp, fPos, SEEK_SET);
	for (ch = fgetc(fp), sk = 0; sk < nFieldsToSkip; ++sk)
	{
		while (!isspace(ch))
			ch = fgetc(fp);
		while (isspace(ch))
		{
			if (ch == '\n')
				break;
			ch = fgetc(fp);
		}
		if (ch == '\n')
			break; // end of line with no entries
	}
	nSubs=0;
	while (ch != '\n')
	{
		while (!isspace(ch))
		{
			ch = fgetc(fp);
		}
		++nSubs;
		while (isspace(ch))
		{
			if (ch == '\n')
				break;
			ch = fgetc(fp);
		}
	}
return nSubs;
}

int vcfLocalLocus::input(FILE *f,FILEPOSITION *locusPosInFile,analysisSpecs const &spec)
{
	char chrStr[10],qualstr[VCFFIELDLENGTH],*ptr,*qtr,word[VCFFIELDLENGTH],*chrStrPtr,afstr[20];
	int foundOneOK=0;
	hereOK();
	while (!foundOneOK) // make this a loop so I can ignore entries which do not pass
	{
	hereOK();
	*locusPosInFile=ftell(f);
	if (!fgets(locusFile::buff,BUFFSIZE-1,f))
		return 0;
	while (!strchr(locusFile::buff,'\n')) // line was too long to fit into buff
	{
		do {
		if (!fgets(locusFile::buff,BUFFSIZE-1,f))
			return 0;
		} while (!strchr(locusFile::buff,'\n')); // eat remainder of line and discard it
	*locusPosInFile=ftell(f); // make sure calling routine knows we updated position
	if (!fgets(locusFile::buff,BUFFSIZE-1,f)) // hope this new line is not too long
		return 0;
	}

	hereOK();
	format[0]='\0';
	ptr=locusFile::buff;
	if (!scanWord(&ptr,chrStr,9))
		return 0;
	if (!scanWord(&ptr,word,VCFFIELDLENGTH-1))
		goto problemReadingLocus;
	pos=atol(word);
	if (!scanWord(&ptr,id,VCFFIELDLENGTH-1))
		goto problemReadingLocus;
	if (!scanWord(&ptr,ref,MAXALLLENGTH-1))
		; // long allele, like a deletion, just use first MAXALLLENGTH
	if (!scanWord(&ptr,alt,MAXALLLENGTH*MAXALL-1))
		;
	if (!scanWord(&ptr,qualstr,VCFFIELDLENGTH-1))
		goto problemReadingLocus;
	if (spec.numVcfFieldsToSkip < 9) // assume that PASS is blank and do not try to read it
		strcpy(filter, "PASS");
	else
	{
		if (!scanWord(&ptr, filter, VCFFIELDLENGTH - 1))
			goto problemReadingLocus;
		if (spec.skipIfNoPass)
		{
			if (strncmp(filter, "PASS", 4))
				continue; // try and read next line
			// this can be set if all files have same set of loci so ones which do not PASS can be ignored
			// if a locus passes in some files but not others, should use unknownIfNoPass instead.
		}
	}
	hereOK();
	if (!scanWord(&ptr,info,VCFFIELDLENGTH-1))
		;
//		goto problemReadingLocus;
// do not make this is a problem, it just means the info field is very long

	scanWord(&ptr,format,VCFFIELDLENGTH-1); // this one is optional
	if (*format)
	{
		parseFormat();
	}
	goto noProblemReadingLocus;
	problemReadingLocus:
		hereOK();
		dcerror(99,"Problem reading locus information from this line: %s",locusFile::buff);
			return 0;
	noProblemReadingLocus:
			hereOK();
			foundOneOK=1;
		SNP=SNP_MAYBE;
		if ((ptr=strstr(info,"VT="))!=0)
			if (strncmp(ptr+3,"SNP",3))
				SNP=SNP_NO;
			else
				SNP=SNP_YES;
		nAltAlls=0;
#if 0
		if ((ptr=strstr(info,"AF="))!=0)
		{
			char freqs[100],freqStr[20];
			int a;
			sscanf(ptr+3,"%[^;]",freqs);
			a=0;
			ptr=freqs;
			while (1)
			{
				sscanf(ptr,"%[^,]",freqStr);
				alleleFreq[a++]=atof(freqStr);
				if ((ptr=strchr(ptr,','))==0)
					break;
				++ptr;
			}
			// nAltAlls=a; the problem with this is that VCF may have only one MAF even though >1 alt allele, so best not to set it
		}
#endif
		AF=AC=AN=-1;
		sprintf(afstr,";%s=",spec.alleleFreqStr);
		if ((ptr=strstr(info,afstr))!=0)
			AF=atof(ptr+strlen(afstr));
		sprintf(afstr,";%s=",spec.alleleNumberStr);
		if ((ptr=strstr(info,afstr))!=0)
			AN=atof(ptr+strlen(afstr));
		sprintf(afstr,";%s=",spec.alleleCountStr);
		if ((ptr=strstr(info,afstr))!=0)
			AC=atof(ptr+strlen(afstr));
		if (AF==-1 && AC!=-1 && AN!=-1)
			if (AN==0)
				AF=0;
			else
				AF=AC/AN;
		if (AF == -1)
		{	
#if 0
			sprintf(afstr,";AF=");
			if ((ptr=strstr(info,afstr))!=0)
				AF=atof(ptr+strlen(afstr));
			// in 1000G, some variants have AF but not EUR_AF
#endif
			AF=0;
		}

	if (strcmp(qualstr,"."))
		qual=atof(qualstr);
	else
		qual=0;
	if (strncmp(chrStr,"chr",3)==0)
		chrStrPtr=chrStr+3;
	else
		chrStrPtr=chrStr;
	if (chrStrPtr[0]=='X')
		chr=23;
	else if (chrStrPtr[0]=='Y')
		chr=24;
	else chr=atoi(chrStrPtr);
	// worry about other stuff later
	return 1;
	}
return 0;
// I think this never happens
}

int masterLocusFile::outputMergedVcfGenotypes(FILE *fo,analysisSpecs const &spec)
{
	int locusCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c,i,col,couldUseThis,a;
	locusCount=0;
recPos=findFirstInRange(spec);
if (recPos!=0L)
{
	while (1)
	{
		testKey=index.current_key();
		if ((c=atoi(testKey))==0 || c>spec.ec)
			break;
		if (c==spec.ec && atol(testKey+3)>spec.ep)
			break;
		load(tempRecord,recPos);
		for (i=0;i<nLocusFiles;++i)
		{
			if (locusFiles[i]->fileType() != VCFFILE)
			{
				dcerror(1,"Tried to call masterLocusFile::outputMergedVcfGenotypes() but a locusFile does not have type VCFFILE");
				return 0;
			}
				if (strcmp(tempRecord.myLocalLocus[i]->filter,"UNTYPED"))
				{
					couldUseThis=i;
					for (;i<nLocusFiles;++i)
					{
						if (!strcmp(tempRecord.myLocalLocus[i]->filter,"PASS"))
							break;
					}
					if (i==nLocusFiles)
						i=couldUseThis;
				// now we have selected a file which is at least typed
				// and if one exists also has a PASS
				// all genotypes from files without a PASS will get output as unknown
				// as long as unknownIfNoPass is set
					fseek(locusFiles[i]->fp,tempRecord.locusPosInFile[i],SEEK_SET);

					for (col=0;col<3;++col)
					{
						do {
							c=fgetc(locusFiles[i]->fp);
							fputc(c,fo);
						}	while (!isspace(c));
						// assume separated by single tabs then this should work
					}
					for (;col<5;++col)
					{
						do {
							c=fgetc(locusFiles[i]->fp);
						}	while (!isspace(c));
						// just need to output ref then alt alleles comma-separated
					}
					fprintf(fo,"%s\t",tempRecord.alls[0]);
					if (tempRecord.nAlls==1)
						fprintf(fo,"%s\t",tempRecord.alls[0]); // I do not know if this can ever happen
					else
						for (a=1;a<tempRecord.nAlls;++a)
							fprintf(fo,"%s%c",tempRecord.alls[a],(a==tempRecord.nAlls-1?'\t':','));
					for (;col<DEFAULTNUMVCFFIELDSTOSKIP;++col)
					{
						do {
							c=fgetc(locusFiles[i]->fp);
							fputc(c,fo);
						}	while (!isspace(c));
					}
					break;
				}
		}
		for (i=0;i<nLocusFiles;++i)
		{
			tempRecord.outputVcfGenotypes(fo,locusFiles[i]->fp,i,nSubs[i],spec);
		}
		fprintf(fo,"\n");
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;

	}
}
return locusCount;
}

int masterLocus::outputVcfGenotypes(FILE *fo,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec)
{
	int s;
	if (!strcmp(myLocalLocus[whichFile]->filter,"UNTYPED"))
	{
		for (s=0;s<nSubs;++s)
			fprintf(fo,".\t");
	}
	else if (strcmp(myLocalLocus[whichFile]->filter,"PASS") && spec.unknownIfNoPass)
			for (s=0;s<nSubs;++s)
			fprintf(fo,".\t");
	else if (myLocalLocus[whichFile]->myType()==VCFFILE)
		return ((vcfLocalLocus *)(myLocalLocus[whichFile]))->outputVcfGenotypes(fo,f,locusPosInFile[whichFile],nSubs,alleleMapping[whichFile]);
return 1;
}


