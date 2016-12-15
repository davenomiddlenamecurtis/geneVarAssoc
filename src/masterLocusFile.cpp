#include "masterLocusFile.hpp"
#include "vcfLocusFile.hpp" // I need this so I can call constructor
#include "hapsLocusFile.hpp"
#include "geneVarParser.hpp"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <math.h>
#ifndef MSDOS
#include <unistd.h>
#endif

#define BCOPY(src,dest) memcpy(src,dest,sizeof(src))
#define BREAD(buff,fp) fread(buff,sizeof(buff),1,fp)
#define BWRITE(buff,fp) fwrite(buff,sizeof(buff),1,fp)
// these are only for arrays, not for pointers or addresses

char *locusFile::buff;

bool masterLocus::readRiskVarInfo(char *s, int withFreqs)
{
	char rest[1000],buff[1000],*ptr,*bptr;
	*rest='\0';
	ptr=s;
	if (*ptr++!=':')
		return false;

	bptr=buff;
	while (*ptr != ':' && ptr != '\0')
	{
		*bptr++=*ptr++;
	}
	*bptr='\0';
	chr=atoi(buff);
	if (*ptr++!=':')
		return false;

	bptr=buff;
	while (*ptr != ':' && ptr != '\0')
	{
		*bptr++=*ptr++;
	}
	*bptr='\0';
	pos=atol(buff);
	if (*ptr++!=':')
		return false;

	bptr=buff;
	while (*ptr != ':' && ptr != '\0')
	{
		*bptr++=*ptr++;
	}
	*bptr='\0';
	strcpy(masterID,buff);
	if (*ptr++!=':')
		return false;

	bptr=buff;
	while (*ptr != ':' && ptr != '\0')
	{
		*bptr++=*ptr++;
	}
	*bptr='\0';
	worstConsequenceType=consequenceType(atoi(buff));
	if (*ptr++!=':')
		return false;

	bptr=buff;
	while (*ptr != ':' && ptr != '\0')
	{
		*bptr++=*ptr++;
	}
	*bptr='\0';
	strcpy(quickConsequence,buff);
	if (*ptr++!=':')
		return false;

	bptr=buff;
	while (*ptr != ':' && ptr != '\0')
	{
		*bptr++=*ptr++;
	}
	*bptr='\0';
	strcpy(PolyPhen,buff);
	if (*ptr++!=':')
		return false;

	// strcpy(rest,ptr);
	// if (sscanf(s,":%d:%ld:%[^:]:%d:%[^:]:%[^:]: %s",&chr,&pos,masterID,worstConsequenceType,quickConsequence,PolyPhen,rest)<6)
	// fails if string has zero length
	//	return false;
	if (withFreqs && sscanf(ptr," :%d:%d:%d",&genoCount[0],&genoCount[1],&genoCount[2])!=3)
		return false;
	if (withFreqs==2 && sscanf(ptr," :%*d:%*d:%*d:\t:%d:%d:%d:%d:%d:%d:",
			&genoCcCount[0][0],&genoCcCount[0][1],&genoCcCount[0][2],
			&genoCcCount[1][0],&genoCcCount[1][1],&genoCcCount[1][2])!=6)
		return false;
	return true;
}

int masterLocus::writeRiskVarInfo(char *s, int withFreqs)
{
	char comment[100],*ptr;
	strcpy(comment,quickConsequence);
	for (ptr=comment;*ptr;++ptr)
				if (isspace(*ptr))
					*ptr='_';
	sprintf(s,":%d:%ld:%s:%d:%s:%s:\t",chr,pos,getID(),worstConsequenceType,comment,PolyPhen);
	if (withFreqs)
		sprintf(strchr(s,'\0'),":%d:%d:%d:\t",genoCount[0],genoCount[1],genoCount[2]);
	if (withFreqs==2)
		sprintf(strchr(s,'\0'),":%d:%d:%d:%d:%d:%d:\t",
			genoCcCount[0][0],genoCcCount[0][1],genoCcCount[0][2],
			genoCcCount[1][0],genoCcCount[1][1],genoCcCount[1][2]);
	return strlen(s);
}

const char *masterLocus::getID() //access function
{
	char *IDptr;
	int l;
	IDptr=masterID;
#if 0
	if (*masterID=='\0')
		for (l = 0; l < nLocusFiles; ++l)
		{
			if (strlen(myLocalLocus[l]->id)>strlen(IDptr))
				IDptr=myLocalLocus[l]->id;
			// just assume the longest is the best one to use
		}
// masterTD should have been set to the best ID and otherwise may pull an invalid ID from a local locus actually belongs to a different position, I think
#endif
	return IDptr;
#if 0
	char id[VCFFIELDLENGTH];
	int l;
	if (*masterID!='\0')
		return masterID;
	*id='\0';
	for (l = 0; l < nLocusFiles; ++l)
	{
		if (strlen(myLocalLocus[l]->id)>strlen(id))
			strcpy(id,myLocalLocus[l]->id);
		// just assume the longest is the best one to use
	}
	return id;
#endif
}

int masterLocusFile::getTotalSubs()
{
	int s,i;
	for (s=0,i=0;i<nLocusFiles;++i)
		s+=nSubs[i];
	return s;
}

int localLocus::typeSpecificCopy(localLocus *src)
{
	if (src->myType() != this->myType())
	{
		dcerror(1, "Trying to do typeSpecificCopy() between localLocus and one of different type");
		return 0;
	}
	chr=src->chr;
	pos=src->pos;
	BCOPY(id,src->id);
	BCOPY(id,src->id);
	BCOPY(ref,src->ref);
	BCOPY(alt,src->alt);
	BCOPY(alls,src->alls);
	BCOPY(alleleFreq,src->alleleFreq);
	nAltAlls=src->nAltAlls;
	SNP=src->SNP;
	BCOPY(filter,src->filter);
	AF=src->AF;
	BCOPY(PolyPhen,src->PolyPhen);
	return 1;
}

int localLocus::read(FILE *fp)
{
	fread(&chr,sizeof(chr),1,fp);
	fread(&pos,sizeof(pos),1,fp);
	BREAD(id,fp);
	BREAD(ref,fp);
	BREAD(alt,fp);
	BREAD(alls,fp);
	BREAD(alleleFreq,fp);
	fread(&nAltAlls,sizeof(nAltAlls),1,fp);
	fread(&SNP,sizeof(SNP),1,fp);
	BREAD(filter,fp);
	fread(&AF,sizeof(AF),1,fp);
	BREAD(PolyPhen,fp);
	return 1;
}

int localLocus::write(FILE *fp)
{
	fwrite(&chr,sizeof(chr),1,fp);
	fwrite(&pos,sizeof(pos),1,fp);
	BWRITE(id,fp);
	BWRITE(ref,fp);
	BWRITE(alt,fp);
	BWRITE(alls,fp);
	BWRITE(alleleFreq,fp);
	fwrite(&nAltAlls,sizeof(nAltAlls),1,fp);
	fwrite(&SNP,sizeof(SNP),1,fp);
	BWRITE(filter,fp);
	fwrite(&AF,sizeof(AF),1,fp);
	BWRITE(PolyPhen,fp);
	return 1;
}

bool scanWord(char **line, char *word, int maxLength, char token)
// scan a string from a line with a maximum length for the string and end up pointing at the next one
// if the string is too long, fill word with maxlength of it and point to the next word but return false
{
	char *ptr;
	int i;
	bool tooLong = false;
	ptr = *line;
	*word = '\0';
	if (token == '\0')
		while (isspace(*ptr))
			++ptr;
	else
		while (isspace(*ptr)||*ptr==token)
			++ptr;
	if (*ptr=='\0')
		return 0;
	for (i=0;i<maxLength;++i)
	{
		if (isspace(*ptr) || *ptr=='\0' || *ptr==token)
			break;
		*word++=*ptr++;
	}
	*word='\0';
	if (*ptr && !isspace(*ptr) && *ptr!=token)
		tooLong=true;
	if (token == '\0')
	{
		while (*ptr && !isspace(*ptr))
			++ptr;
		while (isspace(*ptr))
			++ptr;
	}
	else
	{
		while (*ptr && !isspace(*ptr) && *ptr!=token)
			++ptr;
		while (isspace(*ptr) || *ptr==token)
			++ptr;
	}

	*line=ptr;
	return !tooLong;
}

int masterLocus::readQueryOutput(FILE *fp)
{
	char line[200],CDS_positionStr[100];
	if (alls[0][0]=='N')
		return 1; // predictor produces no output for these
	if (nAlls<2)
		return 1;
	if (!fgets(line,199,fp))
	{
		dcerror(1,"Could not read from variant_predictor output");
		return 0;
	}
	CDS_positionStr[0]='-';
	CDS_positionStr[0]='\0';
	if (sscanf(line,"%*s %*s %*s %*s %*s %*s %s %*s %s",ensemblConsequence,CDS_positionStr)<1)
	{
		dcerror(1,"Could not find consequence in this line from variant_predictor: \n%s",line);
		return 0;
	}
	if (CDS_positionStr[0]!='-')
		sprintf(strchr(ensemblConsequence,'\0')," %s",CDS_positionStr);
	return 1;
}

int masterLocus::writePredictorQuery(FILE *fp)
{
	long start,end;
	if (alls[0][0]=='N')
		return 1; // predictor produces no output for these
	if (nAlls<2)
		return 1;
	// for now, we will only consider first of alternate alleles
	if (alls[0][0]!='-')
	{
		start=pos;
		end=pos+strlen(alls[1])-1;
	}
	else
	{
		// this may not be quite right
		// assume we have deletion and first allele of alternates gives the length
		start=pos+strlen(alls[1]);
		end=pos;
	}
	fprintf(fp,"%d\t%ld\t%ld\t%s/%s\t+\n",chr,start,end,alls[0],alls[1]);
	return 1;
}

int masterLocus::getQuickFeature(refseqGeneInfo &r)
{
	char thisGeneEffect[100];
	float oldEffectWeight;
	int e;
	if (nAlls<2)
		return 1;
	if (chr==r.getChrNum())
	{
		r.getEffect(pos,alls[0],alls[1]);
		if (quickConsequence[0] == '\0' || !strncmp(quickConsequence, "NULL", 4) || !strncmp(quickConsequence, "INTER", 5))
		{
			strcpy(quickConsequence, r.tellEffect());
			worstConsequenceType=r.tellWorstConsequence();
		}
		else
		{
			strcpy(thisGeneEffect,r.tellEffect());
			for (e=0;e<NCONSEQUENCETYPES;++e)
				if (!strncmp(quickConsequence,consequence[e].str,strlen(consequence[e].str)))
				{
					oldEffectWeight=consequence[e].weight;
					break;
				}
				if (oldEffectWeight < consequence[r.tellWorstConsequence()].weight)
				{
					strcpy(quickConsequence, r.tellEffect());
					worstConsequenceType=r.tellWorstConsequence();
				}
		}

	}

#if 0
	if (quickConsequence[0]=='\0' || !strncmp(quickConsequence,"NULL",4) || !strncmp(quickConsequence,"INTER",5))
		// nothing done yet, can use multiple genes till I find one which involves this variant
		// NB this will not work completely - for example if was upstream for other gene 
	{
		if (chr!=r.getChrNum())
		{
			strcpy(quickConsequence,"INTERGENIC");
		}
		else
		{
			r.getEffect(pos,alls[0],alls[1]);
			strcpy(quickConsequence,r.tellEffect());
		}
	}
#endif
	return 1;
}


int masterLocusFile::getQuickConsequences(refseqGeneInfo &r,analysisSpecs const &spec,int redo)
{
	int locusCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c;
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
		if (!redo) // may have already been set using a different gene
			tempRecord.quickConsequence[0]='\0';
		tempRecord.getQuickFeature(r);
		// some other program calls this and feature may already be set
		save(tempRecord,recPos);
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;
// there seems to be an error in the btree routines so that when > 50,000 returns recPos=-1
// need to fix this sometime
	}
}
return locusCount;
}

int masterLocusFile::getEnsemblConsequences(analysisSpecs const &spec)
{
	int locusCount;
	FILEPOSITION recPos;
	const char *testKey;
	char line[1000];
	int c;
	FILE *fp;
	FILEPOSITION fPos;
	locusCount=0;
recPos=findFirstInRange(spec);
if (recPos!=0L)
{
	fp=fopen("predictorQuery.txt","w");
	while (1)
	{
		testKey=index.current_key();
		if ((c=atoi(testKey))==0 || c>spec.ec)
			break;
		if (c==spec.ec && atol(testKey+3)>spec.ep)
			break;
		load(tempRecord,recPos);
		tempRecord.writePredictorQuery(fp);
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;

	}
	fclose(fp);
	unlink("predictorOutput.txt");
	system("perl variant_effect_predictor.pl -i predictorQuery.txt -o predictorOutput.txt --most_severe");
	fp=fopen("predictorOutput.txt","rb"); // binary mode can use fseek/ftell
	if (fp==NULL)
	{
		dcerror(1,"Could not open output file predictorOutput.txt from variant_effect_predictor.pl");
		return 0;
	}
	do {
		fPos=FTELL(fp);
		if (!fgets(line,999,fp))
		{
			dcerror(1,"Could not read any valid line from predictorOutput.txt");
			return 0;
		}
	} while(line[0]=='#');
	FSEEK(fp,fPos,SEEK_SET); // go back to start of first valid line
	recPos=findFirstInRange(spec);	
	while (1)
	{
		testKey=index.current_key();
		if ((c=atoi(testKey))==0 || c>spec.ec)
			break;
		if (c==spec.ec && atol(testKey+3)>spec.ep)
			break;
		load(tempRecord,recPos);
		tempRecord.readQueryOutput(fp);
		save(tempRecord,recPos);
		recPos=index.get_next();
		if (recPos==0L)
			break;
	}
	fclose(fp);
}
return locusCount;
}

#define MAXLOCIINSCOREASSOCFILE 50000
int useLocus[MAXLOCIINSCOREASSOCFILE];
float locusWeight[MAXLOCIINSCOREASSOCFILE];

int masterLocusFile::writeScoreAssocFiles(char *root, float wf, int *useFreqs, int *suppliedNSubs, int writeNames, int writeComments, int writeScoreFile, analysisSpecs &spec)
{
	return writeScoreAssocFiles(*this,root,wf,useFreqs,suppliedNSubs,writeNames,writeComments,writeScoreFile, spec);
}

int masterLocusFile::loadFirst(analysisSpecs &spec)
{
	currentRecPos=findFirstInRange(spec);
	if (currentRecPos==0L)
		return 0;
	else 
		return load(tempRecord,currentRecPos);
}

int masterLocusFile::loadNext(analysisSpecs &spec)
{
	const char *testKey;
	int c;
	currentRecPos = index.get_next();
	if (currentRecPos == 0L)
		return 0;
	testKey = index.current_key();
	if ((c = atoi(testKey)) == 0 || c > spec.ec)
		return 0;
	if (c == spec.ec && atol(testKey + 3) > spec.ep)
		return 0;
	return 
		load(tempRecord, currentRecPos);
}

int masterLocusFile::writeScoreAssocFiles(masterLocusFile &subFile,char *root, float wf, int *useFreqs, int *suppliedNSubs, int writeNames, int writeComments, int writeScorefile,analysisSpecs &spec)
// allow information about subjects to be provided by a different masterLocusFile 
// however we are assuming both files refer to identical set of subjects
{
	char fn[100],buff[1000],buff2[20],comment[1000],*ptr,alleles[MAXSTR+1],commandString[1000];
	allelePair **a;
	probTriple **p;
	int totalSub,lc,s,l,ss,i,c;
	FILE *fp;
	FILEPOSITION recPos;
	const char *testKey;
	strEntry *subName;
	if (nLocusFiles==1)
		nSubs[0]=subFile.getTotalSubs();
	else if (nLocusFiles==subFile.nLocusFiles)
		for (i=0;i<nLocusFiles;++i)
			nSubs[i]=subFile.nSubs[i];
	else
	{
		dcerror(99,"Incompatible numbers of subjects in subFile in masterLocusFile::writeOldScoreAssocFiles()");
		return 0;
	}
	openLocusFiles();
	for (i=0,totalSub=0;i<subFile.nLocusFiles;++i)
		totalSub+=subFile.nSubs[i];
	subName=(strEntry *)calloc(totalSub,sizeof(strEntry));
// hereOK();
	subFile.outputSubNames(subName,spec);
// hereOK();
	if (spec.useProbs)
	{
		assert((p = (probTriple **)calloc(MAXLOCIINSCOREASSOCFILE, sizeof(probTriple*))) != 0);
		for (l = 0; l < MAXLOCIINSCOREASSOCFILE; ++l)
			p[l] = (probTriple *)calloc(totalSub, sizeof(probTriple));
		lc = outputProbs(p, spec);
	}
	else
	{
		assert((a = (allelePair **)calloc(MAXLOCIINSCOREASSOCFILE, sizeof(allelePair*))) != 0);
		for (l = 0; l < MAXLOCIINSCOREASSOCFILE; ++l)
			a[l] = (allelePair *)calloc(totalSub, sizeof(allelePair));
		lc = outputAlleles(a, spec);
	}
// hereOK();
	sprintf(fn,"%s.dat",root);
	fp=fopen(fn,"w");
	if(spec.subPhenos.size()>0)
	{
		if(spec.phenotypes==NULL) // should have been allocated when IDsAndPhenotypesFileName read
		{
			spec.phenotypes=(int*)malloc(sizeof(int)*totalSub);
			assert(spec.phenotypes!=0);
		}
		TStrIntMap::iterator it;
		for(s=0;s<totalSub;++s)
		{
			it=spec.subPhenos.find(subName[s]);
			if (it==spec.subPhenos.end())
				spec.phenotypes[s]=-1;
			else
				spec.phenotypes[s]=it->second;
		}
	}
	for (s=0,i=0;i<subFile.nLocusFiles;++i)
	for (ss=0;ss<subFile.nSubs[i];++s,++ss)
	{
		if(spec.phenotypes&&spec.phenotypes[s]==-1)
			continue;
		fprintf(fp,"%s\t%d\t",subName[s],spec.phenotypes?spec.phenotypes[s]:subFile.cc[i]);
		for (l=0;l<lc;++l)
			if (spec.useProbs)
			{
				fprintf(fp, "%5.3f %5.3f %5.3f\t",p[l][s][0],p[l][s][1],p[l][s][2]);
			}
			else
			{
				fprintf(fp, "%d %d\t",
					(a[l][s][0]>1) ? 2 : a[l][s][0],
					(a[l][s][1] > 1) ? 2 : a[l][s][1]); // force to be biallelic
			}
		fprintf(fp,"\n");
	}
	fclose(fp);
	sprintf(commandString,"scoreassoc %s %s --numloci %d",spec.useProbs?"--gendatafile":"--gcdatafile",fn,lc);
	outputSAInfo(useLocus,locusWeight,spec);
	sprintf(fn,"%s.lf.par",root);
	fp=fopen(fn,"w");
	for (l=0;l<lc;++l)
			fprintf(fp,"%d ",useLocus[l]);
	fprintf(fp,"\n");
	fclose(fp);
	sprintf(strchr(commandString,'\0')," --locusfilterfile %s",fn);
	if (spec.doRecessiveTest)
	{
		sprintf(strchr(commandString, '\0'), " --dorecessive --minweight %f --ldthreshold %f ", spec.weightThreshold, spec.LDThreshold);
		if (spec.useHaplotypes)
			sprintf(strchr(commandString, '\0'), " --usehaps");
	}
#if 0
	fprintf(fp,"\n%f %d\n%d %d %d %d %d %d %f %f %d %d\n",
		wf,
		wFunc,
		spec.useConsequenceWeights,
		useFreqs[0],
		useFreqs[1],
		writeNames,
		writeComments,
		spec.doRecessiveTest,
		spec.weightThreshold,
		spec.LDThreshold,
		spec.useHaplotypes,
		spec.useTrios);
#endif
	if (spec.useConsequenceWeights)
	{
		sprintf(fn,"%s.lw.par",root);
		fp=fopen(fn,"w");
		for (l=0;l<lc;++l)
			fprintf(fp,"%8.5f ",locusWeight[l]);
		fprintf(fp,"\n");
		fclose(fp);
	sprintf(strchr(commandString,'\0')," --locusweightfile %s",fn);
	}
	for (i=0;i<2;++i)
		if (useFreqs[i])
		{
			float *freqs;
			freqs=new float[lc];
			// outputAltFrequencies(freqs,i,sc,sp,ec,ep);
			outputEurAltFrequencies(freqs,i,spec);
			sprintf(fn,"%s.%s.freq.par",root,i?"case":"cont");
			fp=fopen(fn,"w");
			for (l=0;l<lc;++l)
				fprintf(fp,"%8.6f ",freqs[l]);
			fprintf(fp,"\n");
			for (l=0;l<lc;++l)
				fprintf(fp,"%8d ",suppliedNSubs[i]);  
			// assume that if AF available nSubs is unknown
			// change this later to use nSubs if available
			fprintf(fp,"\n");
			fclose(fp);
			sprintf(strchr(commandString,'\0')," --%sfreqfile %s",i?"case":"cont",fn);
			delete[] freqs;
		}
// no longer dealing with names seperately
	if (writeComments)
	{
		sprintf(fn,"%s.comm.par",root);
		fp=fopen(fn,"w");
		recPos=findFirstInRange(spec);
	if (recPos!=0L)
		while (1)
		{
			testKey=index.current_key();
			if ((c=atoi(testKey))==0 || c>spec.ec)
				break;
			if (c==spec.ec && atol(testKey+3)>spec.ep)
				break;
			load(tempRecord,recPos);
			if (tempRecord.ensemblConsequence[0]!='\0')
				sprintf(comment,"%d:%ld:%s:%s",tempRecord.chr,tempRecord.pos,tempRecord.getID(),tempRecord.ensemblConsequence);
			else if (tempRecord.quickConsequence[0]!='\0')
				sprintf(comment,"%d:%ld:%s:%s",tempRecord.chr,tempRecord.pos,tempRecord.getID(),tempRecord.quickConsequence);
			else
				sprintf(comment,"%d:%ld:%s:",tempRecord.chr,tempRecord.pos,tempRecord.getID());
			for (ptr=comment;*ptr;++ptr)
				if (isspace(*ptr))
					*ptr='_';
			strcat(comment,"-");
			buff[0]='\0';
			for (i=0;i<tempRecord.nAlls;++i)
				sprintf(strchr(buff,'\0'),"%s%c",tempRecord.alls[i],i==tempRecord.nAlls-1?'\0':'/');
			strncpy(alleles,buff,MAXSTR);
			alleles[MAXSTR]='\0';
			strcat(comment,alleles);
			if (tempRecord.PolyPhen[0])
			{
				strcat(comment,":PolyPhen:");
				strcat(comment,tempRecord.PolyPhen);
			}
			fprintf(fp,"%s ",comment);
			recPos=index.get_next();
			if (recPos==0L)
				break;
		}
	fprintf(fp,"\n");
	fclose(fp);
	sprintf(strchr(commandString,'\0')," --locusnamefile %s",fn);
	}

	if (spec.useTrios)
		sprintf(strchr(commandString,'\0')," --triofile %s",spec.triosFn);
	if (spec.nExc > 0)
	{
		sprintf(fn,"%s.filter.par",root);
		fp=fopen(fn,"w");
		for (i = 0; i < spec.nExc; ++i)
			fprintf(fp, "%s\n", spec.exclusionStr[i]);
		fclose(fp);
		sprintf(strchr(commandString,'\0')," --filterfile %s",fn);
	}
	sprintf(strchr(commandString,'\0')," --outfile %s.sao",root);
	if (writeScorefile)
		sprintf(strchr(commandString,'\0')," --scorefile %s.sco",root);
#ifndef MSDOS
	sprintf(fn,"%s.sh",root);
#else
	sprintf(fn,"%s.bat",root);
#endif
	fp=fopen(fn,"w");
	fprintf(fp,"%s\n",commandString);
	fclose(fp);
	if (spec.useProbs)
	{
		for (l = 0; l < MAXLOCIINSCOREASSOCFILE; ++l)
			free(p[l]);
		free(p);
	}
	else
	{
		for (l = 0; l < MAXLOCIINSCOREASSOCFILE; ++l)
			free(a[l]);
		free(a);
	}
	free(subName);
	return 1;
}

int masterLocusFile::outputAltFrequencies(float *freqs,int cc,analysisSpecs const &spec)
{
	int locusCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c,f;
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
		freqs[locusCount]=0;
		for (f=0;f<tempRecord.myLocalLocus[cc]->nAltAlls;++f)
			freqs[locusCount]+=tempRecord.myLocalLocus[cc]->alleleFreq[f];
		// cumulative frequency of all alternative alleles
		// as provided in VCF file
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;
	}
}
return locusCount;
}

int masterLocusFile::outputEurAltFrequencies(float *freqs,int cc,analysisSpecs const &spec)
{
	int locusCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c;
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
#if 0
		freqs[locusCount]=tempRecord.myLocalLocus[cc].eurAF;
#else
		if (strcmp(tempRecord.myLocalLocus[cc]->filter,"UNTYPED"))
			freqs[locusCount]=tempRecord.myLocalLocus[cc]->AF;
		else
			freqs[locusCount]=0;
		// I am not sure if this is right but at least I can break here
#endif
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;

	}
}
return locusCount;
}

int masterLocusFile::outputSAInfo(int *useLocus,float *locusWeight,analysisSpecs const &spec)
{
	int locusCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c,i,doNotUseUntypedFreqfile;
	consequenceType cons;
	geneVarParser weightParser;
	std::list<geneVarParser *> excludeParser;
	if (spec.debug)
		weightParser.debug(stdout); // will set debug for both parsers
	locusCount=0;
	recPos=findFirstInRange(spec);
if (recPos!=0L)
{
	if (spec.weightExpression[0])
		weightParser.parse(spec.weightExpression); // only have to parse once
											// it is essential that geneVarParser::thisGene has been set!!
	if(spec.excludeExpressions.size())
	{
		for (std::list<std::string>::const_iterator it=spec.excludeExpressions.begin();it!=spec.excludeExpressions.end();++it)
		{
			geneVarParser *eP=new geneVarParser;
			eP->parse(it->c_str());
			excludeParser.push_back(eP);
		}
	}
	while (1)
	{
		testKey=index.current_key();
		if ((c=atoi(testKey))==0 || c>spec.ec)
			break;
		if (c==spec.ec && atol(testKey+3)>spec.ep)
			break;
		if (spec.useConsequenceWeights==0 && spec.consequenceThreshold==0 && spec.onlyUseSNPs==0)
		{
			useLocus[locusCount]=1;
			locusWeight[locusCount]=1.0;
		}
		else
		{
			locusWeight[locusCount]=1.0; // default if nothing else changes it
			useLocus[locusCount]=1;
			load(tempRecord,recPos);
			doNotUseUntypedFreqfile=0;
			if (spec.unknownIfUntyped)
			{
				for (i=0;i<nLocusFiles;++i)
					if (holdsFreqs[i] && !strcmp(tempRecord.myLocalLocus[i]->filter,"UNTYPED"))
						doNotUseUntypedFreqfile=1;
			}
			if (doNotUseUntypedFreqfile || (tempRecord.isSNP()!=SNP_YES && spec.onlyUseSNPs==1))
			{
				locusWeight[locusCount]=0.0;
				useLocus[locusCount]=0;
			}
		else if (spec.useConsequenceWeights==0 && spec.consequenceThreshold==0)
		{
			useLocus[locusCount]=1;
			locusWeight[locusCount]=1.0;
		}
			else 
			{
				if (spec.weightExpression[0] && spec.useConsequenceWeights)
				{
					geneVarParser::thisLocus=&tempRecord;
					locusWeight[locusCount]=(double)(*weightParser.eval());
				}
				else if (spec.useEnsembl==1 && spec.useConsequenceWeights)
				{
				if (tempRecord.ensemblConsequence[0]=='\0')
					cons=NULL_CONSEQUENCE;
				else
				{
				for (c=0;c<NCONSEQUENCETYPES;++c)
					{
						if (!strncmp(tempRecord.ensemblConsequence,consequence[c].str,strlen(consequence[c].str)))
						{
							cons=consequence[c].t;
							break;
						}
					}
				if (c==NCONSEQUENCETYPES)
						{
							dcerror(1,"Could not recognise this consequence: %s",tempRecord.ensemblConsequence);
							return 0;
						}
				}
				locusWeight[locusCount]=consequence[c].weight;
				if(spec.consequenceThreshold!=0)
					useLocus[locusCount]=(cons>=spec.consequenceThreshold)?1:0;
				}
			else if (spec.useConsequenceWeights)
			{
				if (tempRecord.quickConsequence[0]=='\0')
					cons=NULL_CONSEQUENCE;
				else
				{
				for (c=0;c<NCONSEQUENCETYPES;++c)
					{
						if (!strncmp(tempRecord.quickConsequence,consequence[c].str,strlen(consequence[c].str)))
						{
							cons=consequence[c].t;
							break;
						}
					}
				if (c==NCONSEQUENCETYPES)
						{
							dcerror(1,"Could not recognise this consequence: %s",tempRecord.ensemblConsequence);
							return 0;
						}
				}
				locusWeight[locusCount]=consequence[c].weight;
				if(spec.consequenceThreshold!=0)
					useLocus[locusCount]=(cons>=spec.consequenceThreshold)?1:0;
			}
			}
		}
		if (excludeParser.size())
		{	
			double rv;
			geneVarParser::thisLocus=&tempRecord;
			geneVarParser::thisWeight=locusWeight[locusCount];
			for(std::list<geneVarParser *>::iterator it=excludeParser.begin();it!=excludeParser.end();++it)
			{
				rv=*((*it)->eval());
				if(rv!=0)
					useLocus[locusCount]=0;
				break;
			}
		}
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;
	}
	if(excludeParser.size())
	{
		for(std::list<geneVarParser *>::iterator it=excludeParser.begin();it!=excludeParser.end();++it)
		{
			geneVarParser *gP=*it;
			delete gP;
		}
	}
}
return locusCount;
}

int masterLocusFile::gotoFirstInRange(analysisSpecs &spec)
{
	FILEPOSITION recPos;
	recPos=findFirstInRange(spec);
	if (recPos==0 || load(tempRecord,recPos)==0)
		return 0;
	else
		return 1;
}

int masterLocusFile::gotoNextInRange(analysisSpecs &spec)
{
	const char *testKey;
	int c;
	FILEPOSITION recPos;
	testKey=index.current_key();
	if ((c=atoi(testKey))==0 || c>spec.ec || (c==spec.ec && atol(testKey+3)>spec.ep))
		return 0;
	recPos=index.get_next();
	if (recPos==0 || load(tempRecord,recPos)==0)
		return 0;
	else
		return 1;
}

int masterLocusFile::outputCurrentProbs(probTriple *prob, analysisSpecs &spec)
{
	int i,subCount;
	if (spec.unknownIfUntyped==0)
		{
			spec.altIsCommon=0;
			for (i=0;i<nLocusFiles;++i)
				if (tempRecord.myLocalLocus[i]->AF>0.5)
					spec.altIsCommon=1;
		}
	for (subCount=0,i=0;i<nLocusFiles;++i)
		{
			tempRecord.outputProbs(prob+subCount,locusFiles[i]->fp,i,nSubs[i],spec);
			subCount+=nSubs[i];
		}
	return 1;
}

int masterLocusFile::outputCurrentAlleles(allelePair *all, analysisSpecs &spec)
{
	int i,subCount;
	if (spec.unknownIfUntyped==0)
		{
			spec.altIsCommon=0;
			for (i=0;i<nLocusFiles;++i)
				if (tempRecord.myLocalLocus[i]->AF>0.5)
					spec.altIsCommon=1;
		}
	for (subCount=0,i=0;i<nLocusFiles;++i)
		{
			tempRecord.outputAlleles(all+subCount,locusFiles[i]->fp,i,nSubs[i],spec);
			subCount+=nSubs[i];
		}
	return 1;
}

int masterLocusFile::outputProbs(probTriple **prob, analysisSpecs &spec)
{
	int locusCount, subCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c, i, altIsCommon;
	locusCount = 0;
// hereOK();
	if (gotoFirstInRange(spec))
	do
	{
// hereOK();
		outputCurrentProbs(prob[locusCount++], spec);
	} while (gotoNextInRange(spec));
// hereOK();
	return locusCount;
}

int masterLocusFile::outputAlleles(allelePair **all, analysisSpecs &spec)
{
	int locusCount, subCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c, i, altIsCommon;
	locusCount = 0;
// hereOK();
	if (gotoFirstInRange(spec))
	do
	{
// hereOK();
		outputCurrentAlleles(all[locusCount++], spec);
	} while (gotoNextInRange(spec));
// hereOK();
	return locusCount;
}

int masterLocusFile::outputMergedVCFHeader(FILE *fo)
// no error-checking, fingers crossed
{
	int i,c,col;
	char *ptr;
	for (i=0;i<nLocusFiles;++i)
	{
		if (locusFiles[i]->fileType() != VCFFILE)
		{
			dcerror(1,"Tried to call masterLocusFile::outputMergedVCFHeader() but a locusFile does not have type VCFFILE");
			return 0;
		}
		FSEEK(locusFiles[i]->fp,0L,SEEK_SET);
		while (fgets(locusFile::buff,BUFFSIZE-1,locusFiles[i]->fp),locusFile::buff[1]=='#')
		{
			if (i==0)
				fprintf(fo,"%s",locusFile::buff);
		}
		for (ptr=locusFile::buff,col=0;col<DEFAULTNUMVCFFIELDSTOSKIP;++col)
		{
			while (!isspace(*ptr))
			{
				if (i==0)
					fputc(*ptr,fo);
				++ptr;
			}
			if (i==0)
				fputc(*ptr,fo); // tab, I hope
			++ptr;
		}
		for (col=0;col<nSubs[i];++col)
		{
			while (!isspace(*ptr))
			{
				fputc(*ptr,fo);
				++ptr;
			}
			fputc('\t',fo);
			++ptr; // assume just single tab
		}
	}
	fputc('\n',fo);
	return 1;
}

int masterLocusFile::openLocusFiles()
{
	int i;
	for (i=0;i<nLocusFiles;++i)
	{
		locusFiles[i]->fp=fopen(lfFileNames[i],"rb");
		if (locusFiles[i]->fp==0)
		{
			dcerror(99,"Could not open file %s",lfFileNames[i]);
			closeLocusFiles();
			return 0;
		}
	}
return 1;
}

int masterLocusFile::closeLocusFiles()
{
	int i;
	for (i=0;i<nLocusFiles;++i)
	{
		if (locusFiles[i]!=0 && locusFiles[i]->fp)
		{
			fclose(locusFiles[i]->fp);
			locusFiles[i]->fp=0;
		}
	}
return 1;
}


int masterLocus::outputProbs(probTriple *prob,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec)
{
	int s;
	if (!strcmp(myLocalLocus[whichFile]->filter,"UNTYPED"))
	{
		if (spec.unknownIfUntyped)
			for (s=0;s<nSubs;++s)
				prob[s][0]=prob[s][1]=prob[s][2]=0; // may change this later
		else
			for (s = 0; s < nSubs; ++s)
			{
				if (spec.altIsCommon)
				{
					prob[s][0]=prob[s][1]=0;
					prob[s][2]=1;
				}
				else
				{
					prob[s][2]=prob[s][1]=0;
					prob[s][0]=1;
				}
			}
	}
	else if (strcmp(myLocalLocus[whichFile]->filter,"PASS") && spec.unknownIfNoPass)
			for (s=0;s<nSubs;++s)
				prob[s][0]=prob[s][1]=prob[s][2]=0; // may change this later
	else
		return myLocalLocus[whichFile]->outputProbs(prob,f,locusPosInFile[whichFile],nSubs,alleleMapping[whichFile],spec);
return 1;
}

int masterLocus::outputAlleles(allelePair *all,FILE *f,int whichFile,int nSubs,analysisSpecs const &spec)
{
	int s;
	if (!strcmp(myLocalLocus[whichFile]->filter,"UNTYPED"))
	{
		if (spec.unknownIfUntyped)
			for (s=0;s<nSubs;++s)
				all[s][0]=all[s][1]=0;
		else
			for (s=0;s<nSubs;++s)
				all[s][0]=all[s][1]=spec.altIsCommon?2:1; // mark all as homozygous for common allele, pretty unsatisfactory because there should be a call
	}
	else if (strcmp(myLocalLocus[whichFile]->filter,"PASS") && spec.unknownIfNoPass)
			for (s=0;s<nSubs;++s)
				all[s][0]=all[s][1]=0;
	else
		return myLocalLocus[whichFile]->outputAlleles(all,f,locusPosInFile[whichFile],nSubs,alleleMapping[whichFile],spec);
return 1;
}

int masterLocusFile::print(FILE *fp,analysisSpecs const &spec)
{
	FILEPOSITION recPos;
	const char *testKey;
	int c;
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
		tempRecord.print(fp);
		recPos=index.get_next();
		if (recPos==0L)
			break;

	}
return 1;
}
else
	return 0;
}

int masterLocusFile::printFeatures(FILE *fp,analysisSpecs const &spec,int showFreq)
{
	FILEPOSITION recPos;
	const char *testKey;
	int c;
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
		tempRecord.printFeatures(fp,showFreq);
		recPos=index.get_next();
		if (recPos==0L)
			break;

	}
return 1;
}
else
	return 0;
}

int masterLocus::print(FILE *fp)
{
	int a,f;
	fprintf(fp,"chr = %2d pos = %9ld SNP = %d nLocusFiles = %d\nref = %s\nalt = %s\n",
		chr,pos,SNP,nLocusFiles,ref,alt);
	fprintf(fp,"Ensembl consequence: %s\n",ensemblConsequence);	
	fprintf(fp,"Quick consequence: %s\n",quickConsequence);	
	fprintf(fp,"PolyPhen: %s\n",PolyPhen);	
	fprintf(fp,"nAlls = %d: ",nAlls);
	for (a=0;a<nAlls;++a)
		fprintf(fp,"%s ",alls[a]);
	fprintf(fp,"\nMapping:\n");

	for (f=0;f<nLocusFiles;++f)
	{
		for (a=0;a<MAXALL;++a)
			fprintf(fp,"%2d ",alleleMapping[f][a]);
		fprintf(fp,"\n");
	}
fprintf(fp,"\n");
return 1;
}

int masterLocus::printFeatures(FILE *fp,int showFreq)
{
	int a,f;
	float freq;
	fprintf(fp,"%d_%d_%s",chr,pos,alls[0]);
	for (a=1;a<nAlls;++a)
		fprintf(fp,"/%s",alls[a]);	
	fprintf(fp,"\t%s",ensemblConsequence);
	if (showFreq!=-1)
	{
		freq=0;
		for (f=0;f<myLocalLocus[showFreq]->nAltAlls;++f)
			freq+=myLocalLocus[showFreq]->alleleFreq[f];
		fprintf(fp,"\t%6.4f",freq);
	}
	fprintf(fp,"\t%s\t%s\n",quickConsequence,PolyPhen);
return 1;
}

int masterLocusFile::addLocusFile(char *fn, locusFileType t)
{
	int i;
	for (i=0;i<nLocusFiles;++i)
	{
		if (lfFileNames[i][0]=='\0')
		{
			strcpy(lfFileNames[i],fn);
			currentLocusFile=i;
			fileTypes[i]=t;
			switch (t)
			{
			case VCFFILE:
				locusFiles[i]=(LOCUSFILEPTR)new vcfLocusFile;
				break;
			case SHAPEITHAPSFILE:
				locusFiles[i]=(LOCUSFILEPTR)new hapsLocusFile;
				break;
			default:
				dcerror(99,"Unrecognised locusFileType %d in masterLocusFile::addLocusFile",t);
				return 0;
			}
		break;
		}
		else if (!strcmp(lfFileNames[i],fn))
		{
			currentLocusFile=i;
			dcerror(99,"Cannot add %s to master locus file because it has already been added",fn);
			return 0;
		}
	}
	if (i==nLocusFiles)
	{
		dcerror(99,"Cannot add %s to master locus file because it already contains nLocusFiles=%d",fn,nLocusFiles);
		return 0;
	}
}

int masterLocusFile::readLocusFileEntries(char *fn,analysisSpecs const &spec,int aff)
{
	int i,s,nread;
	char *ptr;
#if 0
	for (i=0;i<nLocusFiles;++i)
	{
		if (lfFileNames[i][0]=='\0')
		{
			strcpy(lfFileNames[i],fn);
			currentLocusFile=i;
			break;
		}
		else if (!strcmp(lfFileNames[i],fn))
		{
			currentLocusFile=i;
			break;
		}
	}
	if (i==nLocusFiles)
	{
		dcerror(99,"Couldn't find entry for VCF file %s in master locus file",fn);
		return 0;
	}
#else
	for (i=0;i<nLocusFiles;++i)
	{
		if (!strcmp(lfFileNames[i],fn))
		{
			currentLocusFile=i;
			break;
		}
	}
	if (i==nLocusFiles)
	{
		dcerror(99,"Couldn't find entry for locus data file %s in master locus file\n",fn);
		return 0;
	}
#endif
	if (locusFiles[currentLocusFile]->fp!=0)
		fclose(locusFiles[currentLocusFile]->fp);
	if ((locusFiles[currentLocusFile]->fp=fopen(fn,"rb"))==0)
	{
		dcerror(99,"Could not open locus data file %s\n",fn);
		return 0;
	}
	cc[currentLocusFile]=aff;
	locusFiles[currentLocusFile]->readHeaderInfo();
	nSubs[currentLocusFile]=locusFiles[currentLocusFile]->getNSUbs();
	nread=0;
	if (tempLocus!=0)
		delete tempLocus;
	if (tempRecord.myLocalLocus[currentLocusFile]!=0)
		delete tempRecord.myLocalLocus[currentLocusFile];
	switch (fileTypes[currentLocusFile])
	{
	case VCFFILE:
		tempLocus=(LOCALLOCUSPTR) new vcfLocalLocus;
		tempRecord.myLocalLocus[currentLocusFile] = (LOCALLOCUSPTR) new vcfLocalLocus;
		break;
	case SHAPEITHAPSFILE:
		tempLocus=(LOCALLOCUSPTR) new hapsLocalLocus;
		tempRecord.myLocalLocus[currentLocusFile] = (LOCALLOCUSPTR) new hapsLocalLocus;
		break;
	default:
		dcerror(99,"In masterLocusFile::readLocusFileEntries() with unimplemented locus data file type %d",fileTypes[currentLocusFile]);
		return 0;
		break;
	}
	while (addLocus(locusFiles[currentLocusFile]->fp,spec))
		++nread;
	delete tempLocus;
	tempLocus=0;
    fclose(locusFiles[currentLocusFile]->fp);
	locusFiles[currentLocusFile]->fp=0;
	return nread;
}

int masterLocusFile::openFiles(char *rfn,char *ifn)
{
	char line[100],*ptr;
	int i;
	if ((recordFile=fopen(rfn,"rb+"))==0)
	{
		if (recordFile=fopen(ifn,"rb+"))
		{
			fclose(recordFile);
			recordFile=0;
			dcerror(99,"Could open index file %s but not master locus file %s. Assume filename misspelt. Else delete index file to start from scratch",ifn,rfn);
			return 0;
		}
		if ((recordFile=fopen(rfn,"wb+"))==0)
		{
			dcerror(99,"Could not open or create master locus file %s",rfn);
			return 0;
		}
		fwrite("Master locus file",strlen("Master locus file"),1,recordFile); // just so first record not at 0
		fwrite(&nLocusFiles,sizeof(nLocusFiles),1,recordFile);
		for (i=0;i<nLocusFiles;++i)
		{
			lfFileNames[i][0]='\0';
			fwrite(lfFileNames[i],sizeof(lfFileNames[i]),1,recordFile);
		}
		for (i = 0; i < nLocusFiles; ++i)
		{
			fileTypes[i]=NOFILETYPE;
			fwrite(&fileTypes[i],sizeof(fileTypes[i]),1,recordFile);
		}
		for (i=0;i<nLocusFiles;++i)
		{
			nSubs[i]=0;
			fwrite(&nSubs[i],sizeof(nSubs[i]),1,recordFile);
		}
		for (i=0;i<nLocusFiles;++i)
		{
			cc[i]=-1;
			fwrite(&cc[i],sizeof(cc[i]),1,recordFile);
		}
		FSEEK(recordFile,0L,SEEK_SET);
		if (!index.make_new(ifn))
		{
			dcerror(99,"Could not create new index file %s",ifn);
			return 0;
		}
		index.add("NO RECORD",0L); // just something to get index started so it is not empty
	}
	else
	{
	fread(line,strlen("Master locus file"),1,recordFile);
	if (strncmp(line,"Master locus file",strlen("Master locus file")))
	{
		dcerror(99,"Master locus file %s does not start with correct characters",rfn);
		return 0;
	}
	fread(&i,sizeof(i),1,recordFile);
	if (i!=nLocusFiles)
	{
		dcerror(99,"Master record file %s contains wrong number of VCF filenames. Expecting %d but found %d.",
			rfn,nLocusFiles,i);
		return 0;
	}
	for (i=0;i<nLocusFiles;++i)
		{
			fread(lfFileNames[i],sizeof(lfFileNames[i]),1,recordFile);
		}
	for (i = 0; i < nLocusFiles; ++i)
		{
			fread(&fileTypes[i],sizeof(fileTypes[i]),1,recordFile);
		}
	for (i=0;i<nLocusFiles;++i)
		fread(&nSubs[i],sizeof(nSubs[i]),1,recordFile);
	for (i=0;i<nLocusFiles;++i)
		fread(&cc[i],sizeof(cc[i]),1,recordFile);
	if (!index.open_old(ifn))
	{
		dcerror(99,"Could open master locus file %s but not index file %s. Assume filename misspelt. Else delete master locus file to start from scratch",rfn,ifn);
		return 0;
	}
	}
	return 1;
}

int masterLocusFile::fill(masterLocus &rec,localLocus *loc,FILEPOSITION locusPosInFile)
{
	char line[MAXALLLENGTH*(MAXALL+1)+1],rest[MAXALLLENGTH*MAXALL+1],*ptr;
	int i,a;
	rec.PolyPhen[0]='\0';
	rec.chr=loc->chr;
	rec.pos=loc->pos;
	rec.SNP=loc->isSNP();
	if (loc->id[0]=='.')
		rec.masterID[0]='\0'; // leave for maybe a matching locus in another file to provide rsID
	else
		strcpy(rec.masterID,loc->id);
	strcpy(rec.ref,loc->ref);
	strcpy(rec.alt,loc->alt);
	strcpy(line,loc->ref);
	if (loc->alt[0]!='.')
	{
		strcat(line,",");
		strcat(line,loc->alt);
		// now have comma-separated list of all alleles
	}
	for (i = 0; i < nLocusFiles; ++i)
	{
		for (a = 0; a<MAXALL; ++a)
			rec.alleleMapping[i][a] = -1;
	}
//	for (rec.nAlls=0,rest[0]='\0';sscanf(line,"%[^,],%s",rec.alls[rec.nAlls],rest)>=1;++rec.nAlls)
	for (rec.nAlls=0,ptr=line;*ptr!='\0';++rec.nAlls)
	{
		if (rec.nAlls>MAXALL)
		{
			dcerror(4,"This variant has too many alleles (max %d): %s\n",MAXALL,line);
			return 0;
		}
		scanWord(&ptr,rec.alls[rec.nAlls],MAXALLLENGTH-1,',');
		if (strlen(rec.alls[rec.nAlls])>1)
			rec.SNP=SNP_NO;
		rec.alleleMapping[currentLocusFile][rec.nAlls]=rec.nAlls;
//		strcpy(line,rest);
//		rest[0]='\0';
	}
	if (rec.nAlls>2)
		rec.SNP=SNP_NO;
	rec.locusPosInFile[currentLocusFile]=locusPosInFile;
	rec.myLocalLocus[currentLocusFile]->typeSpecificCopy(loc);
	// assume that loc is same derived class as rec.myLocalLocus[currentLocusFile]
	// we allocated a localLocus of the right type to deal with the file we are currently reading information from
	if (loc->PolyPhen[0]!='\0')
		strcpy(rec.PolyPhen,loc->PolyPhen);
	for (i = 0; i < nLocusFiles; ++i)
		if (i != currentLocusFile)
		{
			rec.locusPosInFile[i] = 0L; // make it clear that other VCF files may not have this locus
			strcpy(rec.myLocalLocus[i]->filter, "UNTYPED");
		}
	return 1;
}

int masterLocusFile::merge(masterLocus &rec,localLocus *loc,FILEPOSITION locusPosInFile)
{
	char line[MAXALLLENGTH*2],rest[MAXALLLENGTH*2],all[MAXALLLENGTH*2];
	int l,a;
	if (rec.masterID[0]=='\0' && loc->id[0]!='\0' && loc->id[0]!='.')
		strcpy(rec.masterID,loc->id);
	if (rec.SNP==SNP_MAYBE)
		rec.SNP=loc->isSNP();
	strcpy(line,loc->ref);
	if (loc->alt[0]!='.' || loc->alt[1]!='\0')
	{
		strcat(line,",");
		strcat(line,loc->alt);
		// now have comma-separated list of all alleles
	}
	for (a=0,rest[0]='\0';sscanf(line,"%[^,],%s",all,rest)>=1;++a)
	{
		if (strlen(all)>1)
			rec.SNP=SNP_NO;
		for (l=0;l<rec.nAlls;++l)
			if (!strcmp(all,rec.alls[l]))
				break;
		if(strlen(all)>MAXALLLENGTH)
		{
			dcerror(1,"Length of allele for locus at %d:%ld exceeds MAXALLLENGTH of %d so cannot merge(). Need to increase MAXALLLENGTH.\n%s",rec.chr,rec.pos,MAXALLLENGTH,all);
			return 0;
		}
		if (l==rec.nAlls)
			if (rec.nAlls>=MAXALL)
			{
				dcerror(1,"Total number of alleles for locus at %d:%ld exceeds MAXALL of %d so cannot merge(). Need to increase MAXALL.",rec.chr,rec.pos,MAXALL);
				return 0;
			}
			else
				strcpy(rec.alls[rec.nAlls++],all);

		rec.alleleMapping[currentLocusFile][a]=l;
		strcpy(line,rest);
		rest[0]='\0';
	}
	if (rec.nAlls>2)
		rec.SNP=SNP_NO;
	rec.locusPosInFile[currentLocusFile]=locusPosInFile;
	rec.myLocalLocus[currentLocusFile]->typeSpecificCopy(loc);
	if (rec.PolyPhen[0]=='\0' && loc->PolyPhen[0]!='\0')
		strcpy(rec.PolyPhen,loc->PolyPhen);
	return 1;
}

int masterLocusFile::load(masterLocus &rec,FILEPOSITION pos)
{
FSEEK(recordFile,pos,SEEK_SET);
return rec.read(recordFile);
}

int masterLocusFile::save(masterLocus &rec,FILEPOSITION pos)
{
FSEEK(recordFile,pos,SEEK_SET);
return rec.write(recordFile);
}

locusSNP localLocus::isSNP()
{
	// allow alleles to be coded as bases or 12
	char *ptr;
	if (SNP!=SNP_MAYBE)
		; // already sorted
	else if (strlen(ref)!=1 || strchr("CGAT12",ref[0])==0)
		SNP=SNP_NO;
	else
		if (strlen(alt)==1)
			if (alt[0]=='.')
				SNP=SNP_MAYBE;
			else if (strchr("CGAT12",alt[0]))
				SNP=SNP_YES;
			else SNP=SNP_NO;
		else
		{
			ptr=alt;
			while (1)
			{
				if (strchr("CGAT12",*ptr)==0)
				{
					SNP=SNP_NO;
					break;
				}
				else if (ptr[1]=='\0')
				{
					SNP=SNP_YES;
					break;
				}
				else if (ptr[1]!=',')
				{
					SNP=SNP_NO;
					break;
				}
				else
					ptr+=2;
			}
		}
return SNP;
}

localLocus::localLocus()
{
	clear();
}

void localLocus::clear()
{
	int i;
	SNP=SNP_MAYBE;
	nAltAlls=0;
	for (i=0;i<MAXALL;++i)
		alleleFreq[i]=0;
	AF=0;
	strcpy(filter,"UNTYPED");
	id[0]=ref[0]=alt[0]=PolyPhen[0]='\0';
}

void vcfLocalLocus::clear()
{
	localLocus::clear();
	qual=0;
	GQpos=-1;
	GTpos=-1;
	GPpos=-1;
	ADpos=-1;
}

int masterLocusFile::setCurrentFile(char *fn)
{
	for (currentLocusFile=0;currentLocusFile<nLocusFiles;++currentLocusFile)
		if (!strcmp(fn,lfFileNames[currentLocusFile]))
			break;
	if (currentLocusFile==nLocusFiles)
	{
		currentLocusFile=-1;
		dcerror(99,"Could not find entry for VCF file: %s",fn);
		return 0;
	}
	else 
		return 1;
}

masterLocusFile::masterLocusFile(int nLF) : tempRecord(nLF)
{
int i;
nLocusFiles=nLF;
lfFileNames=new filenamestring[nLF];
locusFiles=new LOCUSFILEPTR[nLF];
fileTypes=new locusFileType[nLF];
nSubs=new int[nLF];
cc=new int[nLF];
holdsFreqs=new int[nLF];
for (i = 0; i < nLF; ++i)
{
	locusFiles[i] = 0;
	holdsFreqs[i]=0;
	lfFileNames[i][0]='\0';
}
recordFile=0;
if (locusFile::buff==0)
		locusFile::buff=new char[BUFFSIZE]; // This never gets deleted. Oh well.
tempLocus=0;
}

masterLocusFile::~masterLocusFile()
{
	int i;
	if (recordFile)
	{
		FSEEK(recordFile,0L,SEEK_SET);
		fwrite("Master locus file",strlen("Master locus file"),1,recordFile); // just so first record not at 0
		fwrite(&nLocusFiles,sizeof(nLocusFiles),1,recordFile);
		for (i=0;i<nLocusFiles;++i)
			fwrite(lfFileNames[i],sizeof(lfFileNames[i]),1,recordFile);
		for (i = 0; i < nLocusFiles; ++i)
			fwrite(&fileTypes[i],sizeof(fileTypes[i]),1,recordFile);
		for (i=0;i<nLocusFiles;++i)
			fwrite(&nSubs[i],sizeof(nSubs[i]),1,recordFile);
		for (i=0;i<nLocusFiles;++i)
			fwrite(&cc[i],sizeof(cc[i]),1,recordFile);
		fclose(recordFile);
	}
	closeLocusFiles();
	for (i=0;i<nLocusFiles;++i)
		if (locusFiles[i]!=0)
			delete locusFiles[i];
	delete[] lfFileNames;
	delete[] nSubs;
	delete[] cc;
	delete[] locusFiles;
	delete[] fileTypes;
	delete[] holdsFreqs;
	if (tempLocus!=0)
		delete tempLocus;
}

int masterLocus::read(FILE *fp)
// this will break if elements changed
{
	int i;
	fread(&chr,sizeof(chr),1,fp);
	fread(&pos,sizeof(pos),1,fp);
	fread(&SNP,sizeof(SNP),1,fp);
	fread(&masterID,sizeof(masterID),1,fp);
	fread(&nLocusFiles,sizeof(nLocusFiles),1,fp);
	fread(&ref,sizeof(ref),1,fp);
	fread(&alt,sizeof(alt),1,fp);
	fread(&nAlls,sizeof(nAlls),1,fp);
	fread(&alls,sizeof(alls),1,fp);
	fread(ensemblConsequence,sizeof(ensemblConsequence),1,fp);
	fread(quickConsequence,sizeof(quickConsequence),1,fp);
	fread(&worstConsequenceType,sizeof(worstConsequenceType),1,fp);
	fread(PolyPhen,sizeof(PolyPhen),1,fp);
	for (i=0;i<nLocusFiles;++i)
		fread(&alleleMapping[i],sizeof(alleleMapping[i]),1,fp);
	for (i=0;i<nLocusFiles;++i)
		fread(&locusPosInFile[i],sizeof(locusPosInFile[i]),1,fp);
	for (i=0;i<nLocusFiles;++i)
		myLocalLocus[i]->read(fp);
	return 1;
}

int masterLocus::write(FILE *fp)
// this will break if elements changed
{
	int i;
	fwrite(&chr,sizeof(chr),1,fp);
	fwrite(&pos,sizeof(pos),1,fp);
	fwrite(&SNP,sizeof(SNP),1,fp);
	fwrite(&masterID,sizeof(masterID),1,fp);
	fwrite(&nLocusFiles,sizeof(nLocusFiles),1,fp);
	fwrite(&ref,sizeof(ref),1,fp);
	fwrite(&alt,sizeof(alt),1,fp);
	fwrite(&nAlls,sizeof(nAlls),1,fp);
	fwrite(&alls,sizeof(alls),1,fp);
	fwrite(ensemblConsequence,sizeof(ensemblConsequence),1,fp);
	fwrite(quickConsequence,sizeof(quickConsequence),1,fp);
	fwrite(&worstConsequenceType,sizeof(worstConsequenceType),1,fp);
	fwrite(PolyPhen,sizeof(PolyPhen),1,fp);
	for (i=0;i<nLocusFiles;++i)
		fwrite(&alleleMapping[i],sizeof(alleleMapping[i]),1,fp);
	for (i=0;i<nLocusFiles;++i)
		fwrite(&locusPosInFile[i],sizeof(locusPosInFile[i]),1,fp);
	for (i=0;i<nLocusFiles;++i)
		myLocalLocus[i]->write(fp);
	return 1;
}

masterLocus::masterLocus(int nLF)
{
	int i;
	nLocusFiles=nLF;
alleleMapping=new alleleMap[nLF];
locusPosInFile=new FILEPOSITION[nLF];
myLocalLocus=new LOCALLOCUSPTR[nLF];
for (i=0;i<nLF;++i)
	myLocalLocus[i]=(LOCALLOCUSPTR) new vcfLocalLocus; 
// one day will have to allocate localLocus of correct type to these
ensemblConsequence[0]=quickConsequence[0]=PolyPhen[0]='\0';
worstConsequenceType=NULL_CONSEQUENCE;
masterID[0]='\0';
genoCount[0] = genoCount[1] = genoCount[2] = 0;
}

masterLocus::~masterLocus()
{
int i;
for (i=0;i<nLocusFiles;++i)
	if (myLocalLocus[i]!=0)
		delete myLocalLocus[i]; 
delete[] alleleMapping;
	delete[] locusPosInFile;
	delete[] myLocalLocus;
}

FILEPOSITION masterLocusFile::findFirstInRange(analysisSpecs const &spec)
{
char key[1000];
const char *testKey;
FILEPOSITION recPos,p;
int c;
sprintf(key,"%3d %10ld",spec.sc,spec.sp);
recPos=index.exact_find(key);
if (recPos==0L)
	recPos=index.near_find(key);
if (recPos!=0L)
// go back from here and try to find earliest which will match
{
	while (1)
	{
		recPos=index.get_prev();
		if (recPos==0L)
		{
			recPos=index.get_first();
			break;
		}
		testKey=index.current_key();
		if ((c=atoi(testKey))<spec.sc || (p=atol(testKey+3))<spec.sp)
		{
			recPos=index.get_next();
			break;
		}
	}
	// OK, by now we have gone back too far then come forward one
	// so this one may match but need to check and if not go forward
	while (1)
	{
		testKey=index.current_key();
		if ((c=atoi(testKey))==0)
		{
			recPos=0L;
			break;
			// looks like there were no valid entries for this chromosome
		}
		else if (c>spec.sc || (c==spec.sc && (p=atol(testKey+3))>=spec.sp))
			break;
		// OK, this is the first valid one
		else
		{
			recPos=index.get_next();
			if (recPos==0L)
				break;
			// got to end of index - I think this will never happen 
		}
	}
}
return recPos;
// at this point recPos is 0L or first valid record
// and if it is valid index.getKey() still has key for this record so can do other checks
}

int masterLocusFile::addLocus(FILE *f,analysisSpecs const &spec)
{
FILEPOSITION recPos,locusPosInFile,p;
char key[MAXALLLENGTH*2+100],refTemp[MAXALLLENGTH+1],altTemp[MAXALLLENGTH+1],*ptr;
const char *testKey;
// locusPosInFile=FTELL(f); get this set by input() below
if (!tempLocus->input(f,&locusPosInFile,spec))
	return 0;
locusSNP gotSNP=tempLocus->isSNP();
recPos=findFirstInRange(tempLocus->chr,tempLocus->pos);
if (recPos!=0L)
{
	if (atol(index.current_key()+3)!=tempLocus->pos)
		recPos=0L;
	else
	// we possibly have an old master locus we can use 
    // now we have to see if it is an exact match
	while (1)
	{
		testKey=index.current_key();
		if ((p=atol(testKey+3))!=tempLocus->pos)
		{
			recPos=0L;
			break;
		}
		FSEEK(recordFile,recPos,SEEK_SET);
		tempRecord.read(recordFile);
		if (spec.ignoreAlleles) // treat any variant at same position as if identical
			break;
		if (gotSNP!=SNP_NO)
		{
			if (tempRecord.isSNP()!=SNP_NO)
				break;
			// both look like they could be SNP at this position
			// NB it looks like we are treating all SNPs as the same even if alt allele is different
		}
		else
		{
			if (!strcmp(tempRecord.ref,tempLocus->ref) && !strcmp(tempRecord.alt,tempLocus->alt))
				break;
			// other polymorphism which is exact match
		}
		recPos=index.get_next();
		if (recPos==0L)
			break;
	}
}
// if we have no match then recPos is 0L, else we have read matching record into tempRecord
if (recPos!=0L)
{
	merge(tempRecord,tempLocus,locusPosInFile);
}
else
{
	fill(tempRecord,tempLocus,locusPosInFile);
	FSEEK(recordFile,0L,SEEK_END);
	recPos=FTELL(recordFile);
	sprintf(key,"%3d %10ld",tempLocus->chr,tempLocus->pos);
	index.add(key,recPos);
}
save(tempRecord,recPos);
ptr=tempLocus->ref;
scanWord(&ptr,refTemp,MAXALLLENGTH);
ptr=tempLocus->alt;
scanWord(&ptr,altTemp,MAXALLLENGTH);
sprintf(key,"%s %3d %10ld %s %s",
		lfFileNames[currentLocusFile],tempLocus->chr,tempLocus->pos,refTemp,altTemp);
// now add a key to the same master record but specific for this file
index.add(key,recPos);
// in fact, BT_KSIZ is only 40 so it looks as if key will be truncated to that many characters
return 1;
}


