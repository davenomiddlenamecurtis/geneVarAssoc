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

#if 0
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
#endif

#if 0
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
#endif

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

int masterLocusFile::getTotalSubs   ()
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

int masterLocus::readQueryOutput(FILE *fp, analysisSpecs const &spec )
{
	char line[200],CDS_positionStr[100];
	int all;
	if (alls[0][0]=='N')
		return 1; // predictor produces no output for these
	if (nAlls<2)
		return 1;
	for (all = 0; all < (spec.mergeAltAlleles ? 1 : nAlls - 1); ++all)
	{
		if (!fgets(line, 199, fp))
		{
			dcerror(1, "Could not read from variant_predictor output");
			return 0;
		}
		CDS_positionStr[0] = '-';
		CDS_positionStr[0] = '\0';
		if (sscanf(line, "%*s %*s %*s %*s %*s %*s %s %*s %s", ensemblConsequence[all+1], CDS_positionStr) < 1)
		{
			dcerror(1, "Could not find consequence in this line from variant_predictor: \n%s", line);
			return 0;
		}
		if (CDS_positionStr[0] != '-')
			sprintf(strchr(ensemblConsequence[all+1], '\0'), " %s", CDS_positionStr);
	}
	return 1;
}

int masterLocus::writePredictorQuery(FILE *fp, analysisSpecs const &spec )
{
	long start,end;
	char chrStr[20];
	int all;
	if (alls[0][0]=='N')
		return 1; // predictor produces no output for these
	if (nAlls<2)
		return 1;
	for (all = 0; all < (spec.mergeAltAlleles ? 1 : nAlls - 1);++all)
	{
		start = pos;
			end = pos + strlen(alls[all+1]) - 1; // does this always work ?
			if (chr == 23)
				strcpy(chrStr, "X");
			else
				sprintf(chrStr, "%d", chr);
			fprintf(fp, "%s\t%ld\t%ld\t%s/%s\t+\n", chrStr, start, end, (alls[0][0]=='*'?"-":alls[0]), (alls[0][all + 1] == '*' ? "-" : alls[all + 1]));
	}
	return 1;
}

int masterLocus::getQuickFeature(refseqGeneInfo &r,int a)
{
	char thisGeneEffect[100];
	double oldEffectWeight;
	int e;
	if (nAlls<2)
		return 1;
	if (chr==r.getChrNum())
	{
		r.getEffect(pos,alls[0],alls[a]);
		if (quickConsequence[a][0] == '\0' || !strncmp(quickConsequence[a], "NULL", 4) || !strncmp(quickConsequence[a], "INTER", 5))
		{
			strcpy(quickConsequence[a], r.tellEffect());
			worstConsequenceType[a]=r.tellWorstConsequence();
		}
		else
		{
			strcpy(thisGeneEffect,r.tellEffect());
			for (e=0;e<NCONSEQUENCETYPES;++e)
				if (!strncmp(quickConsequence[a],consequence[e].str,strlen(consequence[e].str)))
				{
					oldEffectWeight=consequence[e].weight;
					break;
				}
				if (oldEffectWeight < consequence[r.tellWorstConsequence()].weight)
				{
					strcpy(quickConsequence[a], r.tellEffect());
					worstConsequenceType[a]=r.tellWorstConsequence();
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
	int c,all;
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
		for (all = 0; all < (spec.mergeAltAlleles ? 1 : tempRecord.nAlls - 1); ++all)
		{
			if (!redo) // may have already been set using a different gene
				tempRecord.quickConsequence[all+1][0] = '\0';
			tempRecord.getQuickFeature(r,all+1);
		}
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
		tempRecord.writePredictorQuery(fp,spec);
		++locusCount;
		recPos=index.get_next();
		if (recPos==0L)
			break;

	}
	fclose(fp);
	unlink("predictorOutput.txt");
	sprintf(line," %s -i predictorQuery.txt -o predictorOutput.txt --pick_allele_gene --force_overwrite",spec.vepCommand);
	checkSystem();
	system(line);
	fp=fopen("predictorOutput.txt","rb"); // binary mode can use fseek/ftell
	if (fp==NULL)
	{
		dcerror.kill();
		dcerror(1,"Could not open output file predictorOutput.txt from variant_effect_predictor.pl");
		return 0;
	}
	do {
		fPos=FTELL(fp);
		if (!fgets(line,999,fp))
		{
			dcerror.kill();
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
		tempRecord.readQueryOutput(fp, spec);
		save(tempRecord,recPos);
		recPos=index.get_next();
		if (recPos==0L)
			break;
	}
	fclose(fp);
}
return locusCount;
}

// #define MAXLOCIINSCOREASSOCFILE 50000
int useLocus[MAXLOCIINSCOREASSOCFILE];

int masterLocusFile::writeScoreAssocFiles(char *root, float wf, int *useFreqs, int *suppliedNSubs, int writeNames, int writeComments, int writeScoreFile, int writeRecScoreFile, analysisSpecs &spec)
{
	return writeScoreAssocFiles(*this,root,wf,useFreqs,suppliedNSubs,writeNames,writeComments,writeScoreFile, writeRecScoreFile, spec);
}

int masterLocusFile::loadFirst(analysisSpecs &spec)
{
	currentRecPos=findFirstInRange(spec);
	if (currentRecPos==0L)
		return 0;
	else 
		return load(tempRecord,currentRecPos);
}

int masterLocusFile::writeFlatFile(masterLocusFile& subFile, char* fn, int totalSub, strEntry* subName,analysisSpecs& spec,int *useLocus)
{
	int lc,s,ss,i,all,l,ll,lll;
	FILE* fp;
	long* subPos;
	allelePair* a;
	probTriple* p;
	assert(subPos = (long *)calloc(totalSub, sizeof(long)));
	lc= 0;
	if (gotoFirstInRange(spec))
		do
		{
			nAlls[lc++] = tempRecord.nAlls;
		} while (gotoNextInRange(spec));

	if (spec.useProbs)
		assert(p= (probTriple*)calloc(totalSub, sizeof(probTriple)));
	else
		assert(a= (allelePair*)calloc(totalSub, sizeof(allelePair)));
	fp = fopen(fn, "wb"); // binary because I had problems with fseek() and text files
	for (s = 0, i = 0; i < subFile.nLocusFiles; ++i)
		for (ss = 0; ss < subFile.nSubs[i]; ++s, ++ss)
		{
			if (spec.phenotypes)
			{
				if (spec.phenotypes[s] == MISSINGPHENOTYPE)
					continue;
			}
			if (spec.isQuantitative)
				fprintf(fp, "%s\t%9.5f\t", subName[s], spec.phenotypes[s]);
			else
			{
				int cc_pheno = spec.phenotypes ? spec.phenotypes[s] : subFile.cc[i];
				if (cc_pheno != 0 && cc_pheno != 1)
					continue; // no longer output subjects with unknown phenotype
				fprintf(fp, "%s\t%d\t", subName[s], cc_pheno);
			}
			subPos[s] = ftell(fp);
			for (l = 0,ll=0; l < lc; ++l)
			{
				if (spec.useProbs)
				{
					if (useLocus[ll++])
						fprintf(fp, "0.000 0.000 0.000\t");
//						fprintf(fp, "%5.3f %5.3f %5.3f\t", 0, 0, 0);
				}
				else
				{
					if (spec.mergeAltAlleles)
					{
						if (useLocus[ll++])
							fprintf(fp, "0 0\t");
//							fprintf(fp, "%d %d\t", 0, 0);
					}
					else
						for (all = 1; all < nAlls[l]; ++all)
						{
							if (useLocus[ll++])
								fprintf(fp, "0 0\t");
//								fprintf(fp, "%d %d\t", 0, 0);
						}
				}
			}
			fprintf(fp, "\n");
		}
// subPos[s] is the position to write the next genotype / probs
	l = 0;
	ll = 0;
	if (gotoFirstInRange(spec))
		do
		{
			if (spec.useProbs)
				outputCurrentProbs(p, spec);
			else
				outputCurrentAlleles(a, spec);
			for (s = 0, i = 0; i < subFile.nLocusFiles; ++i)
				for (ss = 0; ss < subFile.nSubs[i]; ++s, ++ss)
				{
					if (spec.phenotypes)
					{
						if (spec.phenotypes[s] == MISSINGPHENOTYPE)
							continue;
					}
					if (!spec.isQuantitative)
					{
						int cc_pheno = spec.phenotypes ? spec.phenotypes[s] : subFile.cc[i];
						if (cc_pheno != 0 && cc_pheno != 1)
							continue;
					}
					fseek(fp, subPos[s], SEEK_SET);
					if (spec.useProbs)
					{
						if (useLocus[ll])
							fprintf(fp, "%5.3f %5.3f %5.3f\t", p[s][0], p[s][1], p[s][2]);
					}
					else
					{
							if (spec.mergeAltAlleles)
							{
								if (useLocus[ll])
									fprintf(fp, "%d %d\t",
										(a[s][0] > 1) ? 2 : a[s][0],
										(a[s][1] > 1) ? 2 : a[s][1]); // force to be biallelic
							}
							else
								for (all = 1,lll=ll; all < nAlls[l]; ++all)
								{
									if (useLocus[lll++])
										fprintf(fp, "%d %d\t",
										(a[s][0] == all + 1) ? 2 : a[s][0] == 0 ? 0 : 1,
										(a[s][1] == all + 1) ? 2 : a[s][1] == 0 ? 0 : 1);
								}
						}
						subPos[s] = ftell(fp);
				}
			ll += (spec.useProbs||spec.mergeAltAlleles) ? 1 : (nAlls[l] - 1);
			++l;
		} while (gotoNextInRange(spec));
	fclose(fp);
	if (spec.useProbs)
		free(p);
	else
		free(a);
	free(subPos);
	return lc;
}

int masterLocusFile::writeTransposedFile(masterLocusFile& subFile, char* fn, int totalSub, strEntry* subName, analysisSpecs& spec, int* useLocus)
{
	int lc, s, ss, i, all, l, ll, lll;
	FILE* fp;
	allelePair* a;
	probTriple* p;
	int *useSub;
	assert((useSub = (int*)calloc(totalSub, sizeof(int)))!=0);
	lc = 0;
	if (gotoFirstInRange(spec))
		do
		{
			nAlls[lc++] = tempRecord.nAlls;
		} while (gotoNextInRange(spec));

		if (spec.useProbs)
			assert(p = (probTriple*)calloc(totalSub, sizeof(probTriple)));
		else
			assert(a = (allelePair*)calloc(totalSub, sizeof(allelePair)));
		fp = fopen(fn, "w");
		for (s = 0, i = 0; i < subFile.nLocusFiles; ++i)
			for (ss = 0; ss < subFile.nSubs[i]; ++s, ++ss)
			{
				if (spec.phenotypes && spec.isQuantitative)
				{
					if (spec.phenotypes[s] == MISSINGPHENOTYPE)
						continue;
				}
				else
				{
					int cc_pheno = spec.phenotypes ? spec.phenotypes[s] : subFile.cc[i];
					if (cc_pheno != 0 && cc_pheno != 1)
						continue; // no longer output subjects with unknown phenotype
				}
				fprintf(fp, "%s\t", subName[s]);
				useSub[s] = 1;
			}
		fprintf(fp, "\n");
		for (s = 0, i = 0; i < subFile.nLocusFiles; ++i)
			for (ss = 0; ss < subFile.nSubs[i]; ++s, ++ss)
			{
				if (!useSub[s])
					continue;
				if (spec.isQuantitative)
					fprintf(fp, "%9.5f\t", spec.phenotypes[s]);
				else
				{
					int cc_pheno = spec.phenotypes ? spec.phenotypes[s] : subFile.cc[i];
					fprintf(fp, "%d\t", cc_pheno);
				}
			}
		fprintf(fp, "\n");


		l = 0;
		ll = 0;
		if (gotoFirstInRange(spec))
			do
			{
				if (spec.useProbs)
					outputCurrentProbs(p, spec);
				else
					outputCurrentAlleles(a, spec);
				for (all = 1; all < ((spec.useProbs || spec.mergeAltAlleles) ? 1 : nAlls[l]); ++all)
				{
					if (useLocus[ll + all - 1])
					{
						for (s = 0, i = 0; i < subFile.nLocusFiles; ++i)
							for (ss = 0; ss < subFile.nSubs[i]; ++s, ++ss)
							{
								if (!useSub[s])
									continue;
								if (spec.useProbs)
									fprintf(fp, "%5.3f %5.3f %5.3f\t", p[s][0], p[s][1], p[s][2]);
								else if (spec.mergeAltAlleles)
									fprintf(fp, "%d %d\t",
										(a[s][0] > 1) ? 2 : a[s][0],
										(a[s][1] > 1) ? 2 : a[s][1]); // force to be biallelic
								else
									fprintf(fp, "%d %d\t",
										(a[s][0] == all + 1) ? 2 : a[s][0] == 0 ? 0 : 1,
										(a[s][1] == all + 1) ? 2 : a[s][1] == 0 ? 0 : 1);
							}
						fprintf(fp, "\n");
					}
				}
				ll += (spec.useProbs || spec.mergeAltAlleles) ? 1 : (nAlls[l] - 1);
				++l;
			} while (gotoNextInRange(spec));
			fclose(fp);
			if (spec.useProbs)
				free(p);
			else
				free(a);
			free(useSub);
			return lc;
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
#define MAXCOMMENTLENGTH 250000 // can be 200K if multiple transcripts
char comment[MAXCOMMENTLENGTH];
// may be long VEP output
int masterLocusFile::writeScoreAssocFiles(masterLocusFile &subFile,char *root, float wf, int *useFreqs, int *suppliedNSubs, int writeNames, int writeComments, int writeScorefile, int writeRecScorefile, analysisSpecs &spec)
// allow information about subjects to be provided by a different masterLocusFile 
// however we are assuming both files refer to identical set of subjects
{
	char fn[100],buff[MAXALL*MAXALLLENGTH],buff2[1000],*ptr,alleles[MAXSTR+1],commandString[5000],posStr[100];
	allelePair **a;
	probTriple **p;
	int totalSub,lc,s,l,ss,i,c,nValid,all,numSplitLoci,numIncludedLoci,ll,numWeights,numRecWeights,w;
	float** locusWeights;
	FILE *fp;
	FILEPOSITION recPos;
	const char *testKey;
	geneVarParser commentParser;
	geneVarParser::mergeAltAlleles = spec.mergeAltAlleles; // need to make sure this gets set
	geneVarParser::multilineVEP = spec.multilineVEP; // need to make sure this gets set
	strEntry *subName;
	if (spec.commentExpression[0])
		commentParser.parse(spec.commentExpression); // only have to parse once
	if (nLocusFiles==1)
		nSubs[0]=subFile.getTotalSubs();
	else if (nLocusFiles==subFile.nLocusFiles)
		for (i=0;i<nLocusFiles;++i)
			nSubs[i]=subFile.nSubs[i];
	else
	{
		dcerror(99,"Incompatible numbers of subjects in subFile in masterLocusFile::writeScoreAssocFiles()");
		return 0;
	}
	openLocusFiles();
	for (i=0,totalSub=0;i<subFile.nLocusFiles;++i)
		totalSub+=subFile.nSubs[i];
	if (totalSub == 0)
	{
		dcerror(99, "The total number of valid subjects is zero in masterLocusFile::writeScoreAssocFiles()");
		return 0;
	}
	assert(subName=(strEntry *)calloc(totalSub,sizeof(strEntry)));
// hereOK();
	subFile.outputSubNames(subName,spec);
// hereOK();
	nValid=countNumberInRange(spec);
	if (!spec.useFlatFile && !spec.useTransposedFile)
	{
		if (spec.useProbs)
		{
			assert((p = (probTriple**)calloc(nValid, sizeof(probTriple*))) != 0);
			for (l = 0; l < nValid; ++l)
				assert(p[l] = (probTriple*)calloc(totalSub, sizeof(probTriple)));
			lc = outputProbs(p, spec);
		}
		else
		{
			assert((a = (allelePair**)calloc(nValid, sizeof(allelePair*))) != 0);
			for (l = 0; l < nValid; ++l)
				assert(a[l] = (allelePair*)calloc(totalSub, sizeof(allelePair)));
			lc = outputAlleles(a, spec);
		}
	}
// hereOK();
	if(spec.subPhenos.size()>0)
	{
		if(spec.phenotypes==NULL) // should have been allocated when IDsAndPhenotypesFileName read
		{
			spec.phenotypes=(float*)malloc(sizeof(float)*totalSub);
			assert(spec.phenotypes!=0);
		}
		TStrFloatMap::iterator it;
		for(s=0;s<totalSub;++s)
		{
			it=spec.subPhenos.find(subName[s]);
			if (it==spec.subPhenos.end())
				spec.phenotypes[s]=MISSINGPHENOTYPE;
			else
				spec.phenotypes[s]=it->second;
		}
	}
	numWeights = spec.weightExpressions.size();
	numRecWeights = spec.recWeightExpressions.size();
	if (numWeights + numRecWeights == 0)
		numWeights = 1;
	assert((locusWeights = (float**)calloc(numWeights+numRecWeights, sizeof(float*))) != 0);
	for (w = 0; w < numWeights+numRecWeights; ++w)
		assert((locusWeights[w] = (float*)calloc(MAXLOCIINSCOREASSOCFILE, sizeof(float))) != 0);
	numSplitLoci = outputSAInfo(useLocus, locusWeights, spec);
	// useLocus has numSplitLoci entries
	checkSystem();
	sprintf(fn, "%s.dat", root);
	if (spec.useTransposedFile)
		lc = writeTransposedFile(subFile, fn, totalSub, subName, spec, useLocus);
	else if (spec.useFlatFile)
		lc=writeFlatFile(subFile, fn, totalSub, subName, spec,useLocus);
	else
	{
		fp = fopen(fn, "w");
		for (s = 0, i = 0; i < subFile.nLocusFiles; ++i)
			for (ss = 0; ss < subFile.nSubs[i]; ++s, ++ss)
			{
				if (spec.phenotypes)
				{
					if (spec.phenotypes[s] == MISSINGPHENOTYPE)
						continue;
				}
				if (spec.isQuantitative)
					fprintf(fp, "%s\t%9.5f\t", subName[s], spec.phenotypes[s]);
				else
				{
					int cc_pheno = spec.phenotypes ? spec.phenotypes[s] : subFile.cc[i];
					if (cc_pheno != 0 && cc_pheno != 1)
						continue; // no longer output subjects with unknown phenotype
					fprintf(fp, "%s\t%d\t", subName[s], cc_pheno);
				}
				for (l =0,ll=0; l< lc; ++l)
					if (spec.useProbs)
					{
						if (useLocus[ll++])
							fprintf(fp, "%5.3f %5.3f %5.3f\t", p[l][s][0], p[l][s][1], p[l][s][2]);
					}
					else
					{
						if (spec.mergeAltAlleles)
						{
							if (useLocus[ll++])
								fprintf(fp, "%d %d\t",
									(a[l][s][0] > 1) ? 2 : a[l][s][0],
									(a[l][s][1] > 1) ? 2 : a[l][s][1]); // force to be biallelic
						}
						else
							for (all = 1; all < nAlls[l]; ++all)
							{
								if (useLocus[ll++])
									fprintf(fp, "%d %d\t",
										(a[l][s][0] == all + 1) ? 2 : a[l][s][0] == 0 ? 0 : 1,
										(a[l][s][1] == all + 1) ? 2 : a[l][s][1] == 0 ? 0 : 1);
							}
					}
				fprintf(fp, "\n");
			}
		fclose(fp);
	}
	checkSystem();
#if 0
	if (!spec.mergeAltAlleles)
	{
		for (numSplitLoci = 0, l = 0; l < lc; ++l)
			numSplitLoci += nAlls[l] - 1;
	}
#endif
	numIncludedLoci = 0;
	for (ll = 0; ll < numSplitLoci; ++ll)
		if (useLocus[ll])
			++numIncludedLoci;
	sprintf(commandString,"scoreassoc %s %s --numloci %d",spec.useProbs?"--gendatafile":"--gcdatafile",
		fn, numIncludedLoci);
	sprintf(fn,"%s.lf.par",root);
	fp=fopen(fn,"w");
	for (l=0;l< numIncludedLoci;++l)
			fprintf(fp,"1\n");  // idea is now we only output valid loci
	fclose(fp);
	sprintf(strchr(commandString,'\0')," --locusfilterfile %s",fn);
	checkSystem();
#if 0
	if (spec.doRecessiveTest)
	{
		sprintf(strchr(commandString, '\0'), " --dorecessive 1 --minweight %f --ldthreshold %f ", spec.recWeightThreshold, spec.LDThreshold);
		if (spec.useHaplotypes)
			sprintf(strchr(commandString, '\0'), " --usehaps 1");
		if (spec.showHapLocusNames)
			sprintf(strchr(commandString, '\0'), " --showhaplocusnames 1");
	}
#endif
	checkSystem();
	if (spec.useConsequenceWeights)
	{
		for (w = 0; w < numWeights; ++w)
		{
			if (numWeights == 1)
				sprintf(fn, "%s.lw.par", root);
			else
				sprintf(fn, "%s.%d.lw.par", root, w);
			fp = fopen(fn, "w");
			for (l = 0; l < numSplitLoci; ++l)
				if (useLocus[l])
					fprintf(fp, "%8.5f\n", locusWeights[w][l]);
			fclose(fp);
			sprintf(strchr(commandString, '\0'), " --locusweightfile %s", fn);
		}
		for (w = 0; w < numRecWeights; ++w)
		{
			if (numRecWeights == 1)
				sprintf(fn, "%s.lrw.par", root);
			else
				sprintf(fn, "%s.%d.lrw.par", root, w);
			fp = fopen(fn, "w");
			for (l = 0; l < numSplitLoci; ++l)
				if (useLocus[l])
					fprintf(fp, "%8.5f\n", locusWeights[numWeights+w][l]);
			fclose(fp);
			sprintf(strchr(commandString, '\0'), " --locusrecweightfile %s", fn);
		}
	}
	checkSystem();
	for (i=0;i<2;++i)
		if (useFreqs[i])
		{
			float *freqs;
			if (!spec.mergeAltAlleles)
			{
				dcerror.kill();
				dcerror(1,"Must merge alt alleles with --merge-alt-alleles 1 if using allele frequencies rather than counts.\n");
			}
			freqs=new float[lc];
			// outputAltFrequencies(freqs,i,sc,sp,ec,ep);
			outputEurAltFrequencies(freqs,i,spec);
			sprintf(fn,"%s.%s.freq.par",root,i?"case":"cont");
			fp=fopen(fn,"w");
			for (l=0;l<lc;++l)
				if (useLocus[l])
					fprintf(fp,"%8.6f ",freqs[l]);
			fprintf(fp,"\n");
			for (l=0;l<lc;++l)
				if (useLocus[l])
					fprintf(fp,"%8d ",suppliedNSubs[i]);
			// assume that if AF available nSubs is unknown
			// change this later to use nSubs if available
			fprintf(fp,"\n");
			fclose(fp);
			sprintf(strchr(commandString,'\0')," --%sfreqfile %s",i?"case":"cont",fn);
			delete[] freqs;
		}
// no longer dealing with names seperately
	checkSystem();
	if (spec.weightNames.size())
	{
		sprintf(fn, "%s.wn.par", root);
		fp = fopen(fn, "w");
		for (std::list<std::string>::const_iterator it = spec.weightNames.begin(); it != spec.weightNames.end(); ++it)
			fprintf(fp, "%s\n", it->c_str());
		fclose(fp);
		sprintf(strchr(commandString, '\0'), " --locusweightnamefile %s", fn);
	}
	if (spec.recWeightNames.size())
	{
		sprintf(fn, "%s.rwn.par", root);
		fp = fopen(fn, "w");
		for (std::list<std::string>::const_iterator it = spec.recWeightNames.begin(); it != spec.recWeightNames.end(); ++it)
			fprintf(fp, "%s\n", it->c_str());
		fclose(fp);
		sprintf(strchr(commandString, '\0'), " --locusrecweightnamefile %s", fn);
	}
	if (writeComments)
	{
		sprintf(fn,"%s.comm.par",root);
		fp=fopen(fn,"w");
		recPos=findFirstInRange(spec);
		ll = 0;
	if (recPos!=0L)
		while (1)
		{
			testKey=index.current_key();
			if ((c=atoi(testKey))==0 || c>spec.ec)
				break;
			if (c==spec.ec && atol(testKey+3)>spec.ep)
				break;
			load(tempRecord,recPos);
			for (all = 0; all < (spec.mergeAltAlleles ? 1 : tempRecord.nAlls - 1); ++all)
			{
				int aa;
				if (useLocus[ll++] == 0)
					continue;
				if (spec.mergeAltAlleles || tempRecord.nAlls == 2)
				{
					sprintf(posStr, "%d:%ld-%s", tempRecord.chr, tempRecord.pos,tempRecord.alls[0]);
					for (aa = 1; aa < tempRecord.nAlls; ++aa)
						sprintf(strchr(posStr, '\0'), ",%s", tempRecord.alls[aa]);
				}
				else
					sprintf(posStr, "%d:%ld.%d-%s,%s", tempRecord.chr, tempRecord.pos,all+1, tempRecord.alls[0], tempRecord.alls[all+1]);
				if (spec.commentExpression[0])
				{
					geneVarParser::thisLocus = &tempRecord;
					geneVarParser::thisAltAllele = all + 1;
					dcexpr_val*rv = commentParser.eval();
					if (spec.mergeAltAlleles||tempRecord.nAlls==2)
						sprintf(comment, "%s:%s", posStr, (char*)(*rv));
					else
						sprintf(comment, "%s:%s", posStr, (char*)(*rv));
					delete rv;
				}
				else if (tempRecord.ensemblConsequence[all+1][0] != '\0')
					sprintf(comment, "%s:%s", posStr, tempRecord.ensemblConsequence[all+1]);
				else if (tempRecord.quickConsequence[all+1][0] != '\0')
					sprintf(comment, "%s:%s", posStr, tempRecord.quickConsequence[all+1]);
				else
					sprintf(comment, "%s:", posStr);
				for (ptr = comment; *ptr; ++ptr)
					if (isspace(*ptr))
						*ptr = '_';
				strcat(comment, "-");
				buff[0] = '\0';
				for (i = 0; i < tempRecord.nAlls; ++i)
					if (spec.mergeAltAlleles)
						sprintf(strchr(buff, '\0'), "%s%s", ((i ==0) ? "" : "/"),tempRecord.alls[i] );
					else
					{
						if (i==0 || i==all+1)
							sprintf(strchr(buff, '\0'), "%s%s", ((i == 0) ? "" : "/"), tempRecord.alls[i] );
					}
				strncpy(alleles, buff, MAXSTR);
				alleles[MAXSTR] = '\0';
				strcat(comment, alleles);
				fprintf(fp, "%s\n", comment);
			}
			recPos=index.get_next();
			if (recPos==0L)
				break;
		}
	fclose(fp);
	sprintf(strchr(commandString,'\0')," --locusnamefile %s",fn);
	}

	if (spec.useTrios)
		sprintf(strchr(commandString,'\0')," --triofile %s",spec.triosFn);
	checkSystem();
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
		sprintf(strchr(commandString, '\0'), " --scorefile %s.sco", root);
	if (writeRecScorefile)
		sprintf(strchr(commandString, '\0'), " --recscorefile %s.recsco", root);
	sprintf(strchr(commandString, '\0'), " --weightfactor %f", spec.wf);
	if (spec.isQuantitative)
		sprintf(strchr(commandString, '\0'), " --isquantitative 1");
	if (spec.useTransposedFile)
		sprintf(strchr(commandString, '\0'), " --transposedata 1");
	for (l = 0; l < spec.nScoreassocArgs; ++l)
		sprintf(strchr(commandString, '\0'), " %s %s", spec.scoreassocArgs[l][0], spec.scoreassocArgs[l][1]);
	checkSystem();
#ifndef MSDOS
	sprintf(fn,"%s.sh",root);
#else
	sprintf(fn,"%s.bat",root);
#endif
	fp=fopen(fn,"w");
	fprintf(fp,"%s\n",commandString);
	fclose(fp);
	if (!spec.useFlatFile)
	{
		if (spec.useProbs)
		{
			for (l = 0; l < nValid; ++l)
				free(p[l]);
			free(p);
		}
		else
		{
			for (l = 0; l < nValid; ++l)
				free(a[l]);
			free(a);
		}
	}
	for (w = 0; w < numWeights; ++w)
		free(locusWeights[w]);
	free(locusWeights);
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

int masterLocusFile::outputSAInfo(int *useLocus,float **locusWeights,analysisSpecs const &spec)
{
	int locusCount,splitLocusCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c,i,doNotUseUntypedFreqfile,all,w,numWeights,numRecWeights;
	int cons;
	std::list<geneVarParser*> excludeParser, weightParser, recWeightParser;
	geneVarParser d;
	if (spec.debug)
		d.debug(stdout); // will set debug for all parsers
	splitLocusCount=locusCount=0; // keep track of split loci, when alt alleles are not merged
	recPos=findFirstInRange(spec);
if (recPos!=0L)
{
	// only have to parse once
	// it is essential that geneVarParser::thisGene has been set!!
	if (spec.weightExpressions.size()|| spec.recWeightExpressions.size())
	{
		numWeights = spec.weightExpressions.size();
		if (numWeights)
			for (std::list<std::string>::const_iterator it = spec.weightExpressions.begin(); it != spec.weightExpressions.end(); ++it)
			{
				geneVarParser* wP = new geneVarParser;
				const char* s = it->c_str();
				wP->parse(s);
				weightParser.push_back(wP);
			}
		numRecWeights = spec.recWeightExpressions.size();
		if (numRecWeights)
			for (std::list<std::string>::const_iterator it = spec.recWeightExpressions.begin(); it != spec.recWeightExpressions.end(); ++it)
			{
				geneVarParser* wP = new geneVarParser;
				const char* s = it->c_str();
				wP->parse(s);
				recWeightParser.push_back(wP);
			}
	}
	else
		numWeights = 1;
	if (spec.excludeExpressions.size())
	{
		for (std::list<std::string>::const_iterator it=spec.excludeExpressions.begin();it!=spec.excludeExpressions.end();++it)
		{
			geneVarParser *eP=new geneVarParser;
			const char *s = it->c_str();
			eP->parse(s);
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
		for (all = 1; all < (spec.mergeAltAlleles?2:tempRecord.nAlls); ++all)
		{
			if (splitLocusCount >= MAXLOCIINSCOREASSOCFILE)
			{
				dcerror.kill();
				dcerror(1, "Number of variants exceeds MAXLOCIINSCOREASSOCFILE (%d) in masterLocusFile::outputSAInfo()\nNeed to increase MAXLOCIINSCOREASSOCFILE and recompile\n", MAXLOCIINSCOREASSOCFILE);
				return 0;
			}
			for (w = 0; w < numWeights+numRecWeights; ++w)
				locusWeights[w][splitLocusCount] = 1.0; // default if nothing else changes it
			useLocus[splitLocusCount] = 1;
			load(tempRecord, recPos);
			doNotUseUntypedFreqfile = 0;
			if (spec.unknownIfUntyped)
			{
				for (i = 0; i < nLocusFiles; ++i)
					if (holdsFreqs[i] && !strcmp(tempRecord.myLocalLocus[i]->filter, "UNTYPED"))
						doNotUseUntypedFreqfile = 1;
			}
			if (doNotUseUntypedFreqfile || (tempRecord.isSNP() != SNP_YES && spec.onlyUseSNPs == 1))
			{
				for (w = 0; w < numWeights+numRecWeights; ++w)
					locusWeights[w][splitLocusCount] = 0.0;
				useLocus[splitLocusCount] = 0;
			}
			else
			{
				// we are going to use one of the inbuilt annotations to see if we exceed consequenceThreshold
				// then if weight funcion specified use that weight instead
				if (spec.consequenceThreshold || (spec.useConsequenceWeights&& (spec.weightExpressions.size() + spec.recWeightExpressions.size())== 0))
				{
					if (spec.useEnsembl)
					{
						if (tempRecord.ensemblConsequence[all][0] == '\0')
							cons = NULL_CONSEQUENCE;
						else // should change below to do this with tables in the same way as the parser does
						{
							cons = getConsequenceType(tempRecord.ensemblConsequence[all], 1);
							if (cons == -1)
							{
								dcerror.kill();
								dcerror(1, "Could not recognise this consequence: %s", tempRecord.ensemblConsequence);
								return 0;
							}
						}
						locusWeights[0][splitLocusCount] = e_consequence[cons].weight;
						if (cons < spec.consequenceThreshold)
							useLocus[splitLocusCount] = 0;
					}
					else
					{
						if (tempRecord.quickConsequence[all][0] == '\0')
							cons = NULL_CONSEQUENCE;
						else
						{
							cons = getConsequenceType(tempRecord.quickConsequence[all], 0);
							if (cons == -1)
							{
								dcerror.kill();
								dcerror(1, "Could not recognise this consequence: %s", tempRecord.quickConsequence);
								return 0;
							}
						}
						locusWeights[0][splitLocusCount] = consequence[cons].weight;
						if (cons < spec.consequenceThreshold)
							useLocus[splitLocusCount] = 0;
					}
				}
				if (spec.weightExpressions.size()+ spec.recWeightExpressions.size() && spec.useConsequenceWeights)
				{
					geneVarParser::thisLocus = &tempRecord;
					geneVarParser::thisAltAllele = all;
					w = 0;
					for (std::list<geneVarParser*>::iterator it = weightParser.begin(); it != weightParser.end(); ++it)
					{
						dcexpr_val* rv = (*it)->eval();
						locusWeights[w++][splitLocusCount] = double(*rv);
						delete rv;
					}
					for (std::list<geneVarParser*>::iterator it = recWeightParser.begin(); it != recWeightParser.end(); ++it)
					{
						dcexpr_val* rv = (*it)->eval();
						locusWeights[w++][splitLocusCount] = double(*rv);
						delete rv;
					}
				}
				if (excludeParser.size() && useLocus[splitLocusCount])
				{
					double val;
					dcexpr_val *rv;
					geneVarParser::thisLocus = &tempRecord;
					geneVarParser::thisWeight = locusWeights[0][splitLocusCount];
					geneVarParser::thisAltAllele = all;
					for (std::list<geneVarParser *>::iterator it = excludeParser.begin(); it != excludeParser.end(); ++it)
					{
						rv = (*it)->eval();
						val = (double)(*rv);
						delete rv;
						if (val != 0)
						{
							useLocus[splitLocusCount] = 0;
							break;
						}
					}
				}
			}
			++splitLocusCount;
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
return splitLocusCount;
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
	recPos = index.get_next();
	testKey=index.current_key();
	if ((c=atoi(testKey))==0 || c>spec.ec || (c==spec.ec && atol(testKey+3)>spec.ep))
		return 0;
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
	int c, i, altIsCommon,ll;
	locusCount = 0;
// hereOK();
	if (gotoFirstInRange(spec))
	do
	{
		outputCurrentProbs(prob[locusCount++], spec); 
	} while (gotoNextInRange(spec));
// hereOK();
	return locusCount;
}

int masterLocusFile::countNumberInRange(analysisSpecs &spec)
{
	int nValid=0;
	if (gotoFirstInRange(spec))
		do
		{
			++nValid;
		} while (gotoNextInRange(spec));
	return nValid;
}

int masterLocusFile::outputAlleles(allelePair **all, analysisSpecs &spec)
{
	int locusCount, subCount;
	FILEPOSITION recPos;
	const char *testKey;
	int c, i, altIsCommon;
 	locusCount = 0;
// hereOK();
	checkSystem();
	if (gotoFirstInRange(spec))
	do
	{
// hereOK();
		nAlls[locusCount] = tempRecord.nAlls;
		outputCurrentAlleles(all[locusCount++], spec);
		if (locusCount>MAXLOCIINSCOREASSOCFILE)
		{
			dcerror.kill();
			dcerror(1,"Number of variants exceeds MAXLOCIINSCOREASSOCFILE (%d) in masterLocusFile::outputAlleles()\nNeed to increase MAXLOCIINSCOREASSOCFILE and recompile\n",MAXLOCIINSCOREASSOCFILE);
			return 0;
		}
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
	fprintf(fp,"nAlls = %d: ",nAlls);
	for (a=0;a<nAlls;++a)
		fprintf(fp,"%s ",alls[a]);
	fprintf(fp,"\nAllele mapping:\n");

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
	fprintf(fp,"\t%s\n",quickConsequence);
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
		dcerror(99,"Cannot add %s to master locus file because it already contains %d locusFiles, which is equal to nLocusFiles",fn,nLocusFiles);
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
		// index.add("NO RECORD",0L); // just something to get index started so it is not empty
		// no longer need this and having this breaks near_find()
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
	char line[MAXALLLENGTH*MAXALL],rest[MAXALLLENGTH*MAXALL],all[MAXALLLENGTH*2];
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
		if(l==rec.nAlls)
		{
			if(rec.nAlls>=MAXALL)
			{
				dcerror.kill();
				dcerror(1,"Total number of alleles for locus at %d:%ld exceeds MAXALL of %d so cannot merge(). Need to increase MAXALL.",rec.chr,rec.pos,MAXALL);
				return 0;
			}
			if(strlen(all)>MAXALLLENGTH)
			{
				dcerror(1,"Length of allele for locus at %d:%ld exceeds MAXALLLENGTH of %d so will be truncated.\n%s",rec.chr,rec.pos,MAXALLLENGTH,all);
				strncpy(rec.alls[rec.nAlls],all,MAXALLLENGTH);
				rec.alls[rec.nAlls][MAXALLLENGTH]='\0';
				rec.nAlls++;

			}
			else
				strcpy(rec.alls[rec.nAlls++],all);
		}
		rec.alleleMapping[currentLocusFile][a]=l;
		strcpy(line,rest);
		rest[0]='\0';
	}
	if (rec.nAlls>2)
		rec.SNP=SNP_NO;
	rec.locusPosInFile[currentLocusFile]=locusPosInFile;
	rec.myLocalLocus[currentLocusFile]->typeSpecificCopy(loc);
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
	id[0]=ref[0]=alt[0]='\0';
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
nAlls = new int[MAXLOCIINSCOREASSOCFILE];
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
	delete[] nAlls;
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
	int i,all;
	nLocusFiles=nLF;
alleleMapping=new alleleMap[nLF];
locusPosInFile=new FILEPOSITION[nLF];
myLocalLocus=new LOCALLOCUSPTR[nLF];
for (i=0;i<nLF;++i)
	myLocalLocus[i]=(LOCALLOCUSPTR) new vcfLocalLocus; 
// one day will have to allocate localLocus of correct type to these
for (all = 0; all < MAXALL; ++all)
{
	ensemblConsequence[all][0] = quickConsequence[all][0] = '\0';
	worstConsequenceType[all] = NULL_CONSEQUENCE;
}
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
if (recPos == 0L)
	recPos = index.get_next(); // because near_find() can just jump to beginning
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
			if (spec.dontMergeAlleles && strcmp(tempRecord.alt, tempLocus->alt)) // treat allele variants as two separate loci
			{
				recPos = 0L;
				break;
			}
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
	if (!merge(tempRecord,tempLocus,locusPosInFile))
		return 0;
}
else
{
	fill(tempRecord,tempLocus,locusPosInFile);
	FSEEK(recordFile,0L,SEEK_END);
	recPos=FTELL(recordFile);
	sprintf(key,"%3d %10ld %10.10s",tempLocus->chr,tempLocus->pos,tempLocus->alt);
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


