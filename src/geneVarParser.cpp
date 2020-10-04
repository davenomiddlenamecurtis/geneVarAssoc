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

#include "geneVarParser.hpp"

#include <stdlib.h>

// this will mainly be to parse annotations, assign weights

std::map<std::string,weightTable*> weightTableList;
enum { NSIFTCONSEQUENCES=3 };
consequenceReport sift_consequence[]={
	{ 0, "sift_missing", 1.0 },
	{ 1, "tolerated", 10.0 },
	{ 2, "deleterious", 20.0 }
};

enum { NPOLYPHENCONSEQUENCES=5 };
consequenceReport polyphen_consequence[]={
	{ 0, "polyphen_missing", 1.0 },
	{ 1, "benign", 5.0 },
	{ 2, "unknown", 10.0 },
	{ 3, "possibly_damaging", 20.0 },
	{ 4, "probably_damaging", 20.0 }
};

int weightTable::readFromFile(char *fn,char *n)
{
	char line[1001],buff[1001];
	float w;
	int l;
	FILE *fp;
	fp=fopen(fn,"r");
	if (fp==0)
	{
		dcerror(1,"Could not open weight table file: %s\n",fn);
		return 0;
	}
	tableName=n;
	for (l=0;fgets(line,1000,fp)&&sscanf(line,"%s %f",buff,&w)==2;++l)
	{
		weightMap[buff]=w;
	}
	if (l==0)
		dcerror(1,"No valid lines in weight table file: %s\n",fn);
	return l;
}


void weightTable::init(char *n,consequenceReport consequence[],int nConsequence)
{
	int c;
	tableName=n;
	for (c=0;c<nConsequence;++c)
		weightMap[consequence[c].str]=consequence[c].weight;
	weightTableList[tableName]=this;
}

#define MAXINFOLENGTH 20000
#define MAXINFOLENGTHSTR "20000"
char lineBuff[MAXINFOLENGTH+1],tempBuff[MAXINFOLENGTH+1]; // need these to be big for e.g. VEP on TTN

dcexpr_val *performTabixQuery(const char *fn,int addChr,int lower,char *lookupStr,int convert23toX)
{
	char fnBuff[1000],*ptr,*tptr,queryBuff[1000],chrStr[10],altAll[1000];
	int noEntry,c,f,l;
	dcexpr_val *rv;
	FILE *fq;
	strcpy(fnBuff,fn);
	int chr=geneVarParser::thisLocus->getChr();
	if (chr==23 && convert23toX)
		sprintf(chrStr,"X"); // this is going to be optional
	else
		sprintf(chrStr,"%d",chr);
	if(ptr=strchr(fnBuff,'*'))
	{
		strcpy(ptr,chrStr);
		strcpy(lineBuff,fn);
		ptr=strchr(lineBuff,'*');
		strcat(fnBuff,lineBuff);
	}
	sprintf(queryBuff,"%s:%ld-%s/%s",chrStr,geneVarParser::thisLocus->getPos(), geneVarParser::thisLocus->getAll(0), geneVarParser::thisLocus->getAll(1));
	std::map<std::string,std::string>::const_iterator queryIter=geneVarParser::queryCache.find(queryBuff);
	if (queryIter==geneVarParser::queryCache.end())
	{
		int stest;
		remove("tabixQueryOutput.txt");
		sprintf(lineBuff, "tabix %s %s%s:%ld-%ld > tabixQueryOutput.txt", fnBuff, addChr ? lower ? "chr" : "CHR" : "", chrStr, geneVarParser::thisLocus->getPos(), geneVarParser::thisLocus->getPos());
		checkSystem();
		if ((stest=system(lineBuff))!=0)
		{
			dcerror(1,"Could not execute %s, failed with error %d\n",lineBuff,stest);
		}
		fq=fopen("tabixQueryOutput.txt","r");
		noEntry=1;
		if (fq)
		{
			while (fgets(tempBuff, MAXINFOLENGTH, fq))
			{
				if (sscanf(tempBuff, "%*s %*s %*s %*s %[^ \t,]", altAll) == 1 && !strcmp(altAll, geneVarParser::thisLocus->getAll(1)))
				{
					noEntry = 0;
					sscanf(tempBuff, "%*s %*s %*s %*s %*s %*s %*s %" MAXINFOLENGTHSTR "s", lineBuff);
					break;
				}
			}
			fclose(fq);
		}
		if (noEntry)
			sprintf(lineBuff,"NOVCFLINE_%s_%ld_%s",chrStr, geneVarParser::thisLocus->getPos(), geneVarParser::thisLocus->getAll(1));
		geneVarParser::queryCache[queryBuff]=lineBuff;
	}
	else
	{
		strcpy(lineBuff,queryIter->second.c_str());
	}
	if (strncmp(lineBuff,"NOVCFLINE",strlen("NOVCFLINE")))
	{
		if (!strncmp(lookupStr,lineBuff,strlen(lookupStr)) && lineBuff[strlen(lookupStr)]=='=') // first entry in info field
			ptr=lineBuff+strlen(lookupStr)+1;
		else
		{
			sprintf(tempBuff,";%s=",lookupStr);
			if ((ptr=strstr(lineBuff,tempBuff))!=0)
				ptr+=strlen(tempBuff);
		}
		if (ptr==0)
			sprintf(lineBuff,"NOVCFENTRY_%s_%ld_%s",chrStr,geneVarParser::thisLocus->getPos(),lookupStr);
		else
		{
			tptr=tempBuff;
			while(*ptr && *ptr!=';' && !isspace(*ptr))
				*tptr++=*ptr++;
			*tptr=0;
			strcpy(lineBuff,tempBuff);
		}
	}
	rv=new dcexpr_string(lineBuff);
	return rv;
}

dcexpr_val* vcfAddChrLookup_func(dcvnode* b1, dcvnode* b2)
{
	dcexpr_val* r1, * r2;
	EVAL_BOTH;
	dcexpr_val* rv;
	rv = performTabixQuery((char*)(*r2), 1, 0, (char*)(*r1), 0);
	delete r1; delete r2;
	return rv;
}

dcexpr_val* vcfAddChrLookup23toX_func(dcvnode* b1, dcvnode* b2)
{
	dcexpr_val* r1, * r2;
	EVAL_BOTH;
	dcexpr_val* rv;
	rv = performTabixQuery((char*)(*r2), 1, 0, (char*)(*r1), 1);
	delete r1; delete r2;
	return rv;
}

dcexpr_val* vcfAddLowerChrLookup_func(dcvnode* b1, dcvnode* b2)
{
	dcexpr_val* r1, * r2;
	EVAL_BOTH;
	dcexpr_val* rv;
	rv = performTabixQuery((char*)(*r2), 1, 1, (char*)(*r1), 0);
	delete r1; delete r2;
	return rv;
}

dcexpr_val* vcfAddLowerChrLookup23toX_func(dcvnode* b1, dcvnode* b2)
{
	dcexpr_val* r1, * r2;
	EVAL_BOTH;
	dcexpr_val* rv;
	rv = performTabixQuery((char*)(*r2), 1, 1, (char*)(*r1), 1);
	delete r1; delete r2;
	return rv;
}

dcexpr_val* vcfLookup_func(dcvnode* b1, dcvnode* b2)
{
	dcexpr_val* r1, * r2;
	EVAL_BOTH;
	dcexpr_val* rv;
	rv = performTabixQuery((char*)(*r2), 0, 0, (char*)(*r1), 0);
	delete r1; delete r2;
	return rv;
}

dcexpr_val* vcfLookup23toX_func(dcvnode* b1, dcvnode* b2)
{
	dcexpr_val* r1, * r2;
	EVAL_BOTH;
	dcexpr_val* rv;
	rv = performTabixQuery((char*)(*r2), 0, 0, (char*)(*r1), 1);
	delete r1; delete r2;
	return rv;
}

consequenceReport findWorstConsequence(char *s,consequenceReport *r,int n) // this will be used by the function which uses the CSQ entry
// depends on consequences being listed in order from least to most severe
{
	int c;
	char *ptr;
	for (c=n-1;c>0;--c)
	{
		if ((ptr=strstr(s,r[c].str))!=0)
			return r[c];
	}
	return r[0];
}

dcexpr_string *getAlleleAnnotation(dcexpr_val *r1)
{
	const char *altAllStr;
	char *ptr,*sptr;
	char *CSQEntry = (char*)(*r1);
	dcexpr_string *rv;
	altAllStr = geneVarParser::thisLocus->getAll(geneVarParser::thisAltAllele);
	// allele identifier looks like this CSQ=T| or this ,T| and there can be multiple entries for each allele
	lineBuff[0] = '\0';
	sprintf(tempBuff, "%s|", altAllStr);
	if (!strncmp(CSQEntry, altAllStr,strlen(altAllStr)))
	{
		ptr = CSQEntry;
		sptr = lineBuff;
		while (*ptr && !isspace(*ptr) && *ptr != ',')
			*sptr++ = *ptr++;
		*sptr++ = ','; // always append comma, probably won't hurt
		*sptr = '\0';
	}
	sprintf(tempBuff, ",%s|", altAllStr);
	while ((ptr = strstr(CSQEntry, tempBuff))!=0)
	{
		++ptr; // omit comma
		sptr = strchr(lineBuff,'\0');
		while (*ptr && !isspace(*ptr) && *ptr != ',')
			*sptr++ = *ptr++;
		*sptr++ = ','; // always append comma, probably won't hurt
		*sptr = '\0';
	}
	if (lineBuff[0] == '\0')
	{
		dcerror.warn();
		dcerror(1, 
			"Failed to find annotation at %d:%ld for allele %s in this string:\n%s\n\nIf you use --merge-alt-alleles 0 then there must be an allele-specific annotation in the VEP output and you should not use --pick when running VEP. However it is normal for VEP to rename alleles of indels so will use whole string to select most severe variant.\n",
			geneVarParser::thisLocus->getChr(),geneVarParser::thisLocus->getPos(),altAllStr,CSQEntry);
		strcpy(lineBuff, CSQEntry);
	}
	rv = new dcexpr_string(lineBuff);
	return rv;
}

dcexpr_val *extract_sift_func(dcvnode *b1)
{
	dcexpr_val *r1;
	EVAL_R1;
	char *ptr;
	char *annotation;
	dcexpr_string *alleleSpecificAnnotation = 0;
	if (geneVarParser::mergeAltAlleles == 1)
		annotation = (char*)(*r1);
	else
	{
		alleleSpecificAnnotation = getAlleleAnnotation(r1);
		annotation = (char*)(*alleleSpecificAnnotation);
	}
	dcexpr_string *rv;
	rv = new dcexpr_string(findWorstConsequence(annotation, sift_consequence, NSIFTCONSEQUENCES).str);
	delete r1;
	if (alleleSpecificAnnotation)
		delete alleleSpecificAnnotation;
	return rv;
}

dcexpr_val *extract_polyphen_func(dcvnode *b1)
{
	dcexpr_val *r1;
	EVAL_R1;
	char *annotation;
	dcexpr_string *alleleSpecificAnnotation = 0;
	if (geneVarParser::mergeAltAlleles == 1)
		annotation = (char*)(*r1);
	else
	{
		alleleSpecificAnnotation = getAlleleAnnotation(r1);
		annotation = (char*)(*alleleSpecificAnnotation);
	}
	dcexpr_string *rv;
	rv = new dcexpr_string(findWorstConsequence(annotation, polyphen_consequence, NPOLYPHENCONSEQUENCES).str);
	delete r1;
	if (alleleSpecificAnnotation)
		delete alleleSpecificAnnotation;
	return rv;
}

dcexpr_val *extract_vep_func(dcvnode *b1)
{
	dcexpr_val *r1;
	EVAL_R1;
	char *annotation;
	dcexpr_string *alleleSpecificAnnotation=0;
	if (geneVarParser::mergeAltAlleles == 1)
		annotation = (char*)(*r1);
	else
	{
		alleleSpecificAnnotation = getAlleleAnnotation(r1);
		annotation = (char*)(*alleleSpecificAnnotation);
	}
	dcexpr_string *rv;
	rv=new dcexpr_string(findWorstConsequence(annotation,e_consequence,E_NCONSEQUENCETYPES).str);
	delete r1;
	if (alleleSpecificAnnotation)
		delete alleleSpecificAnnotation;
	return rv;
}

dcexpr_val *getWeight_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_double *rv=new dcexpr_double(0);
	weightTable* tab;
	std::map<std::string,weightTable*>::const_iterator tableIter=weightTableList.find((char*)(*r2));
	if (tableIter!=weightTableList.end())
		tab=tableIter->second;
	else
	{
		char *fn=(char*)(*r2);
		tab=new weightTable;
		if (tab->readFromFile(fn,fn)) // will write error if problem
			weightTableList[fn]=tab;
		else
			tab=0;
	}
	if (tab!=0)
	{
		std::map<std::string,float>::const_iterator weightIter=tab->weightMap.find((char*)(*r1));
		if(weightIter==tab->weightMap.end())
		{
			dcerror(1,"Could not find annotation %s in weight table named %s\n",(char*)(*r1),(char*)(*r2));
		}
		else
			*rv=weightIter->second;
	}
	delete r1; delete r2;
	return rv;
}

dcexpr_val *strcat_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_string *rv=new dcexpr_string((char*)(*r1),strlen((char*)(*r1))+strlen((char*)(*r2)));
	strcat((char*)rv,(char*)(*r2));
	delete r1; delete r2;
	return rv;
}

dcexpr_val *startsWith_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_double *rv=new dcexpr_double(strncmp((char*)(*r1),(char*)(*r2),strlen((char*)(*r2)))?0:1);
	delete r1; delete r2;
	return rv;
}

dcexpr_val *annot_func(dcvnode *b1)
{
	dcexpr_val *r1;
	EVAL_R1;
	char *ptr;
	char *annot_type=(char*)(*r1);
	dcexpr_string *rv;
	if (!strcmp(annot_type,"INBUILT"))
	{
		// assume this is already done
		// geneVarParser::thisLocus->getQuickFeature(*geneVarParser::thisGene);
		rv=new dcexpr_string(geneVarParser::thisLocus->reportQuickConsequence(geneVarParser::thisAltAllele));
		ptr=(char*)(*rv);
		while(*ptr && !isspace(*ptr))
			++ptr;
		*ptr='\0'; // remove e.g. gene name and so on
	}
	else if (!strcmp(annot_type,"VEP"))
	{
		rv=new dcexpr_string(geneVarParser::thisLocus->reportEnsemblConsequence(geneVarParser::thisAltAllele));
		ptr=(char*)(*rv);
		while(*ptr && !isspace(*ptr))
			++ptr;
		*ptr='\0'; // remove e.g. gene name and so on
	}
	else
	{
		dcerror(1,"The ANNOT function cannot accept %s as its argument\n",annot_type);
		rv=new dcexpr_string("UNIMPLEMENTEDANNOTATION");
	}
	delete r1;
	return rv;
}

dcexpr_val *attrib_func(dcvnode *b1)
{
	dcexpr_val *r1;
	EVAL_R1;
	char buff[100],chrStr[20];
	int chr;
	char *attrib_type=(char*)(*r1);
	dcexpr_val *rv;
	if (!strcmp(attrib_type,"WEIGHT"))
		rv=new dcexpr_double(geneVarParser::thisWeight);
	else if (!strcmp(attrib_type,"POS"))
	{
		sprintf(buff,"%ld",geneVarParser::thisLocus->getPos());
		rv=new dcexpr_string(buff);
	}
	else if (!strcmp(attrib_type,"CHR"))
	{
		if ((chr=geneVarParser::thisLocus->getChr())==23)
			strcpy(chrStr,"X");
		else
			sprintf(chrStr,"%d",chr);
		rv=new dcexpr_string(chrStr);
	}
	else if (!strcmp(attrib_type,"COORD"))
	{
		if ((chr=geneVarParser::thisLocus->getChr())==23)
			strcpy(chrStr,"X");
		else
			sprintf(chrStr,"%d",chr);
		sprintf(buff,"%s:%ld",chrStr,geneVarParser::thisLocus->getPos());
		rv=new dcexpr_string(buff);
	}
	else if (!strcmp(attrib_type,"ID"))
	{
		rv=new dcexpr_string(geneVarParser::thisLocus->getID());
	}
	else if (!strcmp(attrib_type,"ISSNP"))
	{
		rv=new dcexpr_double((double)geneVarParser::thisLocus->isSNP());
	}
	else
	{
		dcerror(1,"The ATTRIB function cannot accept %s as its argument\n",attrib_type);
		rv=new dcexpr_string("UNIMPLEMENTEDATTRIBUTE");
	}
	delete r1;
	return rv;
}

bool geneVarParser::parserIsInited=0;
masterLocus *geneVarParser::thisLocus;
int geneVarParser::thisAltAllele;
int geneVarParser::mergeAltAlleles;
refseqGeneInfo *geneVarParser::thisGene=0; // because not used by intVarAssoc
double geneVarParser::thisWeight;
std::map<std::string,std::string> geneVarParser::queryCache;
extern int initGeneVarParser();

geneVarParser::geneVarParser()
{
	if(!parserIsInited)
	{
		initGeneVarParser();
		parserIsInited=1;
	}
}

dcexpr_val *geneVarParser::eval()
{
	dcexpr_val *rv;
	if (express::debugFile)
	{
		if (geneVarParser::thisGene)
			fprintf(express::debugFile, "Evaluating expression using gene %s and variant at %d:%ld:\n",
				geneVarParser::thisGene->getGene(),
				geneVarParser::thisLocus->getChr(), geneVarParser::thisLocus->getPos());
		else
			fprintf(express::debugFile, "Evaluating expression using variant at %d:%ld:\n",
				geneVarParser::thisLocus->getChr(), geneVarParser::thisLocus->getPos());
	}

	rv = express::eval();
	if (express::debugFile)
		fprintf(express::debugFile,"Final result of expression evaluation: %s\n",(char*)(*rv));
	return rv;
}

int initGeneVarParser()
{
	weightTable *wt;
	wt=new weightTable;
	wt->init("DEFAULTSIFTWEIGHTS",sift_consequence,NSIFTCONSEQUENCES);
	weightTableList[wt->tableName]=wt;
	wt=new weightTable;
	wt->init("DEFAULTPOLYPHENWEIGHTS",polyphen_consequence,NPOLYPHENCONSEQUENCES);
	weightTableList[wt->tableName]=wt;

	add_bin_op_next("STARTSWITH",startsWith_func);
	add_bin_op_same("STRCAT",strcat_func);
	add_bin_op_same("GETWEIGHT",getWeight_func);
	add_bin_op_same("VCFLOOKUP", vcfLookup_func);
	add_bin_op_same("VCFADDCHRLOOKUP", vcfAddChrLookup_func);
	add_bin_op_same("VCFADDLOWERCHRLOOKUP", vcfAddLowerChrLookup_func);
	add_bin_op_same("VCFLOOKUP23TOX", vcfLookup23toX_func);
	add_bin_op_same("VCFADDCHRLOOKUP23TOX", vcfAddChrLookup23toX_func);
	add_bin_op_same("VCFADDLOWERCHRLOOKUP23TOX", vcfAddLowerChrLookup23toX_func);
	add_un_op("ANNOT",annot_func);
	add_un_op("ATTRIB",attrib_func);
	add_un_op("GETPOLYPHEN",extract_polyphen_func);
	add_un_op("GETSIFT",extract_sift_func);
	add_un_op("GETVEP",extract_vep_func);
	return 1;
}

