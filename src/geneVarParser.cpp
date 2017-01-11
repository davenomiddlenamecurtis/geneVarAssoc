#include "geneVarParser.hpp"

#include <stdlib.h>q

// this will mainly be to parse annotations, assign weights

std::map<std::string,weightTable*> weightTableList;

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

dcexpr_val *performTabixQuery(const char *fn,int addChr,char *lookupStr)
{
	char fnBuff[1000],*ptr,*tptr,queryBuff[1000],lineBuff[5000],tempBuff[1000],chrStr[10];
	dcexpr_val *rv;
	FILE *fq;
	strcpy(fnBuff,fn);
	int chr=geneVarParser::thisLocus->getChr();
	if (chr==23)
		sprintf(chrStr,"X");
	else
		sprintf(chrStr,"%d",chr);
	if(ptr=strchr(fnBuff,'*'))
	{
		strcpy(ptr,chrStr);
		strcpy(lineBuff,fn);
		ptr=strchr(lineBuff,'*');
		strcat(fnBuff,lineBuff);
	}
	sprintf(queryBuff,"tabix %s %s%s:%ld-%ld",fnBuff,addChr?"CHR":"",chrStr,geneVarParser::thisLocus->getPos(),geneVarParser::thisLocus->getPos());
	std::map<std::string,std::string>::const_iterator queryIter=geneVarParser::queryCache.find(queryBuff);
	if (queryIter==geneVarParser::queryCache.end())
	{
		int stest;
		unlink("tabixQueryOutput.txt");
		sprintf(lineBuff,"%s > tabixQueryOutput.txt",queryBuff);
		checkSystem();
		if ((stest=system(lineBuff))!=0)
		{
			dcerror(1,"Could not execute %s, failed with error %d\n",lineBuff,stest);
		}
		fq=fopen("tabixQueryOutput.txt","r");
		if (fq==0 || fscanf(fq,"%*s %*s %*s %*s %*s %*s %*s %s",lineBuff)!=1)
			sprintf(lineBuff,"NOVCFLINE_%s_%ld",chrStr,geneVarParser::thisLocus->getPos());
		if(fq)
			fclose(fq);
		geneVarParser::queryCache[queryBuff]=lineBuff;
	}
	else
	{
		strcpy(lineBuff,queryIter->second.c_str());
	}
	if (strncmp(lineBuff,"NOVCFLINE",strlen("NOVCFLINE")))
	{
		sprintf(tempBuff,";%s",lookupStr);
		tptr=lineBuff;
		while ((ptr=strstr(tptr,tempBuff))!=0) // possible multiple occurrences, e.g. of AF
		{
			if (ptr[strlen(tempBuff)]=='=')
				break;
			else
				tptr=ptr+strlen(tempBuff);
		}
		if (ptr==0)
			sprintf(lineBuff,"NOVCFENTRY_%s_%ld_%s",chrStr,geneVarParser::thisLocus->getPos(),lookupStr);
		else
		{
			ptr+=strlen(tempBuff)+1;
			tptr=tempBuff;
			while(*ptr && *ptr!=';')
				*tptr++=*ptr++;
			*tptr=0;
			strcpy(lineBuff,tempBuff);
		}
	}
	rv=new dcexpr_string(lineBuff);
	return rv;
}

dcexpr_val *vcfAddChrLookup_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_val *rv;
	rv=performTabixQuery((char*)(*r2),1,(char*)(*r1));
	delete r1; delete r2;
	return rv;
}

dcexpr_val *vcfLookup_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_val *rv;
	rv=performTabixQuery((char*)(*r2),0,(char*)(*r1));
	delete r1; delete r2;
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
		rv=new dcexpr_string(geneVarParser::thisLocus->reportQuickConsequence());
		ptr=(char*)(*rv);
		while(*ptr && !isspace(*ptr))
			++ptr;
		*ptr='\0'; // remove e.g. gene name and so on
	}
	else if (!strcmp(annot_type,"VEP"))
	{
		rv=new dcexpr_string(geneVarParser::thisLocus->reportEnsemblConsequence());
		ptr=(char*)(*rv);
		while(*ptr && !isspace(*ptr))
			++ptr;
		*ptr='\0'; // remove e.g. gene name and so on
	}
	else
	{
		dcerror(1,"The ANNOT function cannot accept %s as its argument",annot_type);
		rv=new dcexpr_string("UNIMPLEMENTEDANNOTATION");
	}
	delete r1;
	return rv;
}

dcexpr_val *attrib_func(dcvnode *b1)
{
	dcexpr_val *r1;
	EVAL_R1;
	char buff[100];
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
		sprintf(buff,"%d",geneVarParser::thisLocus->getChr());
		rv=new dcexpr_string(buff);
	}
	else if (!strcmp(attrib_type,"COORD"))
	{
		sprintf(buff,"%d:%ld",geneVarParser::thisLocus->getChr(),geneVarParser::thisLocus->getPos());
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
		dcerror(1,"The ATTRIB function cannot accept %s as its argument",attrib_type);
		rv=new dcexpr_string("UNIMPLEMENTEDATTRIBUTE");
	}
	delete r1;
	return rv;
}

bool geneVarParser::parserIsInited=0;
masterLocus *geneVarParser::thisLocus;
refseqGeneInfo *geneVarParser::thisGene;
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
		fprintf(express::debugFile,"Evaluating expression using gene %s and this locus:\n",geneVarParser::thisGene->getGene());
		geneVarParser::thisLocus->print(express::debugFile);
	}
	rv=express::eval();
	if (express::debugFile)
		fprintf(express::debugFile,"Final result: %s\n",(char*)(*rv));
	return rv;
}

int initGeneVarParser()
{
	add_bin_op_next("STARTSWITH",startsWith_func);
	add_bin_op_same("STRCAT",strcat_func);
	add_bin_op_same("GETWEIGHT",getWeight_func);
	add_bin_op_same("VCFLOOKUP",vcfLookup_func);
	add_bin_op_same("VCFADDCHRLOOKUP",vcfAddChrLookup_func);
	add_un_op("ANNOT",annot_func);
	add_un_op("ATTRIB",attrib_func);
	return 1;
}

