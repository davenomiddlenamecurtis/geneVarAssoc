#include "geneVarParser.hpp"

// this will mainly be to parse annotations, assign weights

std::map<std::string,weightTable*> weightTableList;

void weightTable::init(char *n,consequenceReport consequence[],int nConsequence)
{
	int c;
	tableName=n;
	for (c=0;c<nConsequence;++c)
		weightMap[consequence[c].str]=consequence[c].weight;
	weightTableList[tableName]=this;
}

dcexpr_val *performTabixQuery(char *buff,const char *fn,int addChr)
{
	char fnBuff[1000],*ptr,lineBuff[5000],chrStr[10];
	dcexpr_val *rv;
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
	unlink("tabixQueryOutput.txt");
	sprintf(lineBuff,"tabix %s %s%s:%ld-%ld",fnBuff,addChr?"CHR":"",chrStr,geneVarParser::thisLocus->getPos(),geneVarParser::thisLocus->getPos());
	return rv;
}

dcexpr_val *vcfAddChrLookup_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_val *rv;
	char fnBuff[1000],lineBuff[5000];
	strcpy(fnBuff,(char*)(*r2));
	rv=performTabixQuery(lineBuff,fnBuff,1);
	delete r1; delete r2;
	return rv;
}

dcexpr_val *vcfLookup_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_val *rv;
	char fnBuff[1000],lineBuff[5000];
	strcpy(fnBuff,(char*)(*r2));
	rv=performTabixQuery(lineBuff,fnBuff,1);
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
	if(tableIter==weightTableList.end())
	{
		dcerror(1,"Could not find a weight table named %s\n",(char*)(*r2));
	}
	else
	{
		tab=tableIter->second;
		std::map<std::string,float>::const_iterator weightIter=tab->weightMap.find((char*)(*r1));
		if(weightIter==tab->weightMap.end())
		{
			dcerror(1,"Could not find effect %s in weight table named %s\n",(char*)(*r1),(char*)(*r2));
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

dcexpr_val *annot_func(dcvnode *b1)
{
	dcexpr_val *r1;
	EVAL_R1;
	char *ptr;
	char *annot_type=(char*)(*r1);
	dcexpr_string *rv;
	if (!strcmp(annot_type,"INBUILT"))
	{
		geneVarParser::thisLocus->getQuickFeature(*geneVarParser::thisGene);
		rv=new dcexpr_string(geneVarParser::thisLocus->reportQuickConsequence());
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
	char *attrib_type=(char*)(*r1);
	dcexpr_val *rv;
	if (!strcmp(attrib_type,"WEIGHT"))
		rv=new dcexpr_double(geneVarParser::thisWeight);
	else
	{
		dcerror(1,"The ATTRIB function cannot accept %s as its argument",attrib_type);
		rv=new dcexpr_string("UNIMPLEMENTEDATTRIBUTE");
	}
	delete r1;
	return rv;
}

int initGeneVarParser()
{
	add_bin_op_next("STRCAT",strcat_func);
	add_bin_op_next("GETWEIGHT",getWeight_func);
	add_un_op("ANNOT",annot_func);
	add_un_op("ATTRIB",attrib_func);
	return 1;
}

bool geneVarParser::parserIsInited=0;
masterLocus *geneVarParser::thisLocus;
refseqGeneInfo *geneVarParser::thisGene;
double geneVarParser::thisWeight;

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
	if(express::debugFile)
	{
		fprintf(express::debugFile,"Evaluating expression using gene %s and this locus:\n",geneVarParser::thisGene->getGene());
		geneVarParser::thisLocus->print(express::debugFile);
	}
	rv=express::eval();
	if(express::debugFile)
		fprintf(express::debugFile,"Final result: %s\n",(char*)(*rv));
	return rv;
}
