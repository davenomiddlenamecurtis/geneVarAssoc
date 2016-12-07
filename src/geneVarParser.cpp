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

dcexpr_val *getWeight_func(dcvnode* b1,dcvnode *b2)
{
	dcexpr_val *r1,*r2;
	EVAL_BOTH;
	dcexpr_double *rv=new dcexpr_double(0);
	weightTable* tab;
	std::map<std::string,weightTable*>::const_iterator tableIter=weightTableList.find((char*)(*r2));
	if (tableIter==weightTableList.end())
	{
		dcerror(1,"Could not find a weight table named %s\n",(char*)(*r2));
	}
	else
	{
		tab=tableIter->second;
		std::map<std::string,float>::const_iterator weightIter=tab->weightMap.find((char*)(*r1));
		if (weightIter==tab->weightMap.end())
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
	if(strcmp(annot_type,"INBUILT"))
	{
		dcerror(1,"The ANNOT function currently only accepts INBUILT as its argument");
		rv=new dcexpr_string("ONLYINBUILTANNOTATIONISIMPLEMENTED");
	}
	else
	{
		geneVarParser::thisLocus->getQuickFeature(*geneVarParser::thisGene);
		rv=new dcexpr_string(geneVarParser::thisLocus->reportQuickConsequence());
		ptr=(char*)(*rv);
		while (*ptr && !isspace(*ptr))
			++ptr;
		*ptr='\0'; // remove e.g. gene name and so on
	}
	delete r1;
	return rv;
}

int initGeneVarParser()
{
	add_bin_op_next("STRCAT",strcat_func);
	add_bin_op_next("GETWEIGHT",getWeight_func);
	add_un_op("ANNOT",annot_func);
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

