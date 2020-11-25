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

#ifndef GENEVARPARSERHPP
#define GENEVARPARSERHPP
#include "dcexpr.hpp"
#include "geneVarUtils.hpp"
#include "getGene.hpp"
#include "consequenceType.hpp"

class geneVarParser : public express {
	static bool parserIsInited;
public:
	static masterLocus *thisLocus;
	static refseqGeneInfo *thisGene;
	static int thisAltAllele;
	static double thisWeight;
	static std::map<std::string,std::string> queryCache;
	static int mergeAltAlleles,multilineVEP;
	geneVarParser();
	~geneVarParser() { ; }
	int init(masterLocus &m, refseqGeneInfo &r) { thisLocus = &m; thisGene = &r; thisAltAllele = 0; }
	dcexpr_val *eval();
};

class weightTable {
public:
	std::string tableName;
	std::map<std::string,float> weightMap;
	weightTable() { ; }
	~weightTable() { ; }
	int readFromFile(char *fn,char *n);
	void init(char *n,consequenceReport consequence[],int nConsequence);
};

extern std::map<std::string,weightTable*> weightTableList;

#endif
