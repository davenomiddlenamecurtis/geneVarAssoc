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
	static double thisWeight;
	std::map<std::string,std::string> queryCache;
	geneVarParser();
	~geneVarParser() { ; }
	int init(masterLocus &m,refseqGeneInfo &r) { thisLocus=&m; thisGene=&r; }
	dcexpr_val *eval();
};

class weightTable {
public:
	std::string tableName;
	std::map<std::string,float> weightMap;
	weightTable() { ; }
	~weightTable() { ; }
	int readFromFile(char *n);
	void init(char *n,consequenceReport consequence[],int nConsequence);
};

extern std::map<std::string,weightTable*> weightTableList;

#endif
