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

#ifndef vcfLocusFilHPP
#define vcfLocusFilHPP

#include "masterLocusFile.hpp"

class vcfLocalLocus : localLocus {
	int GTpos,GQpos,GPpos,ADpos; //position of these variables in genotype entry
	float qual;
	char info[VCFFIELDLENGTH];
	char format[VCFFIELDLENGTH];
	int nFieldsToSkip;
public:
	virtual void clear();
	virtual int outputAlleles(allelePair *all,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec);
	virtual int outputProbs(probTriple *all,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec);
	virtual const char *getInfo() { return info; }
	int outputCalls(strEntry *call,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec);
	int outputVcfGenotypes(FILE *fo,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap);
	void parseFormat();
	virtual int read(FILE *fp);
	virtual int write(FILE *fp);
	virtual int input(FILE *f, FILEPOSITION *locusPosInFile, analysisSpecs const &spec);
	vcfLocalLocus() { nFieldsToSkip = DEFAULTNUMVCFFIELDSTOSKIP; clear(); } // need to call clear() here or only get the localLocus version
	~vcfLocalLocus() { ; }
	virtual locusFileType myType() { return VCFFILE; }
	virtual int typeSpecificCopy(localLocus *src);
};

class vcfLocusFile : locusFile {
	virtual locusFileType fileType() { return VCFFILE; }
	int nFieldsToSkip;
public:
	virtual int readHeaderInfo();
	virtual int outputSubNames(strEntry *subName, analysisSpecs &spec);
	vcfLocusFile() { nFieldsToSkip=DEFAULTNUMVCFFIELDSTOSKIP; }
};

// note that these will not work correctly if vcf file has no header because index cannot store a file position of 0L

#endif