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

#ifndef hapsLocusFileHPP
#define hapsLocusFileHPP

#include "masterLocusFile.hpp"

class hapsLocalLocus : localLocus {
	char rsName[100];
public:
	virtual int outputProbs(probTriple *all, FILE *f, FILEPOSITION filePos, int nSubs, int *alleleMap, analysisSpecs const &spec) { return 0; } // not implemented, though could be
	virtual int outputAlleles(allelePair *all,FILE *f,FILEPOSITION filePos,int nSubs,int *alleleMap,analysisSpecs const &spec);
	virtual int input(FILE *f, FILEPOSITION *locusPosInFile, analysisSpecs const &spec);
	virtual void clear() { rsName[0] = '\0'; }
	virtual int read(FILE *fp);
	virtual int write(FILE *fp);
	virtual int typeSpecificCopy(localLocus *src);
	hapsLocalLocus() { clear(); }
	~hapsLocalLocus() { ; }
	virtual locusFileType myType() { return SHAPEITHAPSFILE; }
};

class hapsLocusFile : locusFile {
	virtual locusFileType fileType() { return SHAPEITHAPSFILE; }
public:
	virtual int readHeaderInfo();
	virtual int outputSubNames(strEntry *subName, analysisSpecs &spec);
	hapsLocusFile() { ; }
};

#endif