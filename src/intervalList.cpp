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

#include "intervalList.h"
#include <stdlib.h>
#include <string.h>

#ifndef MSDOS
#define __cdecl
#endif


int __cdecl intervalCompare(const void *ii1,const void *ii2)
{
	interval *i1,*i2;
	int c1, c2;
	i1=(interval*)ii1;
	i2=(interval*)ii2;
	if ((i1->chr[0] == 'X') != (i2->chr[0] == 'X'))
		return((i1->chr[0] == 'X') ? 1 : -1);
	else if ((c1 = atoi(i1->chr)) != (c2 = atoi(i2->chr)))
		return(c1>c2?1:c1<c2?-1:0);
	else
		return(i1->st>i2->st?1:i1->st<i2->st?-1:0);
}


intervalList::intervalList(int mI)
{
	ints=0;
	nInts=0;
	resize(mI);
}

void intervalList::resize(int mI)
{
	interval *newInts;
	newInts=new interval[mI];
	if (ints!=0)
	{
		for (int i=0;i<nInts;++i)
			newInts[i]=ints[i];
		delete [] ints;
	}
	ints=newInts;
	maxInts=mI;
}

void intervalList::sort()
{
	qsort(ints,nInts,sizeof(interval),intervalCompare);
}

void intervalList::merge()
{
	// any overlapping intervals get merged
	// must have been sorted first
	// ignores 0,1 indexing issues, assumes e.g. all 1-based
	int i,j;
	for (i=0;i<nInts-1;)
	{
		if (!strcmp(ints[i].chr, ints[i+1].chr) && ints[i].en>=ints[i+1].st)
		{
			if (ints[i+1].en>ints[i].en)
				ints[i].en=ints[i+1].en;
			for (j=i+1;j<nInts-1;++j)
				ints[j]=ints[j+1];
			--nInts;
		}
		else
			++i;
	}
}

void intervalList::append(int s, int e)
{
	if (nInts > maxInts - 1)
		resize(maxInts + 10);
	ints[nInts].st = s;
	ints[nInts].en = e;
	++nInts;
}

void intervalList::append(char *c,int s, int e)
{
	if (nInts > maxInts - 1)
		resize(maxInts + 10);
	strcpy(ints[nInts].chr, c);
	ints[nInts].st = s;
	ints[nInts].en = e;
	++nInts;
}

void intervalList::append(intervalList &iL)
{
	if (maxInts<nInts+iL.nInts)
		resize(nInts+iL.nInts);
	for (int i=0;i<iL.nInts;++i)
		ints[nInts+i]=iL.ints[i];
	nInts=nInts+iL.nInts;
}

int intervalList::getAllMatches(char *seq,char *toMatch,int start)
// used by refseqTranscript::get3PrimeMatches() which is not used
// looks like it is to search a 3prime sequence for particular sequence matches
// NB it does not set chr
{
	char *ptr;
	int foundAny,s,e;
	foundAny=0;
	ptr=seq;
	while((ptr=strstr(ptr,toMatch))!=0)
	{
		foundAny=1;
		s=ptr-seq+start;
		e=s+strlen(toMatch)-1;
		append(s,e);
		ptr=ptr+1;
	}
	return foundAny;
}