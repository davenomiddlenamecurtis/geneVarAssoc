#include "intervalList.h"
#include <stdlib.h>
#include <string.h>

#ifndef MSDOS
#define __cdecl
#endif


int __cdecl intervalCompare(const void *ii1,const void *ii2)
{
	interval *i1,*i2;
	i1=(interval*)ii1;
	i2=(interval*)ii2;
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
		if (ints[i].en>=ints[i+1].st)
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

void intervalList::append(int s,int e)
{
	if (nInts>maxInts-1)
		resize(maxInts+10);
	ints[nInts].st=s;
	ints[nInts].en=e;
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