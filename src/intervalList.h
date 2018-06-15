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

class interval {
public:
	int st,en;
};

class intervalList {
	int maxInts;
	void resize(int mI);
public:
	int nInts;
	interval *ints;
	intervalList(int mI=10);
	~intervalList() { if (ints!=0) delete[] ints; }
	void append(intervalList &iL);
	void append(int s,int e);
	void sort();
	void merge();
	int getAllMatches(char *seq,char *toMatch,int start);
};
