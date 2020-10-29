class interval {
public:
	int st,en;
	char chr[3];
	interval() { st = en = 0; chr[0] = '\0'; }
// in many applications chr will simply be ignored
// it is 1 - 22 or "X"
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
	void append(int s, int e);
	void append(char *c,int s, int e);
	void sort();
	void merge();
	int getAllMatches(char *seq,char *toMatch,int start);
};
