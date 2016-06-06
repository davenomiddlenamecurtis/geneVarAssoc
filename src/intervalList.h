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
