#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <assert.h>
#include <map>
using namespace std;

#define MAXIDSTRLEN 5000

class ccCounts {
public:
int count[2];
string IDstr;
ccCounts() { clear(); }
void clear() { count[0]=count[1]=0; IDstr=""; }
void increment(int i) { count[i]++; }
};

std::map <string, ccCounts> CNVcounts;
std::map <string, ccCounts>::iterator oneCount;

int main(int argc,char *argv[])
{
	FILE *fi,*fo;
	ccCounts temp;
	char gene[20],line[20000],geneList[10000],rest[10000],ID[100],vep_annot[100],call_qc[100];
	int cc,cnv_size;
	float common;
	for (cc=0;cc<2;++cc)
	{
		fi=fopen(argv[1+cc],"r");
		assert(fi!=0);
		fgets(line,1999,fi);
		while(fgets(line,1999,fi))
		{
			sscanf(line,"%*d %*d %*d %*d %*f %*s %s %*f %*f %*f %*f %f %*f %*f %*f %d %s %s %s ",ID,&common,&cnv_size,call_qc,vep_annot,geneList); 
			if (common!=0 || strcmp(call_qc,"PASS"))
				continue;
			while (*rest=0,sscanf(geneList,"%[^;];%s",gene,rest)>=1)
			{
				if (!strcmp(gene,"HIGHSENS"))
					printf("error");
				if (!strcmp(gene,"NA"))
					break;
				oneCount=CNVcounts.find(gene);
				if (oneCount == CNVcounts.end())
					temp.clear();
				else
					temp=oneCount->second;
				temp.increment(cc);
				temp.IDstr=temp.IDstr+ID+"\t";
				CNVcounts[gene]=temp;
				strcpy(geneList,rest);
			}
		}
		fclose(fi);
	}
	fo=fopen(argv[3],"w");
	for (oneCount=CNVcounts.begin(); oneCount!=CNVcounts.end();oneCount++)
	{
		string str;
		str=oneCount->first;
		strcpy(gene,str.c_str());
		if (oneCount->second.count[0]+oneCount->second.count[1]<200)
			fprintf(fo,"%s\t%5d\t%5d\t%s\n",gene,oneCount->second.count[0],oneCount->second.count[1],oneCount->second.IDstr.c_str());
	}
	fclose(fo);
	return 0;
}

