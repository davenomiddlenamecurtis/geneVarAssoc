#using <System.dll>

using namespace System;
using namespace System::Collections;
using namespace System::Collections::Specialized;

#include <stdio.h>
#include <string.h>

#if 0
sometimes same variant ends up in more than one gene
#endif

int main(int argc,char *argv[])
{
	StringCollection^ myCol = gcnew StringCollection;
	String^ testString;
	char line[1001],line2[1001],testChars[1000];
	int chr,pos,gc[2][3],single;
	FILE *fi,*fo,*fi2,*fo2;
	if (argc>3)
		single=1;
	else
		single=0;
	fi=fopen(argv[1],"r");
	fo=fopen(argv[2],"w");
	if (single)
	{
		fi2=fopen(argv[3],"r");
		fo2=fopen(argv[4],"w");
	}
	while (fgets(line,1000,fi))
	{
		if (!single)
		{
			sscanf(line, "%*s %d %d %*[^\t] || %d %d %d || %d %d %d ",
				&chr, &pos,
				&gc[0][0], &gc[0][1], &gc[0][2],
				&gc[1][0], &gc[1][1], &gc[1][2]);
			sprintf(testChars, "%d %d %d %d %d %d %d %d",
				chr, pos,
				gc[0][0], gc[0][1], gc[0][2],
				gc[1][0], gc[1][1], gc[1][2]);
		}
		else
		{
			sscanf(line,"%*ld %*s %d %d %d %d %d",
				&chr, &pos,
				&gc[0][0], &gc[0][1], &gc[0][2]);
			fgets(line2,1000,fi2);
			sscanf(line2,"%*ld %*s %d %d %d %d %d",
				&chr, &pos,
				&gc[1][0], &gc[1][1], &gc[1][2]);
			sprintf(testChars, "%d %d %d %d %d %d %d %d",
				chr, pos,
				gc[0][0], gc[0][1], gc[0][2],
				gc[1][0], gc[1][1], gc[1][2]);
		}
		testString=gcnew String(testChars);
		printf(line);
		if (myCol->Contains(testString))
			;
		else
		{
			myCol->Add(testString);
			fprintf(fo,"%s",line);
			if (single)
				fprintf(fo2,"%s",line2);
		}

	}
	fclose(fi);
	fclose(fo);
	if (single)
	{
		fclose(fi2);
		fclose(fo2);
	}
	return 0;
}
