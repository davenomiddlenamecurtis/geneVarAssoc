#include "dcerror.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define FIELDLENGTH 40
#define MAXSUBS 10000

char line[FIELDLENGTH*(MAXSUBS + 20)], right[MAXSUBS + 20][FIELDLENGTH],wrong[MAXSUBS + 20][FIELDLENGTH];
int index[MAXSUBS + 20];

int fillArray(char* l, char arr[MAXSUBS + 20][FIELDLENGTH])
{
	char *sptr, *tptr,ch;
	int i;
	i = 0;
	sptr = l;
	while (1)
	{
		tptr = arr[i];
		while ((ch = *sptr++) && !isspace(ch))
			*tptr++ = ch;
		*tptr = '\0';
		++i;
		if (!ch)
			goto done;
		while ((ch = *sptr) && isspace(ch))
			++sptr;
		if (!ch)
			goto done;
	}
	done:
	return i;
}

int main(int argc, char *argv[])
{
	FILE *fi, *fo;
	int f,ff,numFields;
	if (argc < 3)
	{
		printf("Usage: reorderVCFSubs vcf-file-with-subs-in-right-order input-vcf-file reordered-vcf-file\n\n");
		return 1;
	}
	fi = fopen(argv[1], "r");
	while (fgets(line, FIELDLENGTH*(MAXSUBS + 20), fi) && strncmp(line, "#CHROM", strlen("#CHROM")))
		;
	numFields = fillArray(line, right);
	fclose(fi);
	fi = fopen(argv[2], "r");
	fo = fopen(argv[3], "w");
	while (fgets(line, FIELDLENGTH*(MAXSUBS + 20), fi) && strncmp(line, "#CHROM", strlen("#CHROM")))
		fprintf(fo,"%s",line);
	if (numFields != fillArray(line, wrong))
	{
		dcerror(1, "Number of fields does not match between %s and %s\n", argv[1], argv[2]);
		return 1;
	}
	for (f = 0; f < numFields; ++f)
		for (ff = 0; ff < numFields; ++ff)
			if (!strcmp(wrong[f], right[ff]))
			{
				index[ff] = f;
				break;
			}
	do {
		fillArray(line,wrong);
		for (f = 0; f < numFields; ++f)
			fprintf(fo, "%s%c", wrong[index[f]], f == numFields - 1 ? '\n' : '\t');
	} while (fgets(line, FIELDLENGTH*(MAXSUBS + 20), fi)); // hopefully no blank line at the end
	fclose(fi);
	fclose(fo);
	return 0;
}