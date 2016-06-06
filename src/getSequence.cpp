#include "getSequence.hpp"
#include <stdlib.h>

int main(int argc,char *argv[])
{
	faSequenceFile f;
	char buff[200];
	sprintf(buff,"\\reference\\chr%s.fa",argv[1]);
	f.init(buff);
	f.getSequence(buff,atoi(argv[2]),atoi(argv[3]));
	printf("%s\n",buff);
	return 0;
}