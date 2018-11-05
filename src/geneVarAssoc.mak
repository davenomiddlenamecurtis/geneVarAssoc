# Makefile for geneVaraAssoc and related programs 
# the file makes its targets in ../obj and then copies executables to $DCBIN
# this means it and source files can live in ../src which will only be text files

# Destination for executables, change this if you want
DCBIN = ../bin

C = gcc
CC = g++
MAXVCFFILES = 10
MAXSUB = 20000
CFLAGS = $(OURFLAGS) $(DEBUGFLAG)

HEADERS = consequenceType.hpp dcerror.hpp dcindex.hpp geneVarUtils.hpp getGene.hpp getSequence.hpp intervalList.h masterLocusFile.hpp vcfLocusFile.hpp hapsLocusFile.hpp geneVarParser.hpp dcexpr.hpp
EXE = geneVarAssoc intVarAssoc
ifdef INOBJ
all: ${EXE}
else
all:
	if [ ! -e ../obj ] ; then mkdir ../obj ; fi ; \
	if [ ! -e ${DCBIN} ] ; then mkdir ${DCBIN} ; fi ; \
	cd ../obj; \
	make -f ../src/geneVarAssoc.mak INOBJ=INOBJ CFLAGS=$(CFLAGS) MAXSUB=${MAXSUB} MAXVCFFILES=${MAXVCFFILES} ; \
	cp ${EXE} ${DCBIN} ; \
	echo copied executables to ${DCBIN} ; \
	cd ../src
endif

clean:
	cd ../obj ; \
	rm *.o ; \
	cd ../src

VPATH=../src
	
%.o: ../src/%.cpp $(HEADERS)
	$(CC) $(CFLAGS) -DMAXSUB=${MAXSUB} -DMAXVCFFILES=${MAXVCFFILES} -c $< -o ../obj/$@
	
%.o: ../src/%.c $(HEADERS)
	$(C) $(CFLAGS) -DMAXSUB=${MAXSUB} -DMAXVCFFILES=${MAXVCFFILES} -c $< -o ../obj/$@

groupGetGenotypes: groupGetGenotypes.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o outputGenotypes.o 
	$(CC) -o groupGetGenotypes groupGetGenotypes.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o outputGenotypes.o 

genePhaseRec: genePhaseRec.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o genePhaseRec genePhaseRec.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneTestRec: geneTestRec.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o geneTestRec geneTestRec.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneVarAssoc: geneVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o geneVarParser.o dcexpr.o
	$(CC) -o geneVarAssoc geneVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o geneVarParser.o dcexpr.o

groupVarAssoc: groupVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o groupVarAssoc groupVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneVarAssocAll: geneVarAssocAll.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o geneVarAssocAll geneVarAssocAll.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneVarAssocSome: geneVarAssocSome.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o geneVarAssocSome geneVarAssocSome.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

intVarAssoc: intVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o geneVarParser.o dcexpr.o
	$(CC) -o intVarAssoc intVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o geneVarParser.o dcexpr.o

SNPVarAssoc: SNPVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o SNPVarAssoc SNPVarAssoc.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

showAltSubs: showAltSubs.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o geneVarParser.o dcexpr.o
	$(CC) -o showAltSubs showAltSubs.o  consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o geneVarParser.o dcexpr.o

