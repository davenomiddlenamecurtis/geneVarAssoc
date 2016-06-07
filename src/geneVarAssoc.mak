# Makefile for geneVaraAssoc and related programs 
# the file makes its targets in ../obj and then copies executables to $DCBIN
# this means it and source files can live in ../src which will only be text files

# Destination for executables, change this if you want
DCBIN = ../bin

C = gcc
CC = g++

OURFLAGS = $(CFLAGS) 

HEADERS = btree.h consequenceType.hpp dcerror.hpp dcindex.hpp geneVarUtils.hpp getGene.hpp getSequence.hpp intervalList.h masterLocusFile.hpp vcfLocusFile.hpp hapsLocusFile.hpp

ifdef INOBJ
all: geneVarAssoc groupVarAssoc geneVarAssocAll geneVarAssocSome showAltSubs SNPVarAssoc intVarAssoc
else
all:
	if [ ! -e ../obj ] ; then mkdir ../obj ; fi ; \
	if [ ! -e ${DCBIN} ] ; then mkdir ${DCBIN} ; fi ; \
	cd ../obj; \
	make -f ../src/geneVarAssoc.mak INOBJ=INOBJ ; \
	cp geneVarAssoc groupVarAssoc geneVarAssocAll geneVarAssocSome showAltSubs SNPVarAssoc intVarAssoc ../bin ; \
	echo copied executables to /bin ; \
	cd ../src
endif
# unless you have a folder on path called $DCBIN you should delete that line to copy the executables to it

clean:
	cd ../obj ; \
	rm *.o ; \
	cd ../vcfcode

VPATH=../src
	
%.o: ../src/%.cpp $(HEADERS)
	$(CC) $(OURFLAGS) -c $< -o ../obj/$@
	
%.o: ../src/%.c $(HEADERS)
	$(C) $(OURFLAGS) -c $< -o ../obj/$@

groupGetGenotypes: groupGetGenotypes.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o outputGenotypes.o 
	$(CC) -o groupGetGenotypes groupGetGenotypes.o  btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o outputGenotypes.o 

genePhaseRec: genePhaseRec.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o genePhaseRec genePhaseRec.o  btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneTestRec: geneTestRec.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o geneTestRec geneTestRec.o  btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneVarAssoc: geneVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o geneVarAssoc geneVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

groupVarAssoc: groupVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o groupVarAssoc groupVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneVarAssocAll: geneVarAssocAll.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o geneVarAssocAll geneVarAssocAll.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

geneVarAssocSome: geneVarAssocSome.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o geneVarAssocSome geneVarAssocSome.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

intVarAssoc: intVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o intVarAssoc geneVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

SNPVarAssoc: SNPVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o SNPVarAssoc SNPVarAssoc.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

showAltSubs: showAltSubs.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o
	$(CC) -o showAltSubs showAltSubs.o btree.o consequenceType.o dcerror.o dcindex.o geneVarUtils.o getGeneFuncs.o getSequenceFuncs.o intervalList.o vcfLocusFile.o vcfWriteVars.o masterLocusFile.o hapsLocusFile.o

