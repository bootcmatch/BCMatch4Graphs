#BHEADER**********************************************************************
# # This file is part of BCMatch4Graphs.
# # #
# # # $release: 1.0 $
# # #EHEADER**********************************************************************
include ../make.inc
#
TOP=..
INCDIR=$(TOP)/include
LIBDIR=$(TOP)/lib
CINCLUDES +=-I$(INCDIR)

RUNDIR=./RUNS

all: test

test: clustering_AttributedGraph clustering_Graph 

clean:
	rm -f *.o $(RUNDIR)/clustering_Graph $(RUNDIR)/clustering_AttributedGraph 

##################################################################
# Rules
##################################################################
UTILS=ioutil.o

clustering_Graph:  clustering_Graph.o $(UTILS)
	$(LINKER) $(LINKOPT) -o $@ $<  $(UTILS) -lBCM4Graphs -L$(LIBDIR) $(LIBS) -lm $(LINKOPT)
	mv $@ $(RUNDIR)

clustering_AttributedGraph:  clustering_AttributedGraph.o $(UTILS)
	$(LINKER) $(LINKOPT) -o $@ $<  $(UTILS) -lBCM4Graphs -L$(LIBDIR) $(LIBS) -lm $(LINKOPT)
	mv $@ $(RUNDIR)
