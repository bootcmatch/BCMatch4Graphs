#BHEADER**********************************************************************
# # This file is part of BCMatch4Graphs.
# # #
# # # $release: 1.0 $
# # #EHEADER**********************************************************************
#
TOP=../..
include $(TOP)/make.inc

HEADERS =\
 bcm_graphs.h\
 bcm_metadata_exp.h


FILES =\
 bcm_kmeans.c\
 bcm_modularity.c\
 bcm_entropy.c\
 bcm_metadata.c\
 bcm_kmeansblock.c\
 bcm_augmentation.c
 
INCDIR=$(TOP)/include
CINCLUDES += -I$(INCDIR)
OBJS = $(FILES:.c=.o)

##################################################################
# Targets
##################################################################

all: lib copy

copy:
	/bin/cp -p $(HEADERS) $(INCDIR)

clean:
	rm -f *.o 

##################################################################
# Rules
##################################################################

lib: $(OBJS)
	ar cur $(TOP)/$(LIBNAME) $(OBJS)
	ranlib $(TOP)/$(LIBNAME)

$(OBJS) $(FILES): $(HEADERS)
