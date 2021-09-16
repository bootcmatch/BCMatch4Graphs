#BHEADER**********************************************************************
# # This file is part of BCMatch4Graphs
# #
# # $release: 1.0 $
# #EHEADER**********************************************************************
TOP=..
include $(TOP)/make.inc
##################################################################
# Targets
# ##################################################################
LIBNAME=libBCM.a

all:  libs


libs:  util matvec_rule amg_rule bootamg_rule linsolvers_rule graph_rule


###################################################################
# Rules
# ##################################################################
#

bootamg_rule:
	( cd BOOTAMG; $(MAKE) LIBNAME=$(LIBNAME))

util:
	( cd utilities; $(MAKE)  LIBNAME=$(LIBNAME))

matvec_rule:
	( cd MATVEC; $(MAKE)  LIBNAME=$(LIBNAME))

amg_rule:
	( cd AMG; $(MAKE)  LIBNAME=$(LIBNAME))

linsolvers_rule:
	( cd LINSOLVERS; $(MAKE)  LIBNAME=$(LIBNAME))

graph_rule:
	( cd BCMGRAPHS; $(MAKE)  LIBNAME=$(LIBNAME))


clean: 
	( cd utilities; $(MAKE) clean) 
	( cd MATVEC; $(MAKE) clean) 
	( cd BOOTAMG;  $(MAKE) clean) 
	( cd AMG; $(MAKE) clean)
	( cd LINSOLVERS; $(MAKE) clean)
	( cd BCMGRAPHS; $(MAKE) clean)
