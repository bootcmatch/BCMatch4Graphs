#BHEADER**********************************************************************
# # This file is part of BCMatch4Graphs.
# #
# # $release: 1.0 $
# #EHEADER**********************************************************************
include ./make.inc
##################################################################
# Targets
# ##################################################################
LIBNAME=libBCM4Graphs.a

all: libdir includedir libs

libdir:
	(if test ! -d lib ; then mkdir lib; fi)
includedir:
	(if test ! -d include ; then mkdir include; fi; )

libs: src_rule 
	/bin/mv $(LIBNAME) lib

###################################################################
# Rules
# ##################################################################
#

src_rule:
	( cd SRC; $(MAKE) LIBNAME=$(LIBNAME))

tests:
	( cd Test; $(MAKE) )

clean:
	( cd include; /bin/rm -f *.h)
	( cd lib; /bin/rm -f *.a)
	( cd SRC; $(MAKE) clean)
	( cd Test; $(MAKE) clean)
