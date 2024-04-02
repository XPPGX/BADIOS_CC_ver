#---------------------------------------------------------------------
# This file is part of BADIOS framework
# Copyright (c) 2012, 
# By:    Ahmet Erdem Sariyuce, 
#        Erik Saule,
#        Kamer Kaya,
#        Umit V. Catalyurek
#---------------------------------------------------------------------
# For license info, please see the README.txt and LICENSE.txt files in
# the main directory.
#---------------------------------------------------------------------


include makefile.in

INCLUDES   = -O2 -I./

TARGET     = bc-seq-brandes

#CFILES       = ulib 
#CXXFILES   = main-bc bc-seq-brandes

SRCS	= $(CFILES:%=%.c) $(CXXFILES:%=%.cpp)

all: $(TARGET) 

bc-seq-brandes: bucket.o main-bc.o deg1.o bridge.o artp.o idv.o sdv.o ordering.o kernels.o utl.o bc-seq-brandes.o ulib.o 
	$(LD) $(LDFLAGS) -O2 -o $@ ulib.o bucket.o main-bc.o deg1.o bridge.o artp.o idv.o sdv.o ordering.o kernels.o utl.o bc-seq-brandes.o $(LIBS)

tar:
	#/bin/rm BADIOS -r
	/bin/mkdir BADIOS; 
	/bin/cp *.cpp *.hpp *.c *.h *.txt makefile* BADIOS;
	/bin/tar -zcvf BADIOS.tar.gz BADIOS;
	/bin/rm BADIOS -r

clean: 
	/bin/rm -f *.o core *~ $(TARGET) 

depend: $(SRCS)
	makedepend -Y -- $(INCLUDES) -- $(SRCS) 

#********************************************************************
# Dependencies
#********************************************************************
# DO NOT DELETE THIS LINE -- make depend depends on it.

ulib.o: ulib.h
bc-seq-brandes.o: bc-seq-brandes.h ulib.h graph.hpp graph_impl.hpp timestamp.hpp bucket.h
#main-bc.o: bc-seq-brandes.h