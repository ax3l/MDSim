# Makefile for my MD Simulation
# Author : Axel Huebl
# Version: 0.1
# Date   : Januar 12, 2012

ObjSuf        = o
SrcSuf        = cpp
LibSuf        = ar
OutPutOpt     = -o # keep whitespace after "-o"

# do not use more than -O1 for MPI
CXX           = mpic++
#CXX           = g++
CXXFLAGS      = -O2 -c -Wall
LD            = mpic++
#LD            = g++
LDFLAGS       = -O2

#------------------------------------------------------------------------------

OBJS      = 
LIBS      = memory/memory.ar             \
            physics/physics.ar           \
            communicator/communicator.ar 

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(LibSuf)

all:            $(OBJS) $(LIBS)
	$(LD) $(LDFLAGS) main.cpp $(LIBS)

clean:
	@rm -f $(OBJS) core *.o
	for i in $(LIBS); do cd `dirname $$i` && $(MAKE) clean && cd -; done

distclean:      clean

#.SUFFIXES: .$(SrcSuf)

###

#%.$(ObjSuf): %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -c $<

$(LIBS):
	cd `dirname $@` && $(MAKE)

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
