# Makefile for MD Simulation
# Author : Axel Huebl
# Version: 0.1
# Date   : Januar 12, 2012

ObjSuf        = o
SrcSuf        = cpp
OutPutOpt     = -o # keep whitespace after "-o"

# do not use more than -O1 for MPI
#CXX           = mpic++
CXX           = g++
CXXFLAGS      = -O -c -Wall
#LD            = mpic++
LD            = g++
LDFLAGS       = -O

#------------------------------------------------------------------------------

MDO       = simulation_defines.$(ObjSuf)     \
                container.$(ObjSuf)  \
                rules.$(ObjSuf)      \
                communicator.$(ObjSuf)    \
                pathfinder.$(ObjSuf) \
                mp_lib_mpi.$(ObjSuf)

MDS       = parler.$(SrcSuf)     \
                container.$(SrcSuf)  \
                rules.$(SrcSuf)      \
                communicator.$(SrcSuf)    \
                pathfinder.$(SrcSuf) \
                mp_lib_mpi.$(SrcSuf)

OBJS          = $(MDO)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

#all:            $(MDO)
all:
	$(CXX) $(CXXFLAGS) main.cpp 

clean:
	@rm -f $(OBJS) core

distclean:      clean

.SUFFIXES: .$(SrcSuf)

###

%.$(ObjSuf): %.h %.$(SrcSuf)

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
