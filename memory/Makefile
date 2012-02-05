# Makefile for my MD Simulation
# Author : Axel Huebl
# Version: 0.1
# Date   : Januar 27, 2012

ObjSuf        = o
SrcSuf        = cpp
OutPutOpt     = -o # keep whitespace after "-o"

# do not use more than -O1 for MPI
CXX           = mpic++
#CXX           = g++
CXXFLAGS      = -O2 -c -Wall

#------------------------------------------------------------------------------

ARNAME    = $(shell basename `pwd`).ar

OBJS      = $(shell ls *.$(SrcSuf) 2>/dev/null)

#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

all:            $(OBJS)
	ar rcs $(ARNAME) $(OBJS)

clean:
	@rm -f $(OBJS) $(ARNAME)

###

#%.$(ObjSuf): %.$(SrcSuf)
#	$(CXX) $(CXXFLAGS) -c $<

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
