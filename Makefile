# Basic Makefile

### Compilers
CC  = gcc
CXX = g++

DEBUG_LEVEL    = -g
EXTRA_CCFLAGS  = -W -Wall -std=c++11
CPPFLAGS       = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS        = $(CPPFLAGS)

RM = rm -f

SRCDIR := src
INCDIR := include

### ROOT
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

### RAT
RATROOT := /home/zsoldos/theia/rat-pac
RATLIBS  := -L$(RATROOT)/lib -lRATEvent

### BOOST
BOOSTCFLAGS := -I/usr/include/boost/
BOOSTLIBS   := -lboost_system -lboost_filesystem

CPPFLAGS  += -I$(ROOTSYS)/include -I$(INCDIR) $(ROOTCFLAGS) -I$(RATROOT)/include
CPPFLAGS  +=  $(BOOSTCFLAGS)
EXTRALIBS  = $(ROOTLIBS)
EXTRALIBS += $(RATLIBS)
EXTRALIBS += $(BOOSTLIBS)

SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJS = $(subst .cc,.o,$(SRCS))

.PHONY: all clean 
.DEFAULT_GOAL = EventWrapper

help:
	@grep -h -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: TemplateAnalysis

TemplateAnalysis: TemplateAnalysis.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o TemplateAnalysis TemplateAnalysis.cc $(OBJS) $(EXTRALIBS)
	$(RM) TemplateAnalysis.o $(OBJS)

TemplateAnalysisMT: TemplateAnalysisMT.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o TemplateAnalysisMT TemplateAnalysisMT.cc $(OBJS) $(EXTRALIBS)
	$(RM) TemplateAnalysisMT.o $(OBJS)

EventWrapper: pyevent.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o EventWrapper pyevent.cpp $(OBJS) $(EXTRALIBS)
	$(RM) pyevent.o $(OBJS)

PlotCollectedPE: PlotCollectedPE.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o PlotCollectedPE PlotCollectedPE.cc $(OBJS) $(EXTRALIBS)
	$(RM) PlotCollectedPE.o $(OBJS)

CreateEResMatrix: CreateEResMatrix.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o CreateEResMatrix CreateEResMatrix.cc $(OBJS) $(EXTRALIBS)
	$(RM) CreateEResMatrix.o $(OBJS)

CreatePosMatrix: CreatePosMatrix.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o CreatePosMatrix CreatePosMatrix.cc $(OBJS) $(EXTRALIBS)
	$(RM) CreatePosMatrix.o $(OBJS)

PlotPEVSNHits: PlotPEVSNHits.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o PlotPEVSNHits PlotPEVSNHits.cc $(OBJS) $(EXTRALIBS)
	$(RM) PlotPEVSNHits.o $(OBJS)

clean:
	$(RM) $(OBJS) TemplateAnalysis TemplateAnalysisMT EventWrapper PlotCollectedPE CreateEResMatrix CreatePosMatrix PlotPEVSNHits
