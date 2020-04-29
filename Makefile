# Basic Makefile

### Compilers
CC  = gcc
CXX = g++

DEBUG_LEVEL    = -g
EXTRA_CCFLAGS  = -W -Wall
CPPFLAGS       = $(DEBUG_LEVEL) $(EXTRA_CCFLAGS)
CCFLAGS        = $(CPPFLAGS)

RM = rm -f
MV = mv

SRCDIR := src
INCDIR := include

### ROOT
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

### RAT
RATLIBS  := -L$(RATROOT)/lib -lRATEvent

### BOOST
BOOSTCFLAGS := -I/usr/include/boost/
BOOSTLIBS   := -lboost_system -lboost_filesystem

CPPFLAGS  += -I$(INCDIR) $(ROOTCFLAGS) -I$(RATROOT)/include
CPPFLAGS  +=  $(BOOSTCFLAGS)
EXTRALIBS  = $(ROOTLIBS)
EXTRALIBS += $(RATLIBS)
EXTRALIBS += $(BOOSTLIBS)
EXTRALIBS += -L$(PWD)/lib -lcnpy -lz

SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJS = $(subst .cc,.o,$(SRCS))

.PHONY: all clean 
.DEFAULT_GOAL = EventWrapper

help:
	@grep -h -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: TemplateAnalysis

# TARGETLIB ?= CalibFunctions
# TARGETLIB ?= LL
# TARGETLIB ?= EVFunctions
# TARGETLIB ?= HitFunctions
# TARGETLIB ?= MCFunctions
TARGETLIB ?= AnalyzerFunctions

lib$(TARGETLIB).so: $(TARGETLIB).o
	@echo "Compiling library $@"
	$(CXX) $(CPPFLAGS) -shared $^ -o lib/$@ $(EXTRALIBS)
	$(RM) $^

$(TARGETLIB).o: 
	$(CXX) $(CPPFLAGS) -fPIC -c src/$(TARGETLIB).cc -o $@ $(EXTRALIBS)

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

VtxRecon: VtxRecon.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o VtxRecon VtxRecon.cc $(OBJS) $(EXTRALIBS)
	$(RM) VtxRecon.o $(OBJS)

PECounting: PECounting.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o PECounting PECounting.cc $(OBJS) $(EXTRALIBS)
	$(RM) PECounting.o $(OBJS)

SelectIBD: SelectIBD.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o SelectIBD SelectIBD.cc $(OBJS) $(EXTRALIBS)
	$(RM) SelectIBD.o $(OBJS)

CreatePDF: CreatePDF.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o CreatePDF CreatePDF.cc $(OBJS) $(EXTRALIBS)
	$(RM) CreatePDF.o $(OBJS)

FlattenHits: FlattenHits.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o FlattenHits FlattenHits.cc $(OBJS) $(EXTRALIBS)
	$(RM) FlattenHits.o $(OBJS)

clean:
	$(RM) $(OBJS) TemplateAnalysis TemplateAnalysisMT EventWrapper PlotCollectedPE CreateEResMatrix CreatePosMatrix PlotPEVSNHits VtxRecon SelectIBD CreatePDF FlattenHits
