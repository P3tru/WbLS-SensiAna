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
RATROOT := /home/zsoldos/theia/rat-pac-chess
RATLIBS  := -L$(RATROOT)/lib -lRATEvent

CPPFLAGS  += -I$(ROOTSYS)/include -I$(INCDIR) $(ROOTCFLAGS) -I$(RATROOT)/include
EXTRALIBS  = $(ROOTLIBS)
EXTRALIBS += $(RATLIBS)

SRCS = $(wildcard $(SRCDIR)/*.cc)
OBJS = $(subst .cc,.o,$(SRCS))

.PHONY: all clean 
.DEFAULT_GOAL = EventWrapper

help:
	@grep -h -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: EventWrapper

EventWrapper: pyevent.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o EventWrapper pyevent.cpp $(OBJS) $(EXTRALIBS)
	$(RM) pyevent.o $(OBJS)

PlotCollectedPE: PlotCollectedPE.o $(OBJS)
	$(CXX) $(CPPFLAGS) -o PlotCollectedPE PlotCollectedPE.cc $(OBJS) $(EXTRALIBS)
	$(RM) PlotCollectedPE.o $(OBJS)

clean:
	$(RM) $(OBJS) EventWrapper PlotCollectedPE
