CC = gcc-5
CXX = g++-5

SRCDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp

TARGET = $@

COMPILE_FLAGS = -std=c++11
RCOMPILE_FLAGS = -D NDEBUG -O3 -fopenmp
DCOMPILE_FLAGS = -D DEBUG -g
LINK_FLAGS = -gnu 
RLINK_FLAGS = -fopenmp
DLINK_FLAGS = 

INCLUDES = -I/opt/X11/include -I/usr/X11R6/include -I/usr/include -I/usr/local/include -I/usr/local/include/gsl
GLXLIBS = -L/opt/X11/lib -lglfw3 -framework OpenGL -lglew -lpthread
GSLLIBS = -lgsl -lgslcblas
FFTLIBS = -L/usr/lib64 -lfftw3
LIBS = -L/usr/local/lib $(GLXLIBS) $(GSLLIBS) $(FFTLIBS) -lyaml-cpp

UNAME_S:=$(shell uname -s)

print-%: ; @echo $*=$($*)

SHELL = /bin/bash

.SUFFIXES:

ifneq ($(LIBS),)
	COMPILE_FLAGS += $(shell pkg-config --cflags $(LIBS))
	LINK_FLAGS += $(shell pkg-config --libs $(LIBS))
endif

# Special stuff for intel compiler
CC=$(CXX)
ifeq ($(CC),icpc)
	COMPILE_FLAGS += -Wno-deprecated
	RCOMPILE_FLAGS += -openmp -DBOB_OMP
else
	COMPILE_FLAGS += -Wno-deprecated-declarations -Wno-deprecated
	RCOMPILE_FLAGS += -fopenmp -DBOB_OMP
endif

# Combine compiler and linker flags
ifeq ($(CFG),release)
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(RCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(RLINK_FLAGS)
else
	export CXXFLAGS := $(CXXFLAGS) $(COMPILE_FLAGS) $(DCOMPILE_FLAGS)
	export LDFLAGS := $(LDFLAGS) $(LINK_FLAGS) $(DLINK_FLAGS)
endif

# build information on all sources
ifeq ($(UNAME_S),Darwin)
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' | sort -k 1nr | cut -f2-)
else
	SOURCES = $(shell find $(SRCDIR) -name '*.$(SRCEXT)' -printf '%T@\t%p\n' \
			  										| sort -k 1nr | cut -f2-)
endif

# fallback case
rwildcard = $(foreach d, $(wildcard $1*), $(call rwildcard,$d/,$2) \
									$(filter $(subst *,%,$2), $d))
ifeq ($(SOURCES),)
	SOURCES := $(call rwildcard, $(SRCDIR), *.$(SRCEXT))
endif

# Now we have to figure out which we are building of the program list, since that matters
# for things like not compiling more than one main
MAIN_SOURCES = $(SRCDIR)/cytoscore_main.cpp $(SRCDIR)/configure_cytoscore.cpp

SRCS = $(filter-out $(MAIN_SOURCES), $(SOURCES))
ifeq ($(TARGET),cytoscore)
	SRCS += $(SRCDIR)/cytoscore_main.cpp
else ifeq ($(TARGET),configure_cytoscore)
	SRCS += $(SRCDIR)/configure_cytoscore.cpp
endif

OBJECTS = $(SRCS:$(SRCDIR)/%.$(SRCEXT)=$(OBJDIR)/%.o)
DEPS = $(OBJECTS:.o=.d)

.PHONY: dirs
dirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BINDIR)

.PHONY: clean
clean:
	$(RM) -r $(OBJDIR)
	$(RM) -r $(BINDIR)

cytoscore: dirs $(BINDIR)/$(TARGET)
	@echo "Bulding $(BINDIR)/$(TARGET)"

$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LDFLAGS) $(LIBS) -o $@

# add dependencies
-include $(DEPS)

# source file rules
$(OBJDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MP -MMD -c $< -o $@


