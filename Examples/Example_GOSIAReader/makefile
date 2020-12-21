LIBS = -lCygnus -lSpectrum -lTreePlayer -lMathMore
CPP = g++
CFLAGS = -std=c++0x `root-config --cflags --libs --glibs`
INCDIR = $(shell pwd)/../../include
LIBDIR = $(shell pwd)/../../bin
CURDIR = $(shell pwd)
BINDIR = $(CURDIR)

SRCS := $(wildcard *.cxx)
EXES := $(patsubst %.cxx,$(BINDIR)/%,$(SRCS))
#OBJS := $(patsubst %,%.o,$(EXES))

all : $(EXES)

$(BINDIR)/%: $(CURDIR)/%.cxx
	@mkdir -p $(BINDIR)
	@printf "Creating executable $@ from $<\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -L$(LIBDIR) $(LIBS) $< -o $@

clean:
	@printf "Cleaning up\n" 
	@rm $(EXES)
