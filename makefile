CURDIR = $(shell pwd)
SRCDIR = $(CURDIR)/src
INCDIR = $(CURDIR)/include
BINDIR = $(CURDIR)/bin

ROOT_LIBS = `root-config --glibs` -lSpectrum -lTreePlayer -lMathMore

LIBRS = $(INCDIR) $(ROOT_LIBS) $(GSLLIBS)
INCLUDE = $(INCDIR)

CFLAGS = -std=c++11 -g -fPIC `root-config --cflags` `gsl-config --cflags` -Wno-unused-parameter

PLATFORM:=$(shell uname)
$(info PLATFORM: $(PLATFORM))
 
ifeq ($(PLATFORM),Darwin)
export __APPLE__:= 1
CFLAGS     += -Qunused-arguments
CPP        = clang++
else
export __LINUX__:= 1
CPP        = g++
endif

HEAD = $(wildcard include/*.h)
OBJECTS = $(patsubst include/%.h,bin/build/%.o,$(HEAD))

TARGET = bin/libCygnus.so

main: $(TARGET)
	@printf "Make complete\n"

$(TARGET): $(OBJECTS) bin/DictOutput.cxx 
	@printf "Now compiling shared library$@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -I. -L$(LIBRS) -o $@ -shared bin/DictOutput.cxx $(OBJECTS) 

bin/DictOutput.cxx: $(HEAD)
	@printf "Linking libraries\n"
	@rootcint -f $@ -c -p $(HEAD) bin/build/linkdef.h

bin/build/%.o: src/%.cxx include/%.h
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -L$(LIBRS) -o $@ -c $<
 
bin/build/%.o: src/*/%.cxx include/%.h	
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -L$(LIBRS) -o $@ -c $< 

clean:  
	@printf "Tidying up...\n"
	@rm $(OBJECTS)
	@rm bin/*.*
	
