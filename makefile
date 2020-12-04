CURDIR = $(shell pwd)
SRCDIR = $(CURDIR)/src
INCDIR = $(CURDIR)/include
BINDIR = $(CURDIR)/bin

ROOT_LIBS = `root-config --glibs` -lSpectrum -lTreePlayer -lMathMore

LIBRS = -L$(INCDIR) $(ROOT_LIBS)# $(MINSRC) 
INCLUDE = $(INCDIR)# $(MINDIR)

CFLAGS = -std=c++11 -g -fPIC `root-config --cflags` `gsl-config --cflags` -I$(INCDIR) $(ROOT_LIBS) $(GSLLIBS) -Qunused-arguments

CPP=g++

HEAD = $(wildcard include/*.h)
OBJECTS = $(patsubst include/%.h,bin/build/%.o,$(HEAD))

TARGET = bin/libCygnus.so

main: $(TARGET)
	@printf "Make complete\n"

$(TARGET): $(OBJECTS) bin/DictOutput.cxx 
	@printf "Now compiling shared library$@\n"
	@$(CPP) $(CFLAGS) -o $@ -shared bin/DictOutput.cxx $(OBJECTS) -I. $(LIBRS) $(GSLLIBS)

bin/DictOutput.cxx: $(HEAD)
	@printf "Linking libraries\n"
	@rootcint -f $@ -c -p $(HEAD) bin/build/linkdef.h

bin/build/%.o: src/%.cxx include/%.h
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -o $@ -c $< $(LIBRS) $(GSLLIBS)
 
bin/build/%.o: src/*/%.cxx include/%.h	
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -o $@ -c $< $(LIBRS) $(GSLLIBS)

clean:  
	@printf "Tidying up...\n"
	@rm $(OBJECTS)
	@rm bin/*.*
	
