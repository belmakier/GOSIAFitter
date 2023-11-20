CURDIR = $(shell pwd)
SRCDIR = $(CURDIR)/src
INCDIR = $(CURDIR)/include
BINDIR = $(CURDIR)/bin

ROOT_LIBS = `root-config --glibs` -lSpectrum -lTreePlayer -lMathMore
FORTRAN_LIBS = -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/ -lgfortran

LIBRS = $(ROOT_LIBS) $(GSLLIBS) $(FORTRAN_LIBS)
INCLUDE = $(INCDIR)

CFLAGS = -std=c++17 -g -fPIC `root-config --cflags` `gsl-config --cflags` -Wno-unused-parameter

PLATFORM:=$(shell uname)
$(info PLATFORM: $(PLATFORM))
 
ifeq ($(PLATFORM),Darwin)
export __APPLE__:= 1
CFLAGS     += -Qunused-arguments
CPP        = clang++
FORTRAN_LIBS = -L/opt/local/lib/gcc12  -lgfortran
else
export __LINUX__:= 1
CPP        = g++
endif

HEAD = $(wildcard include/*.h)
OBJECTS = $(patsubst include/%.h,lib/%.so,$(HEAD))

$(info   OBJECTS = $(OBJECTS))

TARGET = bin/libGOSIAFitter.so

main: $(TARGET) bin/gosia 
	@printf "Make complete\n"

$(TARGET): $(OBJECTS) bin/DictOutput.cxx lib bin obj/gosia.o
	@printf "Now compiling shared library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -I. -o $@ -shared bin/DictOutput.cxx $(OBJECTS) obj/gosia.o $(LIBRS)

bin/DictOutput.cxx: $(HEAD)
	@printf "Linking libraries\n"
	@rootcint -f $@ -c -p $(HEAD) lib/linkdef.h

lib bin:
	@mkdir -p bin lib

obj/gosia.o: src/gosia_20231120.1.f include/Gosia.h
	@printf "Now compiling object gosia.o\n"
	@gfortran -fPIC -std=legacy -m64 -g -o obj/gosia.o -c src/gosia_20231120.1.f $(FORTRAN_LIBS)

lib/%.so: src/%.cxx include/%.h 
	@printf "Now compiling library $@\n"
	@$(CPP) $(CFLAGS) -I$(INCDIR) -o $@ -shared -c $< $(LIBRS) 

bin/gosia: bin src/gosia_20081208.18.f
	@printf "Compiling GOSIA \n" 
	@gfortran src/gosia_20081208.18.f -o bin/gosia 

clean:  
	@printf "Tidying up...\n"
	@rm $(OBJECTS)
	@rm -r bin/*
	@rm -r obj/*

install:
	cp bin/libGOSIAFitter.so $(HOME)/.local/lib/
	rm -r $(HOME)/.local/include/GOSIAFitter/
	mkdir $(HOME)/.local/include/GOSIAFitter/
	cp include/*.h $(HOME)/.local/include/GOSIAFitter/
