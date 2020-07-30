CC = g++
INCLUDES= $(shell find -name "*.h")
ROOTCONFIG = `root-config --cflags`
ROOTGLIBS = `root-config --glibs`
ROOTLIBS = `root-config --libs` -lMathMore
GSLCONFIG = `gsl-config --cflags`
GSLLIBS = `gsl-config --libs`
SOURCES = $(shell find -name "*.cpp")

all: test

test:
	$(CC) $(INCLUDES) $(ROOTCONFIG) $(GSLCONFIG) $(SOURCES) -o test $(ROOTGLIBS) $(ROOTLIBS) $(GSLLIBS)

clean:
	$(RM) test
