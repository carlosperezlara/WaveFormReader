all:
	rootcint -f Dict_maketree.cpp -c `root-config --cflags` -p LinkDef.h
	g++ -o maketree maketree.cc Aux.cc `root-config --cflags --glibs`  -rpath /Users/cperez/external/root_install/lib


#CXX = $(shell root-config --cxx)
#LD = $(shell root-config --ld)
#INC = $(shell pwd)
#Repo = $(shell git rev-parse --show-toplevel)
#SRC = $(Repo)/src
#SRC = $(Repo)

#CPPFLAGS := $(shell root-config --cflags)
#LDFLAGS := $(shell root-config --glibs) 
#CPPFLAGS += -fpermissive

#TARGET = maketree
#SRC = maketree.cc Aux.cc

#OBJ = $(SRC:.cc=.o)
#OBJ2 = $(SRC2:.cc=.o)
#OBJ3 = $(SRC3:.cc=.o)
#OBJ4 = $(SRC4:.cc=.o)

#all : $(TARGET) 

#$(TARGET) : $(OBJ)
#	@echo $@
#	$(LD) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)



#%.o : %.cc
#	@echo $@
#	$(CXX) $(CPPFLAGS) -o $@ -c $<
#clean :
#	rm -f *.o $(TARGET) *~
