CC = g++ # -g # the -g can be removed it's for debug 
VERSION = -std=gnu++14
OPT = -O3 # -O3 # optimisation de ouf
OBJECTS = main.cpp
CFLAGS = -std=gnu++14 -pedantic -Wall

SRC=$(wildcard *.cpp) # the list of source files 
OBJS = $(SRC:.cpp=.o) # the list of object files (les .cpp sont convertis en .o)
HEAD=$(filter-out main.hpp  , $(SRC:.cpp=.hpp))  #converti tout les .cpp en .h sauf main.cpp et certains autres...

MAIN = Sys


all: $(MAIN)


$(MAIN) : $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS)
# creates the executable from the .o files 


main.o: $(SRC) $(HEAD)



.cpp.o:
	$(CC) $(VERSION) -fPIC $(OPT) -c $< -o $@ $(CFLAGS)
# $(CC) -fPIC $(OPT) -o $@ -c $< $(CFLAGS)
# generic rule enabling to create the .o files from .cpp files
# $@ refers to the name of the cible 
# $< refers to the name of the fist dependency 


.PHONY: clean 

cleanmake:
	rm -f *.o 
# suppress all .o files 






