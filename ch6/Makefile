# Simple makefile to build practice solutions
#http://mrbook.org/blog/tutorials/make/ 
CC=clang++
CFLAGS=-g -O2 -std=c++11 -Wall -pedantic
INC=../../mottuInc
all: dynProg

dynProg: main.cpp
	$(CC) $(CFLAGS) main.cpp -o dynProg -I$(INC)

examples: examples.cpp
	$(CC) $(CFLAGS) examples.cpp -o examples -I$(INC)

#practice: main.o Solution.o
#	$(CC) main.o Solution.o -o practice

#Solution.o: Solution.cpp
#	$(CC) $(CFLAGS) Solution.cpp

#main.o: main.cpp
#	$(CC) $(CFLAGS) main.cpp

clean:
	rm *.o dynProg examples
 

