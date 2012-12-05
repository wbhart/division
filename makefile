CC=gcc
CFLAGS=-O2
INCLUDE=-I. -I/home/wbhart/mpir-2.6.0
LIBS=-L/home/wbhart/mpir-2.6.0/.libs

all: division.o t-division.c
	$(CC) $(CFLAGS) t-division.c division.o -o t-division $(INCLUDE) $(LIBS) -static -lmpir

division.o: division.c division.h
	$(CC) $(CFLAGS) -c division.c -o division.o $(INCLUDE)

.PHONY: all
