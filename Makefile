all: vamos
LIB=libssw.so
CC=g++

vamos: vamos.cpp
	$(CC) -g -fsanitize=address $< -o $@
