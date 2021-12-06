all: vamos
LIB=libssw.so
CC=g++

vamos: vamos-beta.cpp
	$(CC) -g -fsanitize=address $< -o $@
