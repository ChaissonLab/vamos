PROG=vamos
LIBS=-lhts -labpoa -lz 
PROF=/home/cmb-16/mjc/shared/lib/
DEBUG?=""
OPT?=""
CXX=g++ -std=c++17 
CFLAGS=
asan?=""
tsan?=""

ifneq ($(DEBUG), "")
CFLAGS=-g -O0 
else
CFLAGS=-O3 -DNDEBUG #-W -Wall -pedantic -fopenmp -lpthread -lrt
endif

ifneq ($(asan), "")
  CFLAGS+=-fsanitize=address
  LIBS+=-fsanitize=address
endif

ifneq ($(tsan), "")
  CFLAGS+=-fsanitize=thread
  LIBS+=-fsanitize=thread
endif

ifneq ($(OPT), "")
STATIC=-L $(PROF) -lprofiler
endif

all:$(PROG)

vamos: main.o io.o vcf.o vntr.o
	$(CXX) $(CFLAGS) -o $@ $^ -L $(CONDA_PREFIX)/lib $(LIBS) -L ../static_lib/ -lalglib -ledlib 

main.o: main.cpp io.h vcf.h vntr.h read.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include  

io.o: io.cpp io.h vcf.h vntr.h read.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include 

vcf.o: vcf.cpp vcf.h 
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include 

vntr.o: vntr.cpp naive_anno.cpp bounded_anno.cpp io.h vntr.h read.h 
	$(CXX) $(CFLAGS) -c $< -I ../alglib/src -I ../edlib/include -I ../abpoa/include -I $(CONDA_PREFIX)/include 

# edlib.o:
# 	$(CXX) -c edlib/src/edlib.cpp -o edlib.o -I edlib/include

clean:
	rm -f $(PROG) *.o 
