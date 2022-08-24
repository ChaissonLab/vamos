PROG=vamos
LIBS=-lhts  -lz -lpthread 
PROF=/home1/jingwenr/.conda/envs/trf/lib
DEBUG?="-g "
OPT?=""
CXX=g++ -std=c++17 
CFLAGS=
asan?=""
tsan?=""
ABPOA=abpoa/
ifeq ("$(DEBUG)", "TRUE")
CFLAGS=-g  #-fsanitize=address #-Wall -Wextra
else
CFLAGS=-O3 -DNDEBUG #-W -Wall -pedantic -fopenmp -lpthread -lrt
endif

ifneq ($(asan), "")
  CFLAGS+=-fsanitize=address # boundaries check for memory access
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

abpoa/lib/libabpoa.a:
	pushd abpoa && make -j 4 CC=gcc

vamos: main.o io.o vcf.o vntr.o acc_lookup_table.o abpoa/lib/libabpoa.a
	$(CXX) $(CFLAGS) -o $@ $^ -L$(ABPOA)/lib -labpoa -L $(CONDA_PREFIX)/lib $(LIBS) -L ../static_lib/ -lalglib -ledlib -lrt 

main.o: main.cpp io.h vcf.h vntr.h read.h option.h threads.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include  

io.o: io.cpp io.h vcf.h vntr.h read.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -I ../abpoa/include -I $(CONDA_PREFIX)/include 

vcf.o: vcf.cpp vcf.h 
	$(CXX) $(CFLAGS) -c $<  -I $(CONDA_PREFIX)/include 

phase.o: phase.cpp phase.h read.h
	$(CXX) $(CFLAGS) -c $<  -I $(CONDA_PREFIX)/include

vntr.o: vntr.cpp naive_anno.cpp bounded_anno.cpp bounded_anno_local.cpp io.h vntr.h read.h option.h
	$(CXX) $(CFLAGS) -c $< -I $(ABPOA)/include -I ../alglib/src -I ../edlib/include -I -I $(CONDA_PREFIX)/include

acc_lookup_table.o: acc_lookup_table.cpp acc_lookup_table.h vntr.cpp vntr.h
	$(CXX) $(CFLAGS) -c $< -I ../alglib/src -I ../edlib/include -I ../abpoa/include -I $(CONDA_PREFIX)/include

# edlib.o:
# 	$(CXX) -c edlib/src/edlib.cpp -o edlib.o -I edlib/include

clean:
	rm -f $(PROG) *.o 
