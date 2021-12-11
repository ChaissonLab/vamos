PROG=vamos
LIBS=-lz -lhts
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
CFLAGS=-O0 -DNDEBUG 
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

vamos: main.o vcf.o vntr.o 
	$(CXX) $(CFLAGS) -o $@ $^ -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

main.o: main.cpp io.cpp vcf.h vntr.h read.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

vcf.o: vcf.cpp vcf.h 
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

vntr.o: vntr.cpp io.cpp vntr.h read.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -L $(CONDA_PREFIX)/lib $(LIBS) 

clean:
	rm -f $(PROG) *.o 
