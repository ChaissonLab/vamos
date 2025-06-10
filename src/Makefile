PROG = vamos
# linking static libraries
LIBS       = -lhts -lz -lpthread
#ABPOA_LIBS = -lm -lz -lpthread

ifeq ($(CC),cc)
    CC = gcc
endif
CC ?= gcc
CXX ?= g++ -std=c++11


CCFLAGS = -Wall -g2 -Wno-unused-function -Wno-misleading-indentation

DEBUG ?= "-g -std=c++11 "
OPT   ?= ""
asan  ?= ""
tsan  ?= ""

ifeq ("$(DEBUG)", "TRUE")
	  CFLAGS = -g -std=c++14 #-fsanitize=address  -static-libasan -gdwarf-4  #-Wall -Wextra
else
  CFLAGS = -O2 -DNDEBUG #-W -Wall -pedantic -fopenmp -lpthread -lrt
endif

ifneq ($(asan), "")
  CFLAGS += -fsanitize=address # boundaries check for memory access
  LIBS += -fsanitize=address
endif

ifneq ($(tsan), "")
  CFLAGS += -fsanitize=thread
  LIBS += -fsanitize=thread
endif

ifneq ($(OPT), "")
  STATIC = -L $(PROF) -ltcmalloc
endif

PROF = $(CONDA_PREFIX)/lib

ABPOA_DIR = ../abPOA-1.4.1

all: $(PROG)
ZLPATH=""
ifneq ($(ZLIB_ROOT), "")
   ZLPATH=C_INCLUDE_PATH=$(ZLIB_ROOT)/include
endif

parasail-dir/lib/libparasail.so:
	cd parasail && mkdir -p build && cd build && cmake .. -DPARASAIL_ENABLE_TESTS=OFF -DCMAKE_INSTALL_PREFIX=$PWD/../../parasail-dir/ && make -j4 && make install

$(ABPOA_DIR)/lib/libabpoa.a:
	echo $(CC)
	echo "^^CC"
	cd $(ABPOA_DIR) && make sse2=true PREFIX=. CC="$(CC)" $(ZLPATH) lib/libabpoa.a

testmsa: TestMSA.cpp
	$(CXX) $(CFLAGS) -o $@ $^ -L $(CONDA_PREFIX)/lib $(LIBS) -I $(ABPOA_DIR)/include -L $(ABPOA_DIR)/lib -labpoa -lrt -L$(mchaisso)/software/lib

OBJS=main.o io.o vcf.o vntr.o acc_lookup_table.o phase.o msa.o

vamos:  $(OBJS) $(ABPOA_DIR)/lib/libabpoa.a parasail-dir/lib/libparasail.so
	$(CXX) $(CFLAGS) -o $@ $(OBJS) -L $(CONDA_PREFIX)/lib $(LIBS) -L $(ABPOA_DIR)/lib -labpoa -Lparasail-dir/lib -lparasail -Wl,-rpath=$(PWD)/parasail-dir/lib 

main.o: main.cpp io.h vcf.h vntr.h read.h option.h threads.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include -I parasail-dir/include 

io.o: io.cpp io.h vcf.h vntr.h read.h phase.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include  -I parasail-dir/include -I $(ABPOA_DIR)/include -I $(CONDA_PREFIX)/include

phase.o: phase.cpp phase.h read.h
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include

vcf.o: vcf.cpp vcf.h 
	$(CXX) $(CFLAGS) -c $< -I $(CONDA_PREFIX)/include

vntr.o: vntr.cpp bounded_anno_local.cpp io.h vntr.h read.h option.h msa.h
	$(CXX) $(CFLAGS) -c $< -I $(ABPOA_DIR)/include -I ../alglib/src -I ../edlib/include -I $(CONDA_PREFIX)/include

msa.o: msa.cpp msa.h read.h option.h
	$(CXX) $(CFLAGS) -c $< -I $(ABPOA_DIR)/include

acc_lookup_table.o: acc_lookup_table.cpp acc_lookup_table.h vntr.cpp vntr.h
	$(CXX) $(CFLAGS) -c $< -I ../alglib/src -I ../edlib/include -I $(ABPOA_DIR)/include  -I $(CONDA_PREFIX)/include

# edlib.o:
# 	$(CXX) -c edlib/src/edlib.cpp -o edlib.o -I edlib/include

clean:
	rm -f $(PROG) *.o $(ABPOA_SRC)/*.o $(ABPOA_DIR)/lib/libabpoa.a
