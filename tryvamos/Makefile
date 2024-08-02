all: vamossomatic

vamossomatic: vamossomatic.cpp
	g++ -O2 vamossomatic.cpp -o vamossomatic -I $(CONDA_PREFIX)/include -lz -L $(CONDA_PREFIX)/lib -lboost_iostreams
