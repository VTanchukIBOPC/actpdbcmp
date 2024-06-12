all: actpdbcmp

actpdbcmp: actpdbcmp.o  Cluster2.o Compare.o pdb.o Sobseq.o Nrutil.o OptimBase.o PowellOptim.o
	g++ actpdbcmp.o  Cluster2.o Compare.o pdb.o Sobseq.o Nrutil.o OptimBase.o PowellOptim.o -o actpdbcmp

actpdbcmp.o: actpdbcmp.cpp
	g++ -c actpdbcmp.cpp
Cluster2.o: Cluster2.cpp
	g++ -c Cluster2.cpp
Compare.o: Compare.cpp
	g++ -c Compare.cpp
pdb.o: pdb.cpp
	g++ -c pdb.cpp
Sobseq.o: Sobseq.cpp
	g++ -c Sobseq.cpp
Nrutil.o: Nrutil.cpp
	g++ -c Nrutil.cpp
OptimBase.o: Nrutil.cpp
	g++ -c OptimBase.cpp
PowellOptim.o: PowellOptim.cpp
	g++ -c PowellOptim.cpp
clean:
	rm -rf *o actpdbcmp

