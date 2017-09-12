
all: run_founder

run_founder: main.o IndicesGeneticos.o  Laguna.o Poblacion.o Genotipo.o BancoHuevos.o OperarGenotipos.o
	g++ -g -o run_founder main.o IndicesGeneticos.o Laguna.o Poblacion.o Genotipo.o BancoHuevos.o OperarGenotipos.o

main.o: main.cpp MersenneTwister.h
	g++ -c -Wall main.cpp MersenneTwister.h

IndicesGeneticos.o: IndicesGeneticos.cpp IndicesGeneticos.h
	g++ -c -Wall IndicesGeneticos.cpp IndicesGeneticos.h

Laguna.o: Laguna.cpp Laguna.h
	g++ -c -Wall Laguna.cpp Laguna.h

Poblacion.o: Poblacion.cpp Poblacion.h
	g++ -c -Wall Poblacion.cpp Poblacion.h

Genotipo.o: Genotipo.cpp Genotipo.h
	g++ -c -Wall Genotipo.cpp Genotipo.h

BancoHuevos.o: BancoHuevos.cpp BancoHuevos.h
	g++ -c -Wall BancoHuevos.cpp BancoHuevos.h

OperarGenotipos.o: OperarGenotipos.cpp OperarGenotipos.h
	g++ -c -Wall OperarGenotipos.cpp OperarGenotipos.h

clean:
	rm -rf *.o *.gch *.out
