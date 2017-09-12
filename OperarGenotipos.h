#ifndef OPERARGENOTIPOS_H_
#define OPERARGENOTIPOS_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <bitset>
#include <limits>
#include <ctime>
#include <cstdlib>

#include "MersenneTwister.h"

using namespace std;

double randomUnif ();
double randomNormal(double media, double varianza);
void frecuenciasLoci(float frecuencia, float desviacionEstandar,int numeroGenes, float *frecs);
void frecuenciasLociUniforme(int numeroGenes, float *frecs);

// Convierte un número entero indicando el valor del genotipo de un gameto
// en un vector de booleanos

template <typename T> std::string xxx_to_bin(const T& value,
		const int posiciones);

void convertirIntEnBool(int genotipoInt, int tamanyoHaploide, bool *vectorbool);
int binstring2int(string str);
// Convierte un entero indicando un genotipo en un string de 1 y 0 de 2n
string convertirIntEnStringDiploide(int id, int nGenes);

// Convierte un entero indicando un genotipo en un string de 1 y 0 de n
string convertirIntEnStringHaploide(int id, int nGenes);

// Acumula las probabibilidades gaméticas y manda los ceros al final de la
// matriz
//int ordenarAcumular(double **a, int size);
int ordenar(double ***a, int size);
void shearSort(double ***temp, int nCeros, int size);
void sortPart1(double ***temp , int Lo, int Hi, int Nx,bool Up);
void sortPart2(double ***temp , int Lo, int Hi, int Nx,bool Up);
void acumularProb(double ***a, double ***temp, int nCeros, int size);

void ordenarEclosionados(long *a, int size);
void shearSortLong(long *temp, int size);
void sortPart1Long(long *temp, int Lo, int Hi, int Nx,bool Up);
void sortPart2Long(long *temp, int Lo, int Hi, int Nx,bool Up);

// Busqueda Binaria en un array double
int binarySearch(double **list, int size, double searchItem, int ultimo);


#endif /*OPERARGENOTIPOS_H_*/
