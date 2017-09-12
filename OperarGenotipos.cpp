#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <bitset>
#include <limits>
#include <ctime>
#include <cstdlib>

#include "OperarGenotipos.h"
#include "MersenneTwister.h"
#include "Laguna.h"

using namespace std;

double randomUnif() {
	double aleat = (double) rand() / RAND_MAX;
	while (aleat > 1)
		aleat = (double) rand() / RAND_MAX;
	return aleat;
}

// Establece las frecuencias al azar de los loci a partir de una frecuencia
// media y su desviacion estandar siguiendo una distribucion normal
void frecuenciasLoci(float media, float desviacionEstandar, int numeroGenes,
		float *frecs) {
	MTRand rnd;
	float temp;
	for (int i = 0; i < numeroGenes; i++) {
		temp = (float) rnd.randNorm(media, desviacionEstandar);
		if ((temp < 0) || (temp > 1))
			i--;
		else
			frecs[i] = temp;
	}
}

// Establece las frecuencias al azar siguiendo una distribucion uniforme
void frecuenciasLociUniforme(int numeroGenes, float *frecs) {
	MTRand rnd;
	for (int i = 0; i < numeroGenes; i++) {
		frecs[i] = (float) rnd.rand();
		//		cout<<frecs[i]<<" "<<endl;
	}
}

// Convierte un número entero indicando el valor del genotipo de un gameto
// en un vector de booleanos

template<typename T> std::string xxx_to_bin(const T& value,
		const int posiciones) {
	string mystr;
	const std::bitset<std::numeric_limits<T>::digits + 1> bs(value);
	int tam = std::numeric_limits<T>::digits + 1;
	if (tam < posiciones) {
		mystr = string('0', tam - posiciones) + bs.to_string();
	} else if (tam > posiciones) {
		mystr = bs.to_string();
		mystr.erase(0, tam - posiciones);
	} else
		mystr = bs.to_string();
	return (mystr);
}

int binstring2int(string str) {
	bitset<50> bs(str);
	return (int) bs.to_ulong();

}
void convertirIntEnBool(int genotipoInt, int tamanyoHaploide, bool *genotipo) {
	string val = xxx_to_bin(genotipoInt, tamanyoHaploide);
	for (unsigned int i = 0; i < val.length(); i++) {
		genotipo[i] = (val[i] == '0') ? false : true;
	}
}

// Convierte un entero indicando un genotipo en un string de 1 y 0 de 2n
string convertirIntEnStringDiploide(int id, int nGenes) {
	string combinacion;
	combinacion = xxx_to_bin(id, 2 * nGenes);
	return combinacion;
}

// Convierte un entero indicando un genotipo en un string de 1 y 0 de n
string convertirIntEnStringHaploide(int id, int nGenes) {
	string combinacion = xxx_to_bin(id, nGenes);
	return combinacion;
}

// Acumula las probabibilidades gameticas y manda los ceros al final de la matriz
int ordenar(double ***a, int size) {

	int nCeros = 0;
	double **temp = new double *[size];
	for (int i = 0; i < size; i++) {
		temp[i] = new double[2];
		temp[i][0] = 0;
		temp[i][1] = i;
	}

	for (int i = 0; i < size; i++) {
		if (((*a)[i][0]) != 0) {
			temp[i - nCeros][0] = (*a)[i][0];
			temp[i - nCeros][1] = (*a)[i][1];
		} else {
			nCeros++;
		}
	}

	shearSort(&temp, nCeros, size);
	acumularProb(a, &temp, nCeros, size);

	for (int i = 0; i < size; i++)
		delete[] temp[i];
	delete[] temp;
	return nCeros;
}

// Ordena utilizando el algoritmo Shear short para matrices de tipo DOUBLE
void shearSort(double ***temp, int nCeros, int size) {
	int Log, Rows, Cols;
	int pow = 1, div = 1;

	for (int i = 1; i * i <= (size - nCeros); i++)
		if ((size - nCeros) % i == 0)
			div = i;
	Rows = div;
	Cols = (size - nCeros) / div;
	for (Log = 0; pow <= Rows; Log++)
		pow = pow * 2;

	int h[Rows];
	for (int i = 0; i < Rows; i++)
		h[i] = i * Cols;

	for (int k = 0; k < Log; k++) {
		for (int j = 0; j < Cols / 2; j++) {
			for (int i = 0; i < Rows; i++)
				sortPart1(temp, i * Cols, (i + 1) * Cols, 1, (i % 2 == 0 ? true
						: false));
			for (int i = 0; i < Rows; i++)
				sortPart2(temp, i * Cols, (i + 1) * Cols, 1, (i % 2 == 0 ? true
						: false));
		}
		for (int j = 0; j < Rows / 2; j++) {
			for (int i = 0; i < Cols; i++)
				sortPart1(temp, i, Rows * Cols + i, Cols, true);
			for (int i = 0; i < Cols; i++)
				sortPart2(temp, i, Rows * Cols + i, Cols, true);
		}
	}
	for (int j = 0; j < Cols / 2; j++) {
		for (int i = 0; i < Rows; i++)
			sortPart1(temp, i * Cols, (i + 1) * Cols, 1, true);
		for (int i = 0; i < Rows; i++)
			sortPart2(temp, i * Cols, (i + 1) * Cols, 1, true);
	}
	for (int i = 0; i < Rows; i++)
		h[i] = -1;
}

void sortPart1(double ***temp, int Lo, int Hi, int Nx, bool Up) {
	for (int j = Lo; j + Nx < Hi; j += 2 * Nx)
		if ((Up && (*temp)[j][0] > (*temp)[j + Nx][0]) || (!Up && (*temp)[j][0]
				< (*temp)[j + Nx][0])) {
			double T = (*temp)[j][0];
			double posicion = (*temp)[j][1];
			(*temp)[j][0] = (*temp)[j + Nx][0];
			(*temp)[j][1] = (*temp)[j + Nx][1];
			(*temp)[j + Nx][0] = T;
			(*temp)[j + Nx][1] = posicion;

		}
}

void sortPart2(double ***temp, int Lo, int Hi, int Nx, bool Up) {
	for (int j = Lo + Nx; j + Nx < Hi; j += 2 * Nx)
		if ((Up && (*temp)[j][0] > (*temp)[j + Nx][0]) || (!Up && (*temp)[j][0]
				< (*temp)[j + Nx][0])) {
			double T = (*temp)[j][0];
			double posicion = (*temp)[j][1];
			(*temp)[j][0] = (*temp)[j + Nx][0];
			(*temp)[j][1] = (*temp)[j + Nx][1];
			(*temp)[j + Nx][0] = T;
			(*temp)[j + Nx][1] = posicion;
		}
}

// Acumula las probabibilidades gameticas de la matriz ordenada de manera
// ascendente por sort. Esto permitira buscar luego el genotipo en lugar de
// tener que acumularlo.
void acumularProb(double ***a, double *** temp, int nCeros, int size) {
	for (int i = 0; i < (size - nCeros); i++) {
		(*a)[i][0] = 0;
		(*a)[i][1] = 0;
		if (i == 0)
			(*a)[i][0] = (*temp)[i][0];
		else
			(*a)[i][0] = (*temp)[i][0] + (*a)[i - 1][0];
		(*a)[i][1] = (*temp)[i][1];
	}
}

// Busqueda Binaria en un array double
int binarySearch(double **list, int size, double searchItem, int ultimo) {
	int first = 0;
	int last = size - ultimo - 1;
	int mid;

	bool found = false;
	if (last == 1) {
		mid = (first + last) / 2;
		if (list[mid][0] >= searchItem)
			found = true;
		else
			mid += 1;

	} else {
		// Loop until found or end of list.
		while (first <= last && !found) {
			// Find the middle.
			mid = (first + last) / 2;
			// Compare the middle item to the search item.
			if (mid != 0) {
				if (list[mid][0] >= searchItem && list[mid - 1][0]
						<= searchItem)
					found = true;
				else {
					if (list[mid][0] > searchItem)
						last = mid - 1;
					else
						first = mid + 1;
				}
			} else
				found = true;
		}
	}
	return (int) list[mid][1];
}

// Ordena una matriz [][2] con los individuos a eclosionar de manera creciente
void ordenarEclosionados(long *a, int size) {
	shearSortLong(a, size);
}

// Ordena utilizando el algoritmo Shear short para matrices de tipo DOUBLE
void shearSortLong(long *temp, int size) {
	int Log, Rows, Cols;
	int pow = 1, div = 1;

	for (int i = 1; i * i <= (size); i++)
		if ((size) % i == 0)
			div = i;
	Rows = div;
	Cols = (size) / div;
	for (Log = 0; pow <= Rows; Log++)
		pow = pow * 2;

	int h[Rows];
	for (int i = 0; i < Rows; i++)
		h[i] = i * Cols;

	for (int k = 0; k < Log; k++) {
		for (int j = 0; j < Cols / 2; j++) {
			for (int i = 0; i < Rows; i++)
				sortPart1Long(temp, i * Cols, (i + 1) * Cols, 1,
						(i % 2 == 0 ? true : false));
			for (int i = 0; i < Rows; i++)
				sortPart2Long(temp, i * Cols, (i + 1) * Cols, 1,
						(i % 2 == 0 ? true : false));
		}
		for (int j = 0; j < Rows / 2; j++) {
			for (int i = 0; i < Cols; i++)
				sortPart1Long(temp, i, Rows * Cols + i, Cols, true);
			for (int i = 0; i < Cols; i++)
				sortPart2Long(temp, i, Rows * Cols + i, Cols, true);
		}
	}
	for (int j = 0; j < Cols / 2; j++) {
		for (int i = 0; i < Rows; i++)
			sortPart1Long(temp, i * Cols, (i + 1) * Cols, 1, true);
		for (int i = 0; i < Rows; i++)
			sortPart2Long(temp, i * Cols, (i + 1) * Cols, 1, true);
	}
	for (int i = 0; i < Rows; i++)
		h[i] = -1;
}

void sortPart1Long(long *temp, int Lo, int Hi, int Nx, bool Up) {
	for (int j = Lo; j + Nx < Hi; j += 2 * Nx)
		if ((Up && temp[j] > temp[j + Nx]) || (!Up && temp[j] < temp[j + Nx])) {
			long T = temp[j];
			temp[j] = temp[j + Nx];
			temp[j + Nx] = T;

		}
}

void sortPart2Long(long *temp, int Lo, int Hi, int Nx, bool Up) {
	for (int j = Lo + Nx; j + Nx < Hi; j += 2 * Nx)
		if ((Up && temp[j] > temp[j + Nx]) || (!Up && temp[j] < temp[j + Nx])) {
			long T = temp[j];
			temp[j] = temp[j + Nx];
			temp[j + Nx] = T;
		}
}

