#include <vector>
#include <cmath>
#include <iostream>

#include "BancoHuevos.h"

using namespace std;

BancoHuevos::BancoHuevos(int tamanyoHaplotipo)
{
	nHuevosSize = (int)(pow(2, 2*tamanyoHaplotipo));
	nHuevos= new long[nHuevosSize];
	
	for(int i=0; i < nHuevosSize; i++) nHuevos[i]=0;
}

BancoHuevos::~BancoHuevos()
{
	delete [] nHuevos;
}

long BancoHuevos::calcularTamanyoBanco() {
		long tamanyoActual = 0;
		for (int i = 0; i < nHuevosSize; i++) {
			tamanyoActual += nHuevos[i];
		}
		return tamanyoActual;
}

void BancoHuevos::setNHuevos(long *huevos) {
	nHuevos=huevos;
}

long *BancoHuevos::getNHuevos() {
	return nHuevos;
}

int BancoHuevos::size() {
	return nHuevosSize;
}


