#ifndef BANCOHUEVOS_H_
#define BANCOHUEVOS_H_

#include <vector>
#include <cmath>
#include <iostream>

#include "BancoHuevos.h"

using namespace std;

class BancoHuevos {
	private:
		long *nHuevos;
		int nHuevosSize;
	public:
		BancoHuevos(int tamanyoHaplotipo);
		virtual ~BancoHuevos();
	
		long calcularTamanyoBanco();
		void setNHuevos(long *huevos);
		long *getNHuevos();
		int size();

};

#endif /*BANCOHUEVOS_H_*/
