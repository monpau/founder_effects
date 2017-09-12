#ifndef LAGUNA_H_
#define LAGUNA_H_

#include <cmath>
#include <iostream>
#include <vector>
#include <cstring>

#include "OperarGenotipos.h"
#include "Genotipo.h"
#include "Poblacion.h"

using namespace std;

class Laguna {
private:
	string localidad;
	int hidroperiodo;
	long capacidadCarga;
	PoblacionFuente *poblacionFuente;
	Poblacion *pob;

public:
	Laguna(long capacidCarga, int hidroperiod, string local, Poblacion *p, PoblacionFuente *pf);
	Laguna(long capacidCarga, int hidroperiod, string local, PoblacionFuente *p);
	virtual ~Laguna();

	// *** Get and Set
	string getLocalidad();
	void setLocalidad(string localidad);
	int getHidroperiodo();
	void setHidroperiodo(int hidroperiodo);
	long getCapacidadCarga();
	void setCapacidadCarga(long capacidadCarga);
//	int getFundadores();
//	void setFundadores(int fundadores);
//	int getEmigrantes();
//	void setEmigrantes(int emigrantes);
	Poblacion *getPoblacion();
	PoblacionFuente *getPoblacionFuente(); 
	void setPoblacionFuente(PoblacionFuente *poblacion);

	// *** Metodos
	void fundarPoblacionFuente();
	void colonizar (Laguna *lagunaDestino);
	void migrarFuente(Laguna *lagunaDestino);
	void migrarEntre ();
	void eclosionarMigrantes (Laguna *lagunaDestino);
	void establecerTasaCrecimientoGenotipo();
	void crecerAsexual();
};
#endif /*LAGUNA_H_*/
