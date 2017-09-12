#ifndef POBLACION_H_
#define POBLACION_H_

#include <cstring>
#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>

#include "OperarGenotipos.h"
#include "Genotipo.h"
#include "MersenneTwister.h"
#include "BancoHuevos.h"

class PoblacionFuente {
private:
	int numeroLociNeutro;
	int numeroLociSeleccion;
	float mediaAleloNeutro;
	float mediaAleloSeleccion;
	float varianzaAleloNeutro;
	float varianzaAleloSeleccion;
	double *probabilidadesGenotipicas;
	// Sirve para calcular la Fst, lo calculo al fundar la poblacion fuente,
	// [indice genotipo][numero de genes]
	int **nAlelosP;

public:
	// Constructor de la población fuente
	PoblacionFuente();
	PoblacionFuente(int nLociSeleccion, int nLociNeutro,
			float mediaAleloSeleccion, float mediaAleloNeutro,
			float varAleloSeleccion, float varAleloNeutro);
	PoblacionFuente(int nLociSeleccion, int nLociNeutro);
	PoblacionFuente(const PoblacionFuente &p);
	virtual ~PoblacionFuente();

	int getNumeroLociNeutro();
	int getNumeroLociSeleccion();
	float getMediaAleloNeutro();
	float getMediaAleloSeleccion();
	float getVarianzaAleloNeutro();
	float getVarianzaAleloSeleccion();
	int **getNAlelosP();
	void setNAlelosP(int **alelosP);
	double *getProbabilidadesGenotipicas();
	void setProbabilidadesGenotipicas(double *probGenotipicas);
	void print();
};

class Poblacion {
private:
	int numeroGenes;
	long tamanoPoblacional;
	float tasaEclosion;
	float tasaSupervivencia;
	float tasaSanos;
	float tasaCrecimientoPoblacional;
	float componenteEficacia;
	double *recombinacion;
	int *ligamiento;
	int numeroHuevosDiapausicos;
	float tasaMixis;
	float tasaSexAllocation;
	vector<Genotipo> genotipos;
	// Almacena todos los genotipos en forma de String, servirá en la
	// reproduccion sexual.
	string *genotipoString;
	string *genotipoGametoString;
	// Esta matriz sirve para las combinaciones de gemetos
	int **genotipoCombinacionesGametos;
	// Determina si hay o no banco de huevos: true (por defecto) tiene banco de
	// huevos; false no tiene banco antiguo.
	bool existeBanco;
	BancoHuevos *bancoJoven;
	BancoHuevos *bancoViejo;
	// Estas dos matrices sirven para almacenar los datos de eclosionados y
	// supervivientes de los bancos de huevos
	long *auxEclosionados;
	int fundadores;
	int migrantes;
	int *auxMigrados;

public:

	// Constructor para las poblaciones sumidero
	Poblacion();
	Poblacion(int numGenes,
			const float tasaCrecimientoPoblacional,
			const float componenteEficacia, const float tasaEclosion,
			const float tasaSupervivencia, const float tasaSanos,
			double **recombinacion, int** ligamiento,
			const int huevosDiapausicos, const float tasaMixis,
			const float tasaSexAllocation, string **genotipoString,
			string **genotipoGametoString, int ***genotipoCombinacionesGameto,
			int fundadores, int migrantes);
	Poblacion (const Poblacion &p);
	virtual ~Poblacion();

	// *** Get and Set

	long getTamanoPoblacional();
	void setTamanoPoblacional(long tamanoPob);
	//			int getNumeroLociNeutro();
	//			int getNumeroLociSeleccion();
	float getTasaCrecimientoPoblacional();
	vector<Genotipo> *getGenotipos();
	Genotipo *getGenotipo(int id);
	float getComponenteEficacia();
	float getTasaEclosion();
	float getTasaSupervivencia();
	float getTasaSanos();
	void setGenotipos(vector<Genotipo> myv);
	BancoHuevos* getBancoJoven();
	void setBancoJoven(BancoHuevos *bncJoven);
	BancoHuevos* getBancoViejo();
	void setBancoViejo(BancoHuevos *bancoViej);
	string *getGenotipoString();
	void setGenotipoString(string *genotipoStr);
	string *getGenotipoGametoString();
	void setGenotipoGametoString(string *genotipoGametoStr);
	bool isExisteBanco();
	void setExisteBanco(bool existeBnc);
	long *getAuxEclosionados();
	void setAuxEclosionados(long *auxEclosionads);
	int **getGenotipoCombinacionesGametos();
	void setGenotipoCombinacionesGametos(int **genotipoCombinacionesGamets);
	void setTasaSanos(float tSanos);
	int getFundadores();
	void setFundadores(int);
	int getMigrantes();
	void setMigrantes(int);
	int *getAuxMigrados();
	void setAuxMigrados(int *auxMigrad);

	// *** Methods

	// Calcula el tamanyo poblacional total en el momento
	double calcularTamanyoActual();
	// Inicialiar arrayList de genotipos. Todos tienen individuos 0
	void inicializarArrays();
	// Establece la fase gametica del individuo padre
	void establecerFaseGametica(int genotipoPadre, string *resultado);
	// Pool de gametos (piscinita de gametitos)
	void generarPoolGametos(double ***probs);
	// Reproducirse sexualmente
	void reproduccionSexual();
	// Eclosi�n, recolonizacion columna y sedimentacion
	void eclosionEnElBanco();
	void eclosionEnElBancoDeterminista();
	void eclosionEnElBancoAleatorio();
	// Supervivencia
	void sobrevivirEnElBanco();
	void sobrevivirEnElBancoDeterminista();
	void sobrevivirEnElBancoAleatorio();
	void imprimirBanco ();
	void imprimirColumnaAgua ();
};

#endif /*POBLACION_H_*/
