#include <cstring>
#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>

#include "Poblacion.h"
#include "OperarGenotipos.h"
#include "Genotipo.h"
#include "MersenneTwister.h"
#include "BancoHuevos.h"

using namespace std;

//Constructor de la poblacion fuente
PoblacionFuente::PoblacionFuente() {
	numeroLociNeutro = 0;
	numeroLociSeleccion = 0;
	mediaAleloNeutro = 0;
	mediaAleloSeleccion = 0;
	varianzaAleloNeutro = 0;
	varianzaAleloSeleccion = 0;
	nAlelosP = NULL;
}

PoblacionFuente::PoblacionFuente(int nLociSeleccion, int nLociNeutro,
		float fAleloSeleccion, float fAleloNeutro, float vAleloSeleccion,
		float vAleloNeutro) {
	numeroLociNeutro = nLociNeutro;
	numeroLociSeleccion = nLociSeleccion;
	mediaAleloNeutro = fAleloNeutro;
	mediaAleloSeleccion = fAleloSeleccion;
	varianzaAleloNeutro = vAleloNeutro;
	varianzaAleloSeleccion = vAleloSeleccion;
	int ids = (int) pow(2.0, (2.0 * (numeroLociNeutro + numeroLociSeleccion)));
	int genes = numeroLociNeutro + numeroLociSeleccion;
	probabilidadesGenotipicas = new double[ids];
	nAlelosP = new int *[ids];
	for (int i = 0; i < ids; i++) {
		nAlelosP[i] = new int[genes];
	}
}

PoblacionFuente::PoblacionFuente(int nLociSeleccion, int nLociNeutro) {
	numeroLociNeutro = nLociNeutro;
	numeroLociSeleccion = nLociSeleccion;
	mediaAleloNeutro = 0;
	mediaAleloSeleccion = 0;
	varianzaAleloNeutro = 0;
	varianzaAleloSeleccion = 0;
	int ids = (int) pow(2.0, (2.0 * (numeroLociNeutro + numeroLociSeleccion)));
	int genes = numeroLociNeutro + numeroLociSeleccion;
	probabilidadesGenotipicas = new double[ids];
	nAlelosP = new int *[ids];
	for (int i = 0; i < ids; i++) {
		nAlelosP[i] = new int[genes];
	}
}

PoblacionFuente::~PoblacionFuente() {
	if (nAlelosP != NULL) {
		for (int i = 0; i < (int) pow(2.0, (2.0 * (numeroLociNeutro
				+ numeroLociSeleccion))); i++) {
			delete[] nAlelosP[i];
		}
		delete[] nAlelosP;
	}
	delete[] probabilidadesGenotipicas;
}

float PoblacionFuente::getMediaAleloNeutro() {
	return mediaAleloNeutro;
}

float PoblacionFuente::getMediaAleloSeleccion() {
	return mediaAleloSeleccion;
}

float PoblacionFuente::getVarianzaAleloNeutro() {
	return varianzaAleloNeutro;
}

float PoblacionFuente::getVarianzaAleloSeleccion() {
	return varianzaAleloSeleccion;
}

int PoblacionFuente::getNumeroLociNeutro() {
	return numeroLociNeutro;
}

int PoblacionFuente::getNumeroLociSeleccion() {
	return numeroLociSeleccion;
}

void PoblacionFuente::print() {
	for (int a = 0; a < (int) pow(2.0, (2.0 * (numeroLociNeutro
			+ numeroLociSeleccion))); a++) {
		for (int b = 0; b < (numeroLociNeutro + numeroLociSeleccion); b++)
			cout << nAlelosP[a][b] << " ";
		cout << probabilidadesGenotipicas[a] << endl;
	}
}

int **PoblacionFuente::getNAlelosP() {
	return nAlelosP;
}

void PoblacionFuente::setNAlelosP(int **alelosP) {
	nAlelosP = alelosP;
}

double *PoblacionFuente::getProbabilidadesGenotipicas() {
	return probabilidadesGenotipicas;
}

void PoblacionFuente::setProbabilidadesGenotipicas(double *probGenotipicas) {
	probabilidadesGenotipicas = probGenotipicas;
}

Poblacion::Poblacion() {

	numeroGenes = 0;
	tamanoPoblacional = 0;
	tasaEclosion = 0;
	tasaSupervivencia = 0;
	tasaSanos = 0;
	tasaCrecimientoPoblacional = 0;
	componenteEficacia = 0;
	recombinacion = NULL;
	ligamiento = NULL;
	numeroHuevosDiapausicos = 0;
	tasaMixis = 0;
	tasaSexAllocation = 0;
	//genotipos = NULL;
	genotipoString = NULL;
	genotipoGametoString = NULL;
	genotipoCombinacionesGametos = NULL;
	existeBanco = false;
	bancoJoven = NULL;
	bancoViejo = NULL;
	auxEclosionados = NULL;

	fundadores = 0;
	migrantes = 0;
	//AuxMigrados almacena los individuos que migran de esta población
	auxMigrados = NULL;
}

//Constructor de la poblacion sumidero
Poblacion::Poblacion(int numGenes, const float tCrecimientoPoblacional,
		float cEficacia, float tEclosion, float tSupervivencia, float tSanos,
		double **recomb, int **ligam, int hDiapausicos, float tMixis,
		float tSexAllocation, string **gString, string **gGametoString,
		int ***gCombinacionesGameto, int funda, int migra) {

	numeroGenes = numGenes;
	tasaCrecimientoPoblacional = tCrecimientoPoblacional;
	tasaSanos = tSanos;
	componenteEficacia = cEficacia;
	tasaEclosion = tEclosion;
	tasaSupervivencia = tSupervivencia;
	recombinacion = (*recomb);
	ligamiento = (*ligam);
	numeroHuevosDiapausicos = hDiapausicos;
	tasaMixis = tMixis;
	existeBanco = true;
	bancoJoven = new BancoHuevos(numGenes);
	bancoViejo = new BancoHuevos(numGenes);
	tasaSexAllocation = tSexAllocation;
	genotipoString = (*gString);
	genotipoGametoString = (*gGametoString);
	genotipoCombinacionesGametos = (*gCombinacionesGameto);

	long tamanyoMax = (long) (pow(2.0, 2.0 * numeroGenes));

	tamanoPoblacional = 0;
	auxEclosionados = new long[tamanyoMax];
	inicializarArrays();

	fundadores = funda;
	migrantes = migra;
	auxMigrados = new int[migra];
	for (int f = 0; f < migrantes; f++) {
		auxMigrados[f] = -1;
	}

}

Poblacion::~Poblacion() {
	delete[] auxMigrados;
	delete[] auxEclosionados;
	if (bancoViejo != NULL)
		delete bancoViejo;
	delete bancoJoven;
}

// *** Get and Set
long Poblacion::getTamanoPoblacional() {
	return tamanoPoblacional;
}

void Poblacion::setTamanoPoblacional(long tamanoPob) {
	tamanoPoblacional = tamanoPob;
}

float Poblacion::getTasaCrecimientoPoblacional() {
	return tasaCrecimientoPoblacional;
}

vector<Genotipo> *Poblacion::getGenotipos() {
	return &genotipos;
}

Genotipo *Poblacion::getGenotipo(int id) {

	return &genotipos.at(id);
}

float Poblacion::getComponenteEficacia() {
	return componenteEficacia;
}

float Poblacion::getTasaEclosion() {
	return tasaEclosion;
}

float Poblacion::getTasaSupervivencia() {
	return tasaSupervivencia;
}

float Poblacion::getTasaSanos() {
	return tasaSanos;
}

void Poblacion::setTasaSanos(float tSanos) {
	tasaSanos = tSanos;
}

void Poblacion::setGenotipos(vector<Genotipo> myv) {
	genotipos = myv;
}

BancoHuevos *Poblacion::getBancoJoven() {
	return bancoJoven;
}

void Poblacion::setBancoJoven(BancoHuevos *bncJoven) {
	bancoJoven = bncJoven;
}

BancoHuevos *Poblacion::getBancoViejo() {
	return bancoViejo;
}

void Poblacion::setBancoViejo(BancoHuevos *bancoViej) {
	bancoViejo = bancoViej;
}

string *Poblacion::getGenotipoString() {
	return genotipoString;
}

void Poblacion::setGenotipoString(string *genotipoStr) {
	genotipoString = genotipoStr;
}

string *Poblacion::getGenotipoGametoString() {
	return genotipoGametoString;
}

void Poblacion::setGenotipoGametoString(string *genotipoGametoStr) {
	genotipoGametoString = genotipoGametoStr;
}

bool Poblacion::isExisteBanco() {
	return existeBanco;
}

void Poblacion::setExisteBanco(bool existeBnc) {
	existeBanco = existeBnc;
}

long *Poblacion::getAuxEclosionados() {
	return auxEclosionados;
}

void Poblacion::setAuxEclosionados(long *auxEclosionads) {
	auxEclosionados = auxEclosionads;
}

int **Poblacion::getGenotipoCombinacionesGametos() {
	return genotipoCombinacionesGametos;
}

void Poblacion::setGenotipoCombinacionesGametos(
		int **genotipoCombinacionesGamets) {
	genotipoCombinacionesGametos = genotipoCombinacionesGamets;
}

int Poblacion::getFundadores() {
	return fundadores;
}

void Poblacion::setFundadores(int funda) {
	fundadores = funda;
}

int Poblacion::getMigrantes() {
	return migrantes;
}

void Poblacion::setMigrantes(int migran) {
	migrantes = migran;
}

int *Poblacion::getAuxMigrados() {
	return auxMigrados;
}

void Poblacion::setAuxMigrados(int *auxMigrad) {
	auxMigrados = auxMigrad;
}

// *** Methods

//*** Calcula el tamanyo poblacional total en el momento
double Poblacion::calcularTamanyoActual() {
	double tamanyoActual = 0;
	for (unsigned int i = 0; i < genotipos.size(); i++) {
		tamanyoActual += genotipos[i].getIndividuos();
	}
	return tamanyoActual;
}

//*** Inicialiar arrayList de genotipos. Todos tienen individuos 0
void Poblacion::inicializarArrays() {
	long tamanyoMax = (long) (pow(2.0, (2.0 * numeroGenes)));
	for (int i = 0; i < tamanyoMax; i++) {
		auxEclosionados[i] = 0;
		genotipos.push_back(Genotipo(i));
	}
}

//*** Establece la fase gametica del individuo progenitor
void Poblacion::establecerFaseGametica(int genotipoPadre, string *resultado) {
	string genotipoPaterno = genotipoString[genotipoPadre];
	resultado[0] = resultado[1] = "";
	for (unsigned int i = 0; i < genotipoPaterno.length();) {
		resultado[0] += genotipoPaterno[i++];
		resultado[1] += genotipoPaterno[i++];
	}
}

//*** Pool de gametos: establecemos la probabilidad de que un gameto haya surgido de ese genotipo
void Poblacion::generarPoolGametos(double ***probfinal) {

	int max = (int) pow(2.0, 1.0 * numeroGenes); //número de combinaciones haploides
	double tamanyoPoblacional = calcularTamanyoActual();
	string faseGametica[2];

	for (unsigned int i = 0; i < genotipos.size(); i++) { //*Escogemos el primer genotipo padre
		if (genotipos[i].getIndividuos() != 0) {
			establecerFaseGametica(genotipos[i].getId(), faseGametica); //* Vemos que fase gametica tiene: faseGametica es una vector con los gametos en formas de 1 y 0
			double individuosGenotipo = genotipos[i].getIndividuos();
			for (int x = 0; x < max; x++) { //* Recorro todos los posibles gametos
				double frecuenciaRecombinacion = 0.5;
				string genotipoGameto = genotipoGametoString[x]; //A que string corresponde ese valor de x
				int idGameto = binstring2int(genotipoGameto);
				for (unsigned int j = 0; j < genotipoGameto.length()
						&& frecuenciaRecombinacion != 0; j++) {
					char c = genotipoGameto[j];

					//* Buscamos si el gameto puede pertener a ese genotipo y a que fase gametica pertence
					int faseGamt = -1;
					if (c == faseGametica[0][j])
						faseGamt = 0;
					else if (c == faseGametica[1][j])
						faseGamt = 1;

					//* Si el gameto no puede pertenecer a ninguna fase gamética (-1), se le asigna prob 0 y se acabó
					if (faseGamt != -1) {
						//* Si el gameto puede surgir de la genotipo padre seguimos
						if (j == 0) {
							if (faseGametica[0][j] == faseGametica[1][j]) //* Primera posicion, si es homozigoto, la probabilidad sera 1, si no la probabilidad por defeto se mantiene 0.5
								frecuenciaRecombinacion = 1;
						} else {
							if (faseGametica[0][j] == faseGametica[1][j]) //* Para el resto de posiciones, si es homozigoto, probabilidad 1
								frecuenciaRecombinacion *= 1;
							else {
								if (ligamiento[j] == -1) //* Si es heterozigoto y el gen no esta ligado (dependencia == -1) la probabilidad es 0.5
									frecuenciaRecombinacion *= 0.5;
								else {
									if (faseGametica[0][ligamiento[j]]
											== faseGametica[1][ligamiento[j]]) { //* Si esta ligado a un locus homozigoto, la probabilidad es del 0,5
										frecuenciaRecombinacion *= 0.5;
									} else {
										char cAnt =
												genotipoGameto[ligamiento[j]];
										int faseGamtAnt = -1;
										// Comprobamos la fase del gen al que esta ligado
										if (cAnt
												== faseGametica[0][ligamiento[j]])
											faseGamtAnt = 0;
										else
											faseGamtAnt = 1;
										// Si ambos loci estan en fase la probabilidad es la complementaria, si no es la indicada
										if (faseGametica[faseGamtAnt]
												== faseGametica[faseGamt])
											frecuenciaRecombinacion
													*= (1
															- recombinacion[j]);
										else
											frecuenciaRecombinacion *= recombinacion[j];
									}
								}
							}
						}
					} else
						frecuenciaRecombinacion = 0;
				}
				(*probfinal)[idGameto][0] += frecuenciaRecombinacion
						* (individuosGenotipo / tamanyoPoblacional);
				(*probfinal)[idGameto][1] = idGameto;
			}
		}
	}
//	for(int f = 0;f<max; f++){
//		cout<<  convertirIntEnStringHaploide(f,3)<< " "<<(*probfinal)[f][0]<<" "<<(*probfinal)[f][1]<<endl;
//	}
//	cout<<endl;
	//return probfinal;
}

//***Reproducirse sexualmente
void Poblacion::reproduccionSexual() {
	int size = (int) pow(2.0, numeroGenes);

	// Vector para recoger las probabilidades finales de todos los gametos
	double **probGamet = new double *[size];
	for (int i = 0; i < size; i++) {
		probGamet[i] = new double[2];
		probGamet[i][0] = 0;
		probGamet[i][1] = 0;
	}

	generarPoolGametos(&probGamet);
	int nCeros = ordenar(&probGamet, size);

	// Elimino ya aquellos que moriran en el banco de huevos ese ciclo
	long huevosTotal = (long) (calcularTamanyoActual()
			* numeroHuevosDiapausicos * tasaMixis * tasaSexAllocation
			* tasaSanos * tasaSupervivencia);

	MTRand::uint32 seed[MTRand::N];
	for (int s = 0; s < MTRand::N; ++s) {
		seed[s] = time(NULL) * s;
	}
	MTRand rnd(seed);

	for (int i = 0; i < huevosTotal; i++) {
		int genotipoHuevoInt = 0;
		double d1 = rnd.randDblExc();
		double d2 = rnd.randDblExc();
		int a = binarySearch(probGamet, size, d1, nCeros);
		int b = binarySearch(probGamet, size, d2, nCeros);
		genotipoHuevoInt = genotipoCombinacionesGametos[a][b];
		(bancoJoven->getNHuevos())[genotipoHuevoInt] += 1;
	}

	for (int i = 0; i < size; i++) {
		delete[] probGamet[i];
	}
	delete[] probGamet;
}

//***
void Poblacion::eclosionEnElBanco() {
	int noCeros = 0;
	long numHuev = 0;
	for (int i = 0; i < bancoJoven->size(); i++) {
		if (bancoJoven->getNHuevos()[i] != 0) {
			numHuev += bancoJoven->getNHuevos()[i];
			noCeros++;
		}
	}
	if(noCeros>0){
	int aux = (int) numHuev / noCeros * tasaEclosion;
	if (aux >= 1) {
		eclosionEnElBancoDeterminista();
	} else {
		eclosionEnElBancoAleatorio();
	}
	}else{
		cerr<<"La poblacion no eclosiona porque esta extinguida"<<endl;
	}

}

//*** Eclosion, recolonizacion columna y sedimentacion VERSION DETERMINISTA
void Poblacion::eclosionEnElBancoDeterminista() {
	long *bhj = bancoJoven->getNHuevos();
	if (existeBanco == true) {
		long *bhv = bancoViejo->getNHuevos();

		long eclosionadosJoven = 0;
		long eclosionadosViejo = 0;

		for (int j = 0; j < bancoViejo->size(); j++) {
			//De cada genotipo del banco joven eclosiono los que corresponden
			eclosionadosJoven = (long) (bhj[j] * tasaEclosion);
			//Resto del banco joven aquellos que han eclosionado
			bhj[j] -= eclosionadosJoven;
			//De cada genotipo del banco viejo eclosiono los que corresponden
			eclosionadosViejo = (long) (bhv[j] * tasaEclosion);
			//Resto los eclosionados del banco viejo
			bhv[j] -= eclosionadosViejo;
			//Pongo en la columna de agua a los eclosionados de los dos bancos
			genotipos[j].setIndividuos(eclosionadosJoven + eclosionadosViejo);
			//Pongo a los huevos que quedan del banco joven en el banco viejo
			bhv[j] += bhj[j];
			//Vacio el banco joven
			bhj[j] = 0;
		}

	} else {
		long eclosionadosJoven = 0;
		//		cout << "Eclosionan Deterministicamente ";
		for (int j = 0; j < bancoJoven->size(); j++) {
			eclosionadosJoven = (long) (bhj[j] * tasaEclosion);
			bhj[j] -= eclosionadosJoven;
			genotipos[j].setIndividuos(eclosionadosJoven);
			//			if (eclosionadosJoven > 0)
			//				cout << j << " (" << eclosionadosJoven << ") ";
			bhj[j] = 0;
		}
		//		cout << endl;
	}
}

// Eclosion, recolonizacion columna y sedimentacion VERSION ALEATORIA

void Poblacion::eclosionEnElBancoAleatorio() {
	long *bhj = bancoJoven->getNHuevos();
	long *bhv;
	int numGentp = bancoJoven->size();
	int numHuevJ = 0;
	int numHuevV;
	int noCerosJ = 0;
	int noCerosV;
	long **auxEclosionJ = new long *[numGentp];
	long **auxEclosionV = new long *[numGentp];

	for (int f = 0; f < numGentp; f++) {
		auxEclosionJ[f] = new long[2];
		auxEclosionV[f] = new long[2];
		genotipos[f].setIndividuos(0);
		if (bhj[f] != 0) {
			auxEclosionJ[noCerosJ][0] = bhj[f];
			auxEclosionJ[noCerosJ][1] = f;
			numHuevJ += bhj[f];
			noCerosJ++;
		}
	}

	if (existeBanco == true) {
		bhv = bancoViejo->getNHuevos();
		numHuevV = 0;
		noCerosV = 0;
		long *bhv = bancoViejo->getNHuevos();
		for (int i = 0; i < numGentp; i++) {
			if (bhv[i] != 0) {
				auxEclosionV[noCerosV][0] = bhv[i];
				auxEclosionV[noCerosV][1] = i;
				numHuevV += bhv[i];
				noCerosV++;
			}
		}
	}
	int numEclosJ = numHuevJ * tasaEclosion;
	int numEclosV = numHuevV * tasaEclosion;
	MTRand::uint32 seed[MTRand::N];
	for (int s = 0; s < MTRand::N; ++s) {
		seed[s] = time(NULL) * s;
	}
	MTRand rnd(seed);
	bool encontrado = false;
	int numEncontrados = 0;
	int acum = 0;
	int rand = 0;
	int j = 0;
	long aux1 = 0;
	//	cout << "Eclosionan Aleatoriamente " << endl;
	for (int x = 0; x < numEclosJ; x++) {
		rand = rnd.randInt(numHuevJ - numEncontrados - 1) + 1;
		j = 0;
		acum = 0;
		encontrado = false;
		//		cout << "Bicho n�mero " << rand << " de " << numHuevJ - numEncontrados;
		do {
			acum += auxEclosionJ[j][0];
			if (rand <= acum) {
				encontrado = true;
				aux1 = genotipos[(auxEclosionJ[j][1])].getIndividuos();
				genotipos[(auxEclosionJ[j][1])].setIndividuos(aux1 + 1);
				//				cout << " genotipo" << auxEclosionJ[j][1];
				auxEclosionJ[j][0]--;
				bhj[j]--;
				numEncontrados++;
			} else
				j++;
		} while (encontrado == false);
		//		cout << endl;
	}

	if (existeBanco == true) {
		bhv = bancoViejo->getNHuevos();
		numEncontrados = 0;
		for (int x = 0; x < numEclosV; x++) {
			rand = rnd.randInt(numHuevV - numEncontrados - 1) + 1;
			j = 0;
			acum = 0;
			encontrado = false;
			do {
				acum += auxEclosionV[j][0];
				if (rand <= acum) {
					encontrado = true;
					aux1 = genotipos[(auxEclosionV[j][1])].getIndividuos();
					genotipos[(auxEclosionV[j][1])].setIndividuos(aux1 + 1);
					auxEclosionV[j][0]--;
					bhv[j]--;
					numEncontrados++;
				} else
					j++;
			} while (encontrado == false);
		}
	}

	//*** Los huevos del bancoJoven sedimentan en el bancoViejo
	if (existeBanco == true) {
		for (int j = 0; j < bancoViejo->size(); j++) {
			bhv[j] += bhj[j];
			bhj[j] = 0;
		}
	} else {
		for (int j = 0; j < bancoJoven->size(); j++) {
			bhj[j] = 0;
		}
	}

	for (int f = 0; f < numGentp; f++) {
		delete[] auxEclosionJ[f];
		delete[] auxEclosionV[f];
	}
	delete[] auxEclosionJ;
	delete[] auxEclosionV;
}

void Poblacion::sobrevivirEnElBanco() {
	if (existeBanco == true) {
		int noCeros = 0;
		long numHuev = 0;
		for (int i = 0; i < bancoJoven->size(); i++) {
			if (bancoJoven->getNHuevos()[i] != 0) {
				numHuev += bancoJoven->getNHuevos()[i];
				noCeros++;
			}
		}
		if (noCeros > 0) {
			int aux = (int) numHuev / noCeros * tasaSupervivencia;
			if (aux >= 1) {
				sobrevivirEnElBancoDeterminista();
			} else {
				sobrevivirEnElBancoAleatorio();
			}
		}else{
			cerr<<"La poblacion se ha extinguido "<<endl;
		}
	} else {
		setBancoViejo(NULL);
	}
}

//Supervivencia en el banco de huevos DETERMINISTA
void Poblacion::sobrevivirEnElBancoDeterminista() {
	long *bhv = bancoViejo->getNHuevos();
	for (int j = 0; j < bancoViejo->size(); j++) {
		//Dejo en el banco viejo a aquellos que sobreviven
		bhv[j] = (long) (bhv[j] * tasaSupervivencia);
	}
}

//Supervivencia en el banco de huevos Aleatorio
void Poblacion::sobrevivirEnElBancoAleatorio() {
	long *bhv = bancoViejo->getNHuevos();
	int numGentp = bancoJoven->size();
	int numHuev = 0;
	int noCeros = 0;

	long **auxSupervivencia = new long *[numGentp];
	for (int f = 0; f < numGentp; f++) {
		auxSupervivencia[f] = new long[2];
		if (bhv[f] != 0) {
			auxSupervivencia[noCeros][0] = bhv[f];
			auxSupervivencia[noCeros][1] = f;
			numHuev += bhv[f];
			noCeros++;
		}
	}

	int numMuertos = numHuev * tasaSupervivencia;

	MTRand::uint32 seed[MTRand::N];
	for (int s = 0; s < MTRand::N; ++s) {
		seed[s] = time(NULL) * s;
	}
	MTRand rnd(seed);

	bool encontrado = false;
	int numEncontrados = 0;
	int acum = 0;
	int rand = 0;
	int j = 0;
	for (int x = 0; x < numMuertos; x++) {
		rand = rnd.randInt(numHuev - numEncontrados - 1) + 1;
		j = 0;
		acum = 0;
		encontrado = false;
		do {
			acum += auxSupervivencia[j][0];
			if (rand <= acum) {
				encontrado = true;
				auxSupervivencia[j][0]--;
				bhv[j]--;
				numEncontrados++;
			} else
				j++;
		} while (encontrado == false);
	}

	for (int f = 0; f < numGentp; f++) {
		delete[] auxSupervivencia[f];
	}
	delete[] auxSupervivencia;
}

void Poblacion::imprimirBanco() {
	long *bj = bancoJoven->getNHuevos();
	if (existeBanco == true) {
		long *bv = bancoJoven->getNHuevos();
		for (int i = 0; i < bancoJoven->size(); i++) {
			if (bj[i] != 0 || bv[i] != 0) {
				cout << "Genotipo " << i << " " << bj[i] << " + " << bv[i]
						<< endl;
			}
		}
	} else {
		for (int i = 0; i < bancoJoven->size(); i++) {
			if (bj[i] != 0) {
				cout << "Genotipo " << i << " " << bj[i] << endl;
			}
		}

	}
}

void Poblacion::imprimirColumnaAgua() {
	for (int i = 0; i < bancoJoven->size(); i++) {
		if (genotipos.at(i).getIndividuos() != 0)
			cout << "(" << convertirIntEnStringDiploide(i,10) << ") " << genotipos.at(i).getIndividuos() << " ";
	}
	cout << endl;
}

