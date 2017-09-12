#include <cmath>
#include <iostream>
#include <vector>
#include <cstring>

#include "Laguna.h"
#include "OperarGenotipos.h"
#include "Genotipo.h"
#include "Poblacion.h"

using namespace std;

Laguna::Laguna(long capacidCarga, int hidroperiod, string local, Poblacion *p,
		PoblacionFuente *pf) {
	localidad = local;
	capacidadCarga = capacidCarga;
	hidroperiodo = hidroperiod;
	pob = p;
	poblacionFuente = pf;
}

Laguna::Laguna(long capacidCarga, int hidroperiod, string local,
		PoblacionFuente *p) {
	localidad = local;
	capacidadCarga = capacidCarga;
	hidroperiodo = hidroperiod;
	poblacionFuente = p; //PoblacionFuente (p);
}

Laguna::~Laguna() {
}

// Get and Set

string Laguna::getLocalidad() {
	return localidad;
}

void Laguna::setLocalidad(string local) {
	localidad = local;
}

int Laguna::getHidroperiodo() {
	return hidroperiodo;
}

void Laguna::setHidroperiodo(int hidroperiod) {
	hidroperiodo = hidroperiod;
}

long Laguna::getCapacidadCarga() {
	return capacidadCarga;
}

void Laguna::setCapacidadCarga(long capacidCarga) {
	capacidadCarga = capacidCarga;
}

//int Laguna::getFundadores() {
//	return fundadores;
//}
//
//void Laguna::setFundadores(int fundadores) {
//	fundadores = fundadores;
//}
//
//int Laguna::getEmigrantes() {
//	return emigrantes;
//}
//
//void Laguna::setEmigrantes(int emigrantes) {
//	emigrantes = emigrantes;
//}

Poblacion *Laguna::getPoblacion() {
	return pob;
}

PoblacionFuente *Laguna::getPoblacionFuente() {
	return poblacionFuente;
}

void Laguna::setPoblacionFuente(PoblacionFuente *poblac) {
	poblacionFuente = poblac;
}

// Metodos

void Laguna::fundarPoblacionFuente() {

	int genesSeleccion = poblacionFuente->getNumeroLociSeleccion();
	float mediaSeleccion = poblacionFuente->getMediaAleloSeleccion();
	float destSelec = poblacionFuente->getVarianzaAleloSeleccion();
	int genesNeutros = poblacionFuente->getNumeroLociNeutro();
	float mediaNeutro = poblacionFuente->getMediaAleloNeutro();
	float destNeutro = poblacionFuente->getVarianzaAleloNeutro();

	float frecLoci[genesSeleccion + genesNeutros];
	float frecLociSelec[genesSeleccion];
	float frecLociNeutr[genesNeutros];

	//*** Establece la frecuencia de los alelos 0 en la poblacion fuente

	if (mediaNeutro == 0 && mediaSeleccion == 0) {
		if (genesSeleccion == 0) {
			frecuenciasLociUniforme(genesNeutros, frecLoci);

		} else if (genesNeutros == 0) {
			frecuenciasLociUniforme(genesSeleccion, frecLoci);

		} else {
			frecuenciasLociUniforme(genesSeleccion, frecLociSelec);
			frecuenciasLociUniforme(genesNeutros, frecLociNeutr);
			//frecLoci = new float[genesNeutros+genesSeleccion];
			for (int i = 0; i < genesSeleccion; i++)
				frecLoci[i] = frecLociSelec[i];
			for (int i = 0; i < genesNeutros; i++)
				frecLoci[i + genesSeleccion] = frecLociNeutr[i];
		}

	} else {
		if (genesSeleccion == 0) {
			frecuenciasLoci(mediaNeutro, destNeutro, genesNeutros, frecLoci);

		} else if (genesNeutros == 0) {
			frecuenciasLoci(mediaSeleccion, destSelec, genesSeleccion, frecLoci);

		} else {
			frecuenciasLoci(mediaSeleccion, destSelec, genesSeleccion,
					frecLociSelec);
			frecuenciasLoci(mediaNeutro, destNeutro, genesNeutros,
					frecLociNeutr);
			//frecLoci = new float[genesNeutros+genesSeleccion];
			for (int i = 0; i < genesSeleccion; i++)
				frecLoci[i] = frecLociSelec[i];
			for (int i = 0; i < genesNeutros; i++)
				frecLoci[i + genesSeleccion] = frecLociNeutr[i];
		}
	}

	// Establecer el tamanyo maximo de la matriz de genotipos
	int genes = genesNeutros + genesSeleccion;
	int totalGenotipos = (int) pow(2.0, 2.0 * genes);
	bool *genotipoBooleano = new bool[2 * genes];

	// *** Establece la probabilidad de cada uno de los genotipos y rellena NAlelosP con 2 si es homozigoto
	//para el alelo 0; 1 si es heterozigoto y 0 si es homozigoto para el alelo 1
	for (int i = 0; i < totalGenotipos; i++) {
		double prob = 1;
		convertirIntEnBool(i, 2 * genes, genotipoBooleano);

		for (int j = 0; j < genes; j++) {
			if (genotipoBooleano[j * 2] == false && genotipoBooleano[j * 2 + 1]
					== false) {
				prob *= (frecLoci[j] * frecLoci[j]);
				poblacionFuente->getNAlelosP()[i][j] = 2;
			} else if ((genotipoBooleano[j * 2] == false && genotipoBooleano[j
					* 2 + 1] == true) ^ (genotipoBooleano[j * 2] == true
					&& genotipoBooleano[j * 2 + 1] == false)) {
				prob *= (frecLoci[j] * (1 - frecLoci[j]));
				poblacionFuente->getNAlelosP()[i][j] = 1;
			} else {
				prob *= ((1 - frecLoci[j]) * (1 - frecLoci[j]));
				poblacionFuente->getNAlelosP()[i][j] = 0;
			}
		}
		//		for (int g = 0; g<2*genes;g++)
		//			cout<<genotipoBooleano[g];
		//		cout<<" "<<prob<<endl;
		poblacionFuente->getProbabilidadesGenotipicas()[i] = prob;
	}
	delete[] genotipoBooleano;
}

//*** Funcion que establece los individuos que llegan primero desde la laguna fuente
void Laguna::colonizar(Laguna *lagunaDestino) {
	int i, j;
	MTRand rnd;
	int fundadores = lagunaDestino->getPoblacion()->getFundadores();
	double rand, probAcumulada;
	for (i = 0; i < fundadores; i++) {
		rand = rnd.randExc();
		probAcumulada = 0;
		for (j = 0; (rand > probAcumulada); j++) {
			probAcumulada += poblacionFuente->getProbabilidadesGenotipicas()[j];
		}
		j--;
//		cout << "Genotipo fundador " << convertirIntEnStringDiploide(j,3)<< endl;
		int nuevotam =
				lagunaDestino->getPoblacion()->getGenotipo(j)->getIndividuos()
						+ 1;
		lagunaDestino->getPoblacion()->getGenotipo(j)->setIndividuos(nuevotam);
	}
}

//** Función que establece los individuos aleatorios que migran desde la laguna fuente
void Laguna::migrarFuente(Laguna *lagunaDestino) {
	int i, j;
	MTRand rnd;
	int emigrantes = lagunaDestino->getPoblacion()->getMigrantes();
	double rand, probAcumulada;
	for (i = 0; i < emigrantes; i++) {
		rand = rnd.randDblExc();
		probAcumulada = 0;
		for (j = 0; (rand > probAcumulada); j++) {
			probAcumulada += poblacionFuente->getProbabilidadesGenotipicas()[j];
		}
		j--;
		int nuevotam =
				lagunaDestino->getPoblacion()->getGenotipo(j)->getIndividuos()
						+ 1;

		lagunaDestino->getPoblacion()->getGenotipo(j)->setIndividuos(nuevotam);
	}

}

//*** Establece los individuos que migran de una laguna a la otra.
//***   ANTES de la migraci�n los huevos del banco de huevos viejo mueren

void Laguna::migrarEntre() {
	MTRand::uint32 seed[MTRand::N];
	for (int s = 0; s < MTRand::N; ++s) {
		seed[s] = time(NULL) * s;
	}
	MTRand rnd(seed);

	int i, j;
	bool encontrado = false;
	long acumulado = 0;
	long rand;
	long huevAcumuladosJ = 0;

	int emigrantes = getPoblacion()->getMigrantes();
	int excesoMigrantes = 0;
	long huevosJovenes =
			getPoblacion()->getBancoJoven()->calcularTamanyoBanco();
	long huevosTotales = huevosJovenes;
	if (getPoblacion()->isExisteBanco() == true)
		huevosTotales
				+= getPoblacion()->getBancoViejo()->calcularTamanyoBanco();

	// Protecci�n frente a la posibilidad de que haya m�s migrantes que huevos en el banco
	if (huevosTotales < emigrantes) {
		excesoMigrantes = emigrantes - huevosTotales;
		emigrantes = huevosTotales;
	}
	if (huevosTotales != 0) {
		for (i = 0; i < emigrantes; i++) {
			rand = rnd.randInt(huevosTotales - 1) + 1;
			encontrado = false;
			j = 0;
			acumulado = 0;
			if (rand <= (huevosJovenes - huevAcumuladosJ)) {
				do {
					acumulado
							+= getPoblacion()->getBancoJoven()->getNHuevos()[j];
					if (rand <= acumulado) {
						encontrado = true;
						huevAcumuladosJ++;
						huevosTotales--;
						getPoblacion()->getAuxMigrados()[i] = j;
						getPoblacion()->getBancoJoven()->getNHuevos()[j]--;
						//						cout << "migra " << j << " y quedan " << huevosTotales
						//								<< " huevos" << endl;
					} else
						j++;
				} while (encontrado == false);
			} else if (rand > (huevosJovenes - huevAcumuladosJ)
					&& getPoblacion()->isExisteBanco() == true) {
				acumulado = huevosJovenes - huevAcumuladosJ;
				do {
					acumulado
							+= getPoblacion()->getBancoViejo()->getNHuevos()[j];
					if (rand <= acumulado) {
						encontrado = true;
						huevosTotales--;
						getPoblacion()->getAuxMigrados()[i] = j;
						getPoblacion()->getBancoViejo()->getNHuevos()[j]--;
					} else
						j++;
				} while (encontrado == false);
			} else
				cerr << "no he podido migrar " << endl;
		}
	}
	if (excesoMigrantes != 0) {
		for (int r = emigrantes; r < (excesoMigrantes + emigrantes); r++) {
			getPoblacion()->getAuxMigrados()[r] = -1;
		}
	}
}

//** Eclosiona los individuos migrantes
void Laguna::eclosionarMigrantes(Laguna *lagunaDestino) {
	int migrantes = getPoblacion()->getMigrantes();
	int aux = 0;
	double aux2;
	Poblacion *p = lagunaDestino->getPoblacion();
	for (int i = 0; i < migrantes; i++) {
		if (getPoblacion()->getAuxMigrados()[i] != -1) {
			aux = getPoblacion()->getAuxMigrados()[i];
			aux2 = p->getGenotipo(aux)->getIndividuos();
			p->getGenotipo(aux)->setIndividuos(aux2 + 1);
			//			cout << "A�ado el bicho " << aux << endl;

		}
		getPoblacion()->getAuxMigrados()[i] = -1;
	}
}

void Laguna::establecerTasaCrecimientoGenotipo() {
	//	 Alelo 0 domina en loc 1 y el alelo 1 domina en loc 2
	int genSelecc = poblacionFuente->getNumeroLociSeleccion();

	for (unsigned int i = 0; i < pob->getGenotipos()->size(); i++) {
		//Genotipo geno = genos[i];
		if (genSelecc == 0) {
			pob->getGenotipo(i)->setTasaCrecimientoGenotipo(
					pob->getTasaCrecimientoPoblacional());
		} else {
			int nHomozigotosAlelo0 = 0;
			int nHomozigotosAlelo1 = 0;
			for (int j = 0; j < genSelecc; j++) {
				if (poblacionFuente->getNAlelosP()[i][j] == 2)
					nHomozigotosAlelo0 += 1;
				else if (poblacionFuente->getNAlelosP()[i][j] == 0)
					nHomozigotosAlelo1 += 1;
			}
			if (localidad == "localidad1") {
				if (pob->getTasaCrecimientoPoblacional() + ((nHomozigotosAlelo0
						- nHomozigotosAlelo1) * pob->getComponenteEficacia())
						< 0)
					pob->getGenotipo(i)->setTasaCrecimientoGenotipo(0);
				else
					pob->getGenotipo(i)->setTasaCrecimientoGenotipo(
							pob->getTasaCrecimientoPoblacional()
									+ ((nHomozigotosAlelo0 - nHomozigotosAlelo1)
											* pob->getComponenteEficacia()));
			} else if (pob->getTasaCrecimientoPoblacional()
					+ ((nHomozigotosAlelo1 - nHomozigotosAlelo0)
							* pob->getComponenteEficacia()) < 0) {
				pob->getGenotipo(i)->setTasaCrecimientoGenotipo(0);
			} else
				pob->getGenotipo(i)->setTasaCrecimientoGenotipo(
						pob->getTasaCrecimientoPoblacional()
								+ ((nHomozigotosAlelo1 - nHomozigotosAlelo0)
										* pob->getComponenteEficacia()));
		}
	}
}

void Laguna::crecerAsexual() {
	double tamanyoPoblParcial = 0;
	int noCeros = 0;
	int genotipos = getPoblacion()->getGenotipos()->size();
	double **arrayTemp = new double *[genotipos];
	for (int f = 0; f < genotipos; f++)
		arrayTemp[f] = new double[2];
	for (int i = 0; i < genotipos; i++) {
		if (getPoblacion()->getGenotipo(i)->getIndividuos() != 0) {
			arrayTemp[noCeros][0]
					= getPoblacion()->getGenotipo(i)->getIndividuos();
			arrayTemp[noCeros][1] = i;
			tamanyoPoblParcial
					+= getPoblacion()->getGenotipo(i)->getIndividuos();
			noCeros++;
		}

	}

	double contador = 0.0;
	float rg;
	double tamanyoInicial, tamanyoFinal;
	//	cout << "Individuos tras crec asex " << endl;
	for (int j = 0; j < hidroperiodo; j++) {
		if (j != 0)
			tamanyoPoblParcial = contador;
		contador = 0.0;
		for (int i = 0; i < noCeros; i++) {
			rg
					= getPoblacion()->getGenotipo(arrayTemp[i][1])->getTasaCrecimientoGenotipo();
			tamanyoInicial = arrayTemp[i][0];
			tamanyoFinal = (tamanyoInicial * (1 + rg * (1 - (tamanyoPoblParcial
					/ getCapacidadCarga()))));
			if (tamanyoFinal <= 0) {
				arrayTemp[i][0] = 0.0;
				contador += 0;
			} else {
				arrayTemp[i][0] = tamanyoFinal;
				contador += tamanyoFinal;
			}
		}

	}
	for (int i = 0; i < noCeros; i++) {
		long redondeado = arrayTemp[i][0];
		//		cout<<arrayTemp[i][0]<<" "<<redondeado;
		getPoblacion()->getGenotipo(arrayTemp[i][1])->setIndividuos(redondeado);
		//				if(getPoblacion()->getGenotipo(arrayTemp[i][1])->getIndividuos()!=0)
		//				cout <<"Genotipo "<<i<< " ("
		//						<< getPoblacion()->getGenotipo(arrayTemp[i][1])->getIndividuos() << ") "<<endl;
	}
	for (int f = 0; f < genotipos; f++)
		delete[] arrayTemp[f];
	delete[] arrayTemp;
}
