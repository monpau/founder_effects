#include <cmath>
#include <vector>
#include <iostream>

#include "IndicesGeneticos.h"
#include "Laguna.h"
#include "Poblacion.h"
#include "Genotipo.h"

void calcularFstYQst(Laguna *lf, Laguna *l1, Laguna *l2, double **resultado) {
	PoblacionFuente *pf = lf->getPoblacionFuente();
	Poblacion *p1 = l1->getPoblacion();
	Poblacion *p2 = l2->getPoblacion();
	double n1 = p1->calcularTamanyoActual();
	double n2 = p2->calcularTamanyoActual();
	if (n1 == 0 || n2 == 0) {
		(*resultado)[0] = 11111;
		(*resultado)[1] = 11111;
	} else {
		double sumaAlelosSelP1[pf->getNumeroLociSeleccion()];
		double sumaAlelosSelP2[pf->getNumeroLociSeleccion()];
		double sumaAlelosNeutP1[pf->getNumeroLociNeutro()];
		double sumaAlelosNeutP2[pf->getNumeroLociNeutro()];
		for (int j = 0; j < pf->getNumeroLociSeleccion(); j++)
			sumaAlelosSelP1[j] = sumaAlelosSelP2[j] = 0;
		for (int j = 0; j < pf->getNumeroLociNeutro(); j++)
			sumaAlelosNeutP1[j] = sumaAlelosNeutP2[j] = 0;
//		cout << "Genotipos tras eclosi�n (C�lculo de Fst)" << endl;
		for (unsigned int i = 0; i < p1->getGenotipos()->size(); i++) {
			if ((p1->getGenotipo(i)->getIndividuos() != 0) || (p2->getGenotipo(
					i)->getIndividuos() != 0)) {
				for (int j = 0; j < pf->getNumeroLociSeleccion(); j++) {
					sumaAlelosSelP1[j] += pf->getNAlelosP()[i][j]
							* p1->getGenotipo(i)->getIndividuos();
					sumaAlelosSelP2[j] += pf->getNAlelosP()[i][j]
							* p2->getGenotipo(i)->getIndividuos();
				}
				for (int j = 0; j < pf->getNumeroLociNeutro(); j++) {
					sumaAlelosNeutP1[j] += pf->getNAlelosP()[i][j
							+ pf->getNumeroLociSeleccion()]
							* p1->getGenotipo(i)->getIndividuos();
					sumaAlelosNeutP2[j] += pf->getNAlelosP()[i][j
							+ pf->getNumeroLociSeleccion()]
							* p2->getGenotipo(i)->getIndividuos();
				}
//				cout << "Genotipo "<<i<<" "<< convertirIntEnStringDiploide(i,5)<< " "
//						<< p1->getGenotipo(i)->getIndividuos() << " "
//						<< p2->getGenotipo(i)->getIndividuos() << endl;
			}
		}

		if (pf->getNumeroLociSeleccion() > 0) {
			(*resultado)[0] = calcularFst(n1, n2, sumaAlelosSelP1,
					sumaAlelosSelP2, pf->getNumeroLociSeleccion());
			(*resultado)[1] = calcularFst(n1, n2, sumaAlelosNeutP1,
					sumaAlelosNeutP2, pf->getNumeroLociNeutro());
		} else {
			(*resultado)[0] = 99999.0;
			(*resultado)[1] = calcularFst(n1, n2, sumaAlelosNeutP1,
					sumaAlelosNeutP2, pf->getNumeroLociNeutro());
		}
	}
}

double calcularFst(double n1, double n2, double *sumaAlelosP1,
		double *sumaAlelosP2, int nGenes) {

	double nTotAlelos1 = 2 * n1;
	double nTotAlelos2 = 2 * n2;
	double frecP1[nGenes];
	double frecP2[nGenes];
	double frecPTot[nGenes];
	double hS[nGenes];
	double hT[nGenes];
	double fSt[nGenes];
	double mediaHs = 0;
	double mediaHt = 0;
	double mediafSt = 0;
	for (int i = 0; i < nGenes; i++) {
		//		if (n1 > 0)
		frecP1[i] = sumaAlelosP1[i] / nTotAlelos1;
		//		else
		//			frecP1[i] = 0;
		//		if (n2 > 0)
		frecP2[i] = sumaAlelosP2[i] / nTotAlelos2;
		//		else
		//			frecP2[i] = 0;
		frecPTot[i] = (sumaAlelosP1[i] + sumaAlelosP2[i]) / (nTotAlelos1
				+ nTotAlelos2);
		hS[i] = (2 * frecP1[i] * (1 - frecP1[i]) + 2 * frecP2[i] * (1
				- frecP2[i])) / 2;
		hT[i] = 2 * frecPTot[i] * (1 - frecPTot[i]);
		fSt[i] = (hT[i] - hS[i]) / hT[i];
		mediaHs += hS[i];
		mediaHt += hT[i];
//				cout <<frecP1[i] << " (" << sumaAlelosP1[i] << ") " << frecP2[i] << " ("
//						<< sumaAlelosP2[i] << ") "<<endl;
	}
//	for (int r = 0; r < nGenes; r++) {
//		cout<<"Hs "<< hS[r] << " Ht " << hT[r] <<" Fst "<<hS[r]/hT[r]<<endl;
//
//	}

	mediaHs = mediaHs / nGenes;
	mediaHt = mediaHt / nGenes;
	if (mediaHt == 0) {
		mediafSt = 0;

	} else
		mediafSt = (mediaHt - mediaHs) / mediaHt;
//	cout<<"Fst "<<mediafSt<<endl;
	return mediafSt;
}

