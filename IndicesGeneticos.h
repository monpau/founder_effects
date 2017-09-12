#ifndef INDICESGENETICOS_H_
#define INDICESGENETICOS_H_

#include <cmath>
#include <vector>
#include <iostream>

#include "Laguna.h"
#include "Poblacion.h"
#include "Genotipo.h"


void calcularFstYQst(Laguna *lf, Laguna *l1, Laguna *l2, double **resultado);
double calcularFst(double n1, double n2, double *sumaAlelosP1, double *sumaAlelosP2, int nGenes); 

#endif /*INDICESGENETICOS_H_*/
