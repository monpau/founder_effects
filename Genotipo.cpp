#include "Genotipo.h"

Genotipo::Genotipo(int genotipo)
{
	id=genotipo;
	individuos= 0.0;
	tasaCrecimientoGenotipo= 0.0;
}

Genotipo::~Genotipo()
{
}

double Genotipo::getIndividuos() {
	return individuos;
}

void Genotipo::setIndividuos(double nind) {
	individuos = nind;
}

float Genotipo::getTasaCrecimientoGenotipo() {
	return tasaCrecimientoGenotipo;
}

int Genotipo::getId() {
	return id;
}

void Genotipo::Genotipo::setId(int newid) {
	id = newid;
}

void Genotipo::setTasaCrecimientoGenotipo(float tasacrec) {
	tasaCrecimientoGenotipo = tasacrec;
}

