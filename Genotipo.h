#ifndef GENOTIPO_H_
#define GENOTIPO_H_

class Genotipo
{
private:
	
	double individuos;
	int id;
	float tasaCrecimientoGenotipo;

public:
	Genotipo(int genotipo);
	virtual ~Genotipo();

	double getIndividuos();
	void setIndividuos(double nind);
	float getTasaCrecimientoGenotipo();
	int getId();
	void setId(int newid);
	void setTasaCrecimientoGenotipo(float tasacrec);
};

#endif /*GENOTIPO_H_*/
