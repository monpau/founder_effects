/* Author: Javier Montero-Pau
 *
 * Code to perform simulations to explore the genetic impact of founder effects
 * and local adaptation in zooplanktonic populations.
 * Details of the implementation are described in
 *
 * The program requires the following input parameters:
 * [1] Diapausing egg bank (boolean, default TRUE)
 * [2] population growth (r) (float)
 * [3] carrying capacity (float)
 * [4] number of migrants (integer)
 * [5] number of founder (integer)
 * [6] number of genes under selection (integer)
 * [7] number of neutral genes (integer)
 * [8] fitness component (float)
 * [9] genetic recombination rate (float)
 * [10] output file"
 *
 * Other parameters such as length of lake hydroperiod, production, survival and
 * deterioration dynamics of diapausing eggs can also be modified
 */

#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>

#include "Poblacion.h"
#include "OperarGenotipos.h"
#include "IndicesGeneticos.h"
#include "Laguna.h"
#include "BancoHuevos.h"

using namespace std;

int main(int argc, const char* argv[]) {

	if (argc != 11)
		cerr
				<< "[1] Diapausing egg bank (boolean, default TRUE), [2] population growth (r) (float), [3] carrying capacity (float), [4] number of migrants (integer), [5] number of founder (integer), [6] number of genes under selection (integer), [7] number of neutral genes (integer), [8] fitness component (float), [9] genetic recombination rate (float), [10] output file"
				<< endl;

	char reo[10];
	strcpy(reo, argv[10]);
	ofstream myfile;
	myfile.open(reo);

	bool existeBanco = atoi(argv[1]);

	double *resultadoFst = new double[2];
	long capacidadCarga = atol(argv[3]);
	int migrantes = atoi(argv[4]);
	int fundadores = atoi(argv[5]);
	int genSelc = atoi(argv[6]);
	int genNeut = atoi(argv[7]);

		bool seleccion = false;
	if (genSelc > 0)
		seleccion = true;

	float tasaCrecPob = atof(argv[2]);
	float componenteEfic = atof(argv[8]);

	double *recombinacion;
	int *ligamiento;
	float tasaRecombinacion = atof(argv[9]);

	recombinacion = new double[10];
	ligamiento = new int[10];
	int lig [] = {-1,-1,-1,-1,-1,0,1,2,3,4};
	for (int i = 0; i < 10; i++) {
		recombinacion[i] = tasaRecombinacion;
		ligamiento[i] = lig[i];
	}

	//*** Hydroperiod length
	int hidroperiodo = 60

	//*** Diapausing eggs parameters
	//The diapusing egg banks exists by default, to omit it, it has to be set to False
	int huevDiap = 3;
	float tasaMixis = 0.7;
	float tasaSexAllocation = 0.5;
	// hatching rate has been corrected by dividing real rate 0.0347/(1-0.2721)
	float tasaEclosion = 0.0477;
	// survivor rate is 1 - (deterioration rate + loss rate)
	float tasaSupervivencia = 0.7279;
	// percentage of healthy eggs just after production
	float tasaSanos = 1.0;

	//*** auxiliary arrays for calculations
	int nGentpTot = (int) pow(2.0, (2.0 * (genNeut + genSelc)));
	int nGamtTot = (int) pow(2.0, 1.0 * (genNeut + genSelc));

	//*** auxiliary arrays used during sexual reproduction
	// individual genotypes in string
	string *genotipoString = new string[nGentpTot];
	for (int i = 0; i < nGentpTot; i++) {
		genotipoString[i] = convertirIntEnStringDiploide(i, genNeut + genSelc);
	}

	// gamete genotypes in string
	string *genotipoGametoString = new string[nGamtTot];
	for (int i = 0; i < nGamtTot; i++) {
		genotipoGametoString[i] = convertirIntEnStringHaploide(i, (genNeut
				+ genSelc));
	}
	// gametic combinations [gamete1][gamete2].
	int **genotipoCombinacionesGametos = new int*[nGamtTot];
	for (int i = 0; i < nGamtTot; i++)
		genotipoCombinacionesGametos[i] = new int[nGamtTot];
	for (int i = 0; i < nGamtTot; i++) {
		for (int j = 0; j < nGamtTot; j++) {
			string genotipoHuevoString = "";
			string gameto1 = genotipoGametoString[i];
			string gameto2 = genotipoGametoString[j];
			for (int x = 0; x < (genNeut + genSelc); x++) {
				genotipoHuevoString += gameto1[x];
				genotipoHuevoString += gameto2[x];
			}
			genotipoCombinacionesGametos[i][j] = binstring2int(genotipoHuevoString);
		}
	}


	//*** output description
	myfile << "\nDiapausing egg bank " << existeBanco
	 			 << "\nPopulation growth rate "<<tasaCrecPob
				 << "\nCarrying capacity " << capacidadCarga
				 << "\nNumber of migrants " << migrantes
				 << "\nNumber of founders " << fundadores
				 << "\nNumber of genes under selection " << genSelc
				 << "\nNumber of neutral genes "	<< genNeut
				 << "\nFitness component " << componenteEfic
				 << "\nGenetic recombination rate "<< tasaRecombinacion<<endl;

	for (int der = 0; der < 50; der++) {
		// *** source population
		PoblacionFuente *pf = new PoblacionFuente(genSelc, genNeut);

		//			PoblacionFuente *pf = new PoblacionFuente(genSelc, genNeut, mediaAlelSelec, mediaAlelNeut,
		//					varianzaAlelSelec, varianzaAlelNeut);

		// *** sink populations
		Poblacion *p1 = new Poblacion((int) (genNeut + genSelc),
				(float) tasaCrecPob, (float) componenteEfic,
				(float) tasaEclosion, (float) tasaSupervivencia,
				(float) tasaSanos, &recombinacion, &ligamiento, (int) huevDiap,
				(float) tasaMixis, (float) tasaSexAllocation, &genotipoString,
				&genotipoGametoString, &genotipoCombinacionesGametos,
				(int) fundadores, (int) migrantes);
		Poblacion *p2 = new Poblacion((int) (genNeut + genSelc),
				(float) tasaCrecPob, (float) componenteEfic,
				(float) tasaEclosion, (float) tasaSupervivencia,
				(float) tasaSanos, &recombinacion, &ligamiento, (int) huevDiap,
				(float) tasaMixis, (float) tasaSexAllocation, &genotipoString,
				&genotipoGametoString, &genotipoCombinacionesGametos,
				(int) fundadores, (int) migrantes);

		// *** source lake is created
		Laguna lf = Laguna(capacidadCarga, hidroperiodo, "fuente", pf);
		lf.fundarPoblacionFuente();

		// *** sink lakes are created and population growths assigned
		Laguna l1 = Laguna(capacidadCarga, hidroperiodo, "localidad1", p1, pf);
		Laguna l2 = Laguna(capacidadCarga, hidroperiodo, "localidad2", p2, pf);

		l1.establecerTasaCrecimientoGenotipo();
		l2.establecerTasaCrecimientoGenotipo();

		l1.getPoblacion()->setExisteBanco(existeBanco);
		l2.getPoblacion()->setExisteBanco(existeBanco);


		for (int i = 0; i < 4000; ++i) {
			if (i == 0) {
				//*** colonization from source lake
				lf.colonizar(&l1);
				lf.colonizar(&l2);
			}
			//						else if (i > 0 && migrantes > 0) {
			//				// *** Migramos desde la laguna fuente
			//				lf.migrarFuente(&l1);
			//				lf.migrarFuente(&l2);
			//						}
			if (((i+1)%100 == 0)|| i== 1 ) {
				calcularFstYQst(&lf, &l1, &l2, &resultadoFst);
				if (seleccion == false) {
					myfile << resultadoFst[1] << " ";
				} else {
					myfile << resultadoFst[0] << " " << resultadoFst[1] << " ";
				}
			}

			// *** asexual growth
			l1.crecerAsexual();
			l2.crecerAsexual();

			// *** sexual reproduction
			l1.getPoblacion()->reproduccionSexual();
			l2.getPoblacion()->reproduccionSexual();

			// *** survival in the bank
			l1.getPoblacion()->sobrevivirEnElBanco();
			l2.getPoblacion()->sobrevivirEnElBanco();

			// *** migration between sink lakes
			l1.migrarEntre();
			l2.migrarEntre();

			//*** hatching from the bank and sedimentation of the rest of eggs
			l1.getPoblacion()->eclosionEnElBanco();
			l2.getPoblacion()->eclosionEnElBanco();

			//*** hatching of migrants
			l1.eclosionarMigrantes(&l2);
			l2.eclosionarMigrantes(&l1);
		}
		delete pf;
		delete p1;
		delete p2;
	}
	myfile << endl;
	myfile.close();

	for (int a = 0; a < nGamtTot; a++)
		delete[] genotipoCombinacionesGametos[a];

	delete[] genotipoCombinacionesGametos;
	delete[] genotipoGametoString;
	delete[] genotipoString;
	delete[] ligamiento;
	delete[] recombinacion;
	delete[] resultadoFst;
}
