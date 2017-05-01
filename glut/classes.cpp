#include "screencasts.h"

//kereső algoritmushoz molekula class függvények

//alap konstruktor
init_molekula::init_molekula() {
	x = 0;
	y = 0;
	z = 0;

	ter = false;
	kimenet = false;
	actual = new double[pow(2, bemenetek_szama)];
	torolt = false;
}

//konstruktor
void init_molekula::initialize_molekula(int _x, int  _y, int  _z, bool _ter, int _bemenet_szam, bool _kimenet) {
	x = _x;
	y = _y;
	z = _z;

	ter = _ter;
	bemenet_szam = _bemenet_szam;
	kimenet = _kimenet;

	dronpa[x][y][z].van = true;
	dronpa[x][y][z].dip = DEF_DIPOL;
	dronpa[x][y][z].dipA = DEF_DIPOL;
	dronpa[x][y][z].dipB = DEF_DIPOL;
	dronpa[x][y][z].ter = ter;
	torolt = false;
}


//delete molekula
void init_molekula::delete_molekula() {
	dronpa[x][y][z].van = false;
	dronpa[x][y][z].dip = 0;
	dronpa[x][y][z].dipA = 0;
	dronpa[x][y][z].dipB = 0;
	dronpa[x][y][z].qeA = 0;
	dronpa[x][y][z].qeB = 0;
	dronpa[x][y][z].qp1A = 0;
	dronpa[x][y][z].qp1B = 0;
	dronpa[x][y][z].qp2A = 0;
	dronpa[x][y][z].qp2B = 0;
	dronpa[x][y][z].ter = false;
	torolt = true;

	for (int i = 0; i < 6; i++) {
		szomszedok.pop_back();
	}
}

//set szomszédok
void init_molekula::set_szomszedok() {
	std::vector<int> i;
	i.push_back(x + 1);
	i.push_back(y);
	i.push_back(z);
	szomszedok.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x - 1);
	i.push_back(y);
	i.push_back(z);
	szomszedok.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y + 1);
	i.push_back(z);
	szomszedok.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y - 1);
	i.push_back(z);
	szomszedok.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y);
	i.push_back(z + 1);
	szomszedok.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y);
	i.push_back(z - 1);
	szomszedok.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();
}

//tér set
void init_molekula::set_ter_mol() {
	dronpa[x][y][z].ter = true;
}

//dipól lekérése
double init_molekula::get_dipole() {
	return dronpa[x][y][z].dip;
}

//init_dipole set
void init_molekula::set_init_dipole() {
	init_dipole = dronpa[x][y][z].dip;
}

//desired vektor megadása
void init_molekula::set_desired() {
	desired[0] = init_dipole - tolerance;
	desired[1] = init_dipole + tolerance;
}

//tér nagyság megadása
void init_molekula::set_ter(double terMag) {
	dronpa[x][y][z].ter = true;
	dronpa[x][y][z].terMag = terMag;
}

//set actual dipole value
void init_molekula::set_actual() {
	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		actual[i] = init_dipole;
	}
}

//update actual dipole value
void init_molekula::update_actual(int sor) {
	actual[sor] = dronpa[x][y][z].dip;
}

//reset dipole moment
void init_molekula::reset_dipole(double dipole) {
	dronpa[x][y][z].dip = dipole;
	dronpa[x][y][z].dipA = dipole;
	dronpa[x][y][z].dipB = dipole;
	dronpa[x][y][z].qeA = 0;
	dronpa[x][y][z].qeB = 0;
	dronpa[x][y][z].qp1A = 0;
	dronpa[x][y][z].qp1B = 0;
	dronpa[x][y][z].qp2A = 0;
	dronpa[x][y][z].qp2B = 0;
}



/* GENETIC ALGORITHM */

//constructor
DNA::DNA() {
	genes = new double[bemenetek_szama * 2];
	for (int i = 0; i < bemenetek_szama * 2; i++) {
		genes[i] = fRand(-DEF_MAX_TER, DEF_MAX_TER);
	}
	fitness = 1000;
	hasonlit = false;
}

//get fields from genes
double **DNA::getFields() {
	double **fields = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { fields[i] = new double[2]; }

	for (int i = 0; i < bemenetek_szama; i++) {
		fields[i][0] = genes[2 * i];
		fields[i][1] = genes[2 * i + 1];
	}
	return fields;
}

//calculate the fitness
void DNA::calcFitness() {
	SIMULATION(getFields(), false, molekulaSzam);
	fitness = fitness_func(molekulaSzam);
}

// Crossover
DNA DNA::crossover(DNA partner) {
	// A new child
	DNA child;

	int midpoint = int(fRand(0, bemenetek_szama * 2 - 0.00000000001)); // Pick a midpoint

																	   // Half from one, half from the other
	for (int i = 0; i < bemenetek_szama * 2; i++) {
		if (i > midpoint) child.genes[i] = genes[i];
		else              child.genes[i] = partner.genes[i];
	}
	return child;
}

// Based on a mutation probability, picks a new random character
void DNA::mutate(float mutationRate) {
	for (int i = 0; i < bemenetek_szama * 2; i++) {
		if (fRand(0, 1) < mutationRate) {
			genes[i] = fRand(-DEF_MAX_TER, DEF_MAX_TER);
		}
	}
}