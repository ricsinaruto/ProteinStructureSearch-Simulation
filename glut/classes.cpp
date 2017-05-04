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

	/*for (int i = 0; i < 6; i++) {
	szomszedok.pop_back();
	}*/
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

void init_molekula::unset_ter_mol() {
	dronpa[x][y][z].ter = false;
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
	mol_szam = 0;
	vec_len = bemenetek_szama * 2 + 1 + molekulaSzam;
	genes = new double[vec_len];
	for (int i = 0; i < bemenetek_szama * 2; i++) {
		genes[i] = fRand(-DEF_MAX_TER, DEF_MAX_TER);
	}
	for (int i = bemenetek_szama * 2; i < vec_len - 1; i++) {
		if (fRand(0, 1) < 0.1) {
			genes[i] = int(fRand(0, 3 - 0.000000001));
			mol_szam++;
		}
		else genes[i] = 3;
	}

	bool good_one = false;
	while (!good_one) {
		genes[vec_len - 1] = int(fRand(0, molekulaSzam - 0.000000001));
		int temp = bemenetek_szama * 2 + genes[vec_len - 1];
		if (genes[temp] != 3) good_one = true;
	}

	fitness = 100000;
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
	mol_szam = 0;
	std::vector<int> indexek;
	// kitöröljük az előzőeket (mint egy destruktor)
	for (int i = 0; i < molekulaSzam; i++) {
		//kivéve legelső futás
		protein[i].delete_molekula();
		if (genes[bemenetek_szama * 2 + i] != 3) {
			indexek.push_back(i);
			mol_szam++;
		}
	}

	// létrehozzuk az újakat és rárakjuk a bemeneteket
	for (int i = 0; i < mol_szam; i++) {
		protein[i];
		// 10es számrendszer:
		int x = indexek[i] / 100;
		int y = (indexek[i] - x * 100) / 10;
		int z = indexek[i] - x * 100 - y * 10;
		protein[i].initialize_molekula(x + 12, y + 12, z + 12, false, 2, false);
	}

	// miután inicializáltunk kell egy első futás
	futas();

	//bemenetek definiálása
	for (int i = 0; i < mol_szam; i++) {
		protein[i].set_init_dipole();
		protein[i].set_desired();
		protein[i].set_actual();
		if (genes[bemenetek_szama * 2 + indexek[i]] == 0) {
			protein[i].ter = true;
			protein[i].set_ter_mol();
			protein[i].bemenet_szam = 0;
		}
		if (genes[bemenetek_szama * 2 + indexek[i]] == 1) {
			protein[i].ter = true;
			protein[i].set_ter_mol();
			protein[i].bemenet_szam = 1;
		}

		//1 kimenet
		if (genes[vec_len - 1] == indexek[i]) protein[i].kimenet = true;
	}



	SIMULATION(getFields(), false, mol_szam);
	fitness = fitness_func(mol_szam);

	indexek.clear();
}

// Crossover
DNA DNA::crossover(DNA partner) {
	// A new child
	DNA child;


	// Half from one, half from the other
	bool good_one = false;
	while (!good_one) {
		int midpoint1 = int(fRand(0, vec_len - 0.00000000001)); // Pick a midpoint
		int midpoint2 = int(fRand(0, vec_len - 0.00000000001)); // Pick a midpoint

		for (int i = 0; i < vec_len; i++) {
			if ((midpoint1 >midpoint2 && i < midpoint1 && i>midpoint2) || (midpoint1 <midpoint2 && i > midpoint1 && i<midpoint2)) child.genes[i] = genes[i];
			else              child.genes[i] = partner.genes[i];
		}

		//viable gene?
		int temp = bemenetek_szama * 2 + child.genes[vec_len - 1];
		if (child.genes[temp] != 3) good_one = true;
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

	int temp = bemenetek_szama * 2 + genes[vec_len - 1];
	for (int i = bemenetek_szama * 2; i < vec_len - 1; i++) {
		if (fRand(0, 1) < mutationRate) {
			bool good_one = false;
			while (!good_one) {
				if (fRand(0, 1) < 0.1) {
					genes[i] = int(fRand(0, 3 - 0.000000001));
				}
				else genes[i] = 3;

				if (genes[temp] != 3) good_one = true;
			}
		}
	}
	if (fRand(0, 1) < mutationRate) {
		bool good_one = false;
		while (!good_one) {
			genes[vec_len - 1] = int(fRand(0, molekulaSzam - 0.000000001));
			int temp = bemenetek_szama * 2 + genes[vec_len - 1];
			if (genes[temp] != 3) good_one = true;
		}
	}
}