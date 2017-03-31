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
	szomszedok = new double[DEF_PROTEIN_NUMBER];
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
