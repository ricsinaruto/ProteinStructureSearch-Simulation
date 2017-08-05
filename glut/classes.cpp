#include "main.h"

/* PROTEIN CLASS */

// basic constructor
init_molecule::init_molecule() {
	x = 0;
	y = 0;
	z = 0;

	hasField = false;
	isOutput = false;
	actual = new double[pow(2, num_in)];
	deleted = false;
}

// constructor
void init_molecule::initialize_molecule(int _x, int  _y, int  _z, bool _hasField, int _input_num, bool _isOutput) {
	x = _x;
	y = _y;
	z = _z;

	hasField = _hasField;
	input_num = _input_num;
	isOutput = _isOutput;

	dronpa[x][y][z].is = true;
	dronpa[x][y][z].dip = DEF_DIPOL;
	dronpa[x][y][z].dipA = DEF_DIPOL;
	dronpa[x][y][z].dipB = DEF_DIPOL;
	dronpa[x][y][z].field = hasField;
	deleted = false;
}

// delete molecule
void init_molecule::delete_molecule() {
	dronpa[x][y][z].is = false;
	dronpa[x][y][z].dip = 0;
	dronpa[x][y][z].dipA = 0;
	dronpa[x][y][z].dipB = 0;
	dronpa[x][y][z].qeA = 0;
	dronpa[x][y][z].qeB = 0;
	dronpa[x][y][z].qp1A = 0;
	dronpa[x][y][z].qp1B = 0;
	dronpa[x][y][z].qp2A = 0;
	dronpa[x][y][z].qp2B = 0;
	dronpa[x][y][z].field = false;
	deleted = true;
}

// set neighbours
void init_molecule::set_neighbours() {
	std::vector<int> i;

	i.push_back(x + 1);
	i.push_back(y);
	i.push_back(z);
	neighbours.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x - 1);
	i.push_back(y);
	i.push_back(z);
	neighbours.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y + 1);
	i.push_back(z);
	neighbours.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y - 1);
	i.push_back(z);
	neighbours.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y);
	i.push_back(z + 1);
	neighbours.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();

	i.push_back(x);
	i.push_back(y);
	i.push_back(z - 1);
	neighbours.push_back(i);
	for (int j = 0; j < 3; j++) i.pop_back();
}

// set field
void init_molecule::set_field_mol() {
	dronpa[x][y][z].field = true;
}

// unset field
void init_molecule::unset_field_mol() {
	dronpa[x][y][z].field = false;
}

// get dipole value
double init_molecule::get_dipole() {
	return dronpa[x][y][z].dip;
}

// set initial dipole value
void init_molecule::set_init_dipole() {
	init_dipole = dronpa[x][y][z].dip;
}

// set desired vector
void init_molecule::set_desired() {
	desired[0] = init_dipole - tolerance;
	desired[1] = init_dipole + tolerance;
}

// set field magnitude
void init_molecule::set_field(double fieldMag) {
	dronpa[x][y][z].field = true;
	dronpa[x][y][z].fieldMag = fieldMag;
}

// set actual dipole value
void init_molecule::set_actual() {
	for (int i = 0; i < pow(2, num_in); i++) {
		actual[i] = init_dipole;
	}
}

// update actual dipole value
void init_molecule::update_actual(int row) {
	actual[row] = dronpa[x][y][z].dip;
}

//reset dipole moment
void init_molecule::reset_dipole(double dipole) {
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

// constructor
DNA::DNA() {
	mol_num = 0;
	vec_len = num_in * 2 + 1 + num_molecules;
	genes = new double[vec_len];
	for (int i = 0; i < num_in * 2; i++) {
		genes[i] = fRand(-DEF_MAX_FIELD, DEF_MAX_FIELD);
	}
	for (int i = num_in * 2; i < vec_len - 1; i++) {
		if (fRand(0, 1) < 0.1) {
			genes[i] = int(fRand(0, 3 - 0.000000001));
			mol_num++;
		}
		else genes[i] = 3;
	}

	bool good_one = false;
	while (!good_one) {
		genes[vec_len - 1] = int(fRand(0, num_molecules- 0.000000001));
		int temp = num_in * 2 + genes[vec_len - 1];
		if (genes[temp] != 3) good_one = true;
	}

	fitness = 100000;
	similar = false;
}

// get field values from genes
double **DNA::getFields() {
	double **fields = new double*[num_in];
	for (int i = 0; i < num_in; i++) { fields[i] = new double[2]; }

	for (int i = 0; i < num_in; i++) {
		fields[i][0] = genes[2 * i];
		fields[i][1] = genes[2 * i + 1];
	}
	return fields;
}

// calculate the fitness
void DNA::calcFitness() {
	mol_num = 0;
	std::vector<int> indexes;

	// we delete the previous (like a destructor)
	for (int i = 0; i < num_molecules; i++) {
		//except first run
		protein[i].delete_molecule();
		if (genes[num_in * 2 + i] != 3) {
			indexes.push_back(i);
			mol_num++;
		}
	}

	// initialize new ones 
	for (int i = 0; i < mol_num; i++) {
		protein[i];
		// base 10:
		int x = indexes[i] / 100;
		int y = (indexes[i] - x * 100) / 10;
		int z = indexes[i] - x * 100 - y * 10;
		protein[i].initialize_molecule(x + 12, y + 12, z + 12, false, 2, false);
	}

	// after initialization we need a first run
	run();

	// define the inputs
	for (int i = 0; i < mol_num; i++) {
		protein[i].set_init_dipole();
		protein[i].set_desired();
		protein[i].set_actual();

		if (genes[num_in * 2 + indexes[i]] == 0) {
			protein[i].hasField = true;
			protein[i].set_field_mol();
			protein[i].input_num = 0;
		}
		if (genes[num_in * 2 + indexes[i]] == 1) {
			protein[i].hasField = true;
			protein[i].set_field_mol();
			protein[i].input_num = 1;
		}

		// only 1 output
		if (genes[vec_len - 1] == indexes[i]) protein[i].isOutput = true;
	}

	// do a simulation loop
	SIMULATION(getFields(), false, mol_num);
	fitness = fitness_func(mol_num);

	indexes.clear();
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
			if ((midpoint1 >midpoint2 && i < midpoint1 && i>midpoint2) ||
				(midpoint1 <midpoint2 && i > midpoint1 && i<midpoint2)) 
				 child.genes[i] = genes[i];
			else child.genes[i] = partner.genes[i];
		}

		//viable gene?
		int viable = num_in * 2 + child.genes[vec_len - 1];
		if (child.genes[viable] != 3) good_one = true;
	}

	return child;
}

// Based on a mutation probability picks a new random parameter
void DNA::mutate(float mutationRate) {
	// mutate field values separately
	for (int i = 0; i < num_in * 2; i++) {
		if (fRand(0, 1) < mutationRate) {
			genes[i] = fRand(-DEF_MAX_FIELD, DEF_MAX_FIELD);
		}
	}

	// mutate input molecules
	int temp = num_in * 2 + genes[vec_len - 1];
	for (int i = num_in * 2; i < vec_len - 1; i++) {
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

	// pick a new output molecule
	if (fRand(0, 1) < mutationRate) {
		bool good_one = false;
		while (!good_one) {
			genes[vec_len - 1] = int(fRand(0, num_molecules - 0.000000001));
			int temp = num_in * 2 + genes[vec_len - 1];
			if (genes[temp] != 3) good_one = true;
		}
	}
}