#include "main.h"
using namespace std;

// specify outputs, runs if the k key is pressed
void set_outputs() {
	num_in = log2(counter - 1);
	num_out = numbers[0];
	for (int i = 1; i < counter; i++) {
		outputs[numbers[0] - 1][i - 1] = numbers[i];

		//set back to initial value
		numbers[i] = 36;
	}

	// reset
	numbers[0] = 36;
	counter = 0;
}

// this runs if the r key is pressed. it will run a simulation step for the specified structure
// can also save to file, defined in input arguments
void run_old(char* file_name, bool first_open) {
	ofstream fileOut;
	ifstream fileIn;
	int i, j, k, l, n;
	string line;
	n = 50 / dt;	// number of simulation steps
	int count = 0;

	// set file open modes
	if (first_open) fileOut.open("tmp.csv");
	else fileOut.open("tmp.csv", ios::app);
	if (first_open) fileOut << ",";

	for (i = 1; i <= 36; i++) {
		for (j = 1; j <= 36; j++) {
			for (k = 1; k <= 36; k++) {
				if (dronpa[i][j][k].is) {
					if (first_open) fileOut << i - 1 << " " << j - 1 << " " << k - 1 << ",";
					count++;
				}
			}
		}
	}
	if (first_open) fileOut << endl;

	int *iArray = new int[count];
	int *jArray = new int[count];
	int *kArray = new int[count];
	count = 0;

	// filter out the proteins that should be used in simulation
	for (i = 1; i <= 36; i++) {
		for (j = 1; j <= 36; j++) {
			for (k = 1; k <= 36; k++) {
				if (dronpa[i][j][k].is) {
					iArray[count] = i;
					jArray[count] = j;
					kArray[count] = k;
					count++;
				}
			}
		}
	}

	/* Simulation loop */
	for (l = 0; l < n; l++) {
		t += dt;
		fileOut << t << ",";

		for (i = 0; i < count; i++) {
			// step A
			if (l % 2 == 0) {
				U = K / dist*(dronpa[iArray[i] - 1][jArray[i]][kArray[i]].dipB + dronpa[iArray[i] + 1][jArray[i]][kArray[i]].dipB +
					dronpa[iArray[i]][jArray[i] - 1][kArray[i]].dipB + dronpa[iArray[i]][jArray[i] + 1][kArray[i]].dipB +
					dronpa[iArray[i]][jArray[i]][kArray[i] - 1].dipB + dronpa[iArray[i]][jArray[i]][kArray[i] + 1].dipB);

				if (dronpa[iArray[i]][jArray[i]][kArray[i]].field) U += dronpa[iArray[i]][jArray[i]][kArray[i]].fieldMag;

				dronpa[iArray[i]][jArray[i]][kArray[i]].qeA = dronpa[iArray[i]][jArray[i]][kArray[i]].qeB +
					(U / Ce1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qeB / (Ce1*Ce2))*dt;

				if (U > (dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A = dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B;

				if (U < (dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A = dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B;

				dronpa[iArray[i]][jArray[i]][kArray[i]].dip = dronpa[iArray[i]][jArray[i]][kArray[i]].dipA = -100 +
					dronpa[iArray[i]][jArray[i]][kArray[i]].qeA + dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A + 
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A;
			}

			// step B
			else {
				U = K / dist*(dronpa[iArray[i] - 1][jArray[i]][kArray[i]].dipA + dronpa[iArray[i] + 1][jArray[i]][kArray[i]].dipA +
					dronpa[iArray[i]][jArray[i] - 1][kArray[i]].dipA + dronpa[iArray[i]][jArray[i] + 1][kArray[i]].dipA +
					dronpa[iArray[i]][jArray[i]][kArray[i] - 1].dipA + dronpa[iArray[i]][jArray[i]][kArray[i] + 1].dipA);

				if (dronpa[iArray[i]][jArray[i]][kArray[i]].field) U += dronpa[iArray[i]][jArray[i]][kArray[i]].fieldMag;

				dronpa[iArray[i]][jArray[i]][kArray[i]].qeB = dronpa[iArray[i]][jArray[i]][kArray[i]].qeA +
					(U / Ce1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qeA / (Ce1*Ce2))*dt;

				if (U > (dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B = dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A;

				if (U < (dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B = dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A;

				dronpa[iArray[i]][jArray[i]][kArray[i]].dip = dronpa[iArray[i]][jArray[i]][kArray[i]].dipB = -100 +
					dronpa[iArray[i]][jArray[i]][kArray[i]].qeB + dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B + 
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B;
			}
			fileOut << dronpa[iArray[i]][jArray[i]][kArray[i]].dip << ",";
		}
		fileOut << endl;
	}
	fileOut.close();

	fileIn.open("tmp.csv");
	fileOut.open(file_name);
	while (getline(fileIn, line)) {
		for (i = 0; i < line.length(); i++) {
			if (line[i] == '.') fileOut << ".";
			else fileOut << line[i];
		}
		fileOut << endl;
	}
	fileIn.close();
	fileOut.close();

	delete[] iArray;
	delete[] jArray;
	delete[] kArray;
}

// same as the function before, but more fields can be applied at the same time, used for searching algo
void run() {
	int i, j, k, l, n;
	string sor;
	n = 50 / dt;
	int count = 0;

	for (i = 1; i <= 36; i++) {
		for (j = 1; j <= 36; j++) {
			for (k = 1; k <= 36; k++) {
				if (dronpa[i][j][k].is) {
					count++;
				}
			}
		}
	}

	int *iArray = new int[count];
	int *jArray = new int[count];
	int *kArray = new int[count];
	count = 0;

	for (i = 1; i <= 36; i++) {
		for (j = 1; j <= 36; j++) {
			for (k = 1; k <= 36; k++) {
				if (dronpa[i][j][k].is) {
					iArray[count] = i;
					jArray[count] = j;
					kArray[count] = k;
					count++;
				}
			}
		}
	}

	// run simulation step
	for (l = 0; l < n; l++) {
		for (i = 0; i < count; i++) {
			// step A
			if (l % 2 == 0) {
				U = K / dist*(dronpa[iArray[i] - 1][jArray[i]][kArray[i]].dipB + dronpa[iArray[i] + 1][jArray[i]][kArray[i]].dipB +
					dronpa[iArray[i]][jArray[i] - 1][kArray[i]].dipB + dronpa[iArray[i]][jArray[i] + 1][kArray[i]].dipB +
					dronpa[iArray[i]][jArray[i]][kArray[i] - 1].dipB + dronpa[iArray[i]][jArray[i]][kArray[i] + 1].dipB);

				if (dronpa[iArray[i]][jArray[i]][kArray[i]].field) U += dronpa[iArray[i]][jArray[i]][kArray[i]].fieldMag;

				dronpa[iArray[i]][jArray[i]][kArray[i]].qeA = dronpa[iArray[i]][jArray[i]][kArray[i]].qeB +
					(U / Ce1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qeB / (Ce1*Ce2))*dt;

				if (U > (dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A = dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B;

				if (U < (dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A = dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B;

				dronpa[iArray[i]][jArray[i]][kArray[i]].dip = dronpa[iArray[i]][jArray[i]][kArray[i]].dipA = -100 +
					dronpa[iArray[i]][jArray[i]][kArray[i]].qeA + dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A +
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A;
			}

			// step B
			else {
				U = K / dist*(dronpa[iArray[i] - 1][jArray[i]][kArray[i]].dipA + dronpa[iArray[i] + 1][jArray[i]][kArray[i]].dipA +
					dronpa[iArray[i]][jArray[i] - 1][kArray[i]].dipA + dronpa[iArray[i]][jArray[i] + 1][kArray[i]].dipA +
					dronpa[iArray[i]][jArray[i]][kArray[i] - 1].dipA + dronpa[iArray[i]][jArray[i]][kArray[i] + 1].dipA);

				if (dronpa[iArray[i]][jArray[i]][kArray[i]].field) U += dronpa[iArray[i]][jArray[i]][kArray[i]].fieldMag;

				dronpa[iArray[i]][jArray[i]][kArray[i]].qeB = dronpa[iArray[i]][jArray[i]][kArray[i]].qeA +
					(U / Ce1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qeA / (Ce1*Ce2))*dt;

				if (U > (dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B = dronpa[iArray[i]][jArray[i]][kArray[i]].qp1A;

				if (U < (dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A / Cp2)) dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B =
					dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A + (U / Cp1 - dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A /
					(Cp1*Cp2))*dt;
				else dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B = dronpa[iArray[i]][jArray[i]][kArray[i]].qp2A;

				dronpa[iArray[i]][jArray[i]][kArray[i]].dip = dronpa[iArray[i]][jArray[i]][kArray[i]].dipB = -100 +
					dronpa[iArray[i]][jArray[i]][kArray[i]].qeB + dronpa[iArray[i]][jArray[i]][kArray[i]].qp1B + dronpa[iArray[i]][jArray[i]][kArray[i]].qp2B;
			}
		}
	}

	delete[] iArray;
	delete[] jArray;
	delete[] kArray;
}

// compare in1 and in2
bool compare(int in1, int in2, int comp) {
	bool isGood = false;

	if ((in1) < in2 && comp == 0) isGood = true;
	if (in1>(in2) && comp == 1) isGood = true;
	if (comp == 2) isGood = true;

	return isGood;
}

// return the factorial of input number
int factorial(int f) {
	if (f == 0) return 1;
	return(f * factorial(f - 1));
}

// random number generator
double fRand(double fMin, double fMax) {
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

// checks all the rows of a logical function, returns whether they are good
bool checkLogic(int mol_numb) {
	bool good = true;
	int output_it = 0;

	for (int i = 0; i < pow(2, num_in); i++) {
		for (int j = 0; j < mol_numb; j++) {
			if (good && protein[j].isOutput) {
				if (!compare(protein[j].actual[i], protein[j].desired[outputs[output_it][i]], 
					outputs[output_it][i])) good = false;
				output_it++;
			}
		}
		output_it = 0;
	}

	return good;
}

// calculates fitness function
double fitness_func(int mol_numb) {
	double fitness = 0;
	int output_it = 0;
	bool same_dipole = true;

	for (int i = 0; i < pow(2, num_in); i++) {
		for (int j = 0; j < mol_numb; j++) {
			if (protein[j].isOutput) {
				if (protein[j].actual[i] > protein[j].init_dipole + 0.01 || protein[j].actual[i] < protein[j].init_dipole - 0.01)
					same_dipole = false;

				if (outputs[output_it][i]) {
					if (protein[j].actual[i] < protein[j].desired[outputs[output_it][i]] + OVER_FIT) {
						fitness += pow(protein[j].desired[outputs[output_it][i]] + OVER_FIT - protein[j].actual[i], 2);
					}
					output_it++;
				}
				else if (!outputs[output_it][i]) {
					if (protein[j].actual[i] > protein[j].desired[outputs[output_it][i]] - OVER_FIT) {
						fitness += pow(protein[j].desired[outputs[output_it][i]] - OVER_FIT - protein[j].actual[i], 2);
					}
					output_it++;
				}
			}
		}
		output_it = 0;
	}

	// filter the case when an output has no neighbours 
	if (same_dipole) fitness = 10000;
	fitness = fitness / 100;
	fitness = 1 / (1 + log(1 + fitness));
	return fitness;
}

// runs some simulation steps, if save is true then it saves the results to .csv file 
void SIMULATION(double **field_vec, bool save, int mol_numb) {
	for (int i = 0; i < pow(2, num_in); i++) {
		for (int j = 0; j < mol_numb; j++) {
			protein[j].reset_dipole(protein[j].init_dipole);
		}

		//apply field
		for (int j = 0; j < mol_numb; j++) {
			if (protein[j].hasField) {
				protein[j].set_field(field_vec[protein[j].input_num][inputs[i][protein[j].input_num]]);
			}
		}

		// save to .csv
		std::string ok = "graf" + std::to_string(i) + ".csv";
		char* c = &ok[0];
		if (save) run_old(c, true);
		else run();


		for (int j = 0; j < mol_numb; j++) {
			if (protein[j].hasField) {
				protein[j].set_field(0);
			}
		}
		if (save) run_old(c, false);
		else run();

		// save the dipole moments
		for (int j = 0; j <mol_numb; j++) {
			if (protein[j].isOutput) {
				protein[j].update_actual(i);
			}
		}
	}
}

// initialize protein objects
void protein_def() {
	for (int k = 0; k < num_molecules; k++) {
		protein[k];
	}

	// FOR TESTING ============ //
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				protein[k + j * 3 + i * 9].initialize_molecule(18 + i, 18 + j, 18 + k, false, 2, false);
			}
		}
	}

	// first run only once, and save the initial dipole values
	run();
	for (int i = 0; i < num_molecules; i++) {
		protein[i].set_init_dipole();
		protein[i].set_desired();
		protein[i].set_actual();
	}
	// ========================= //
}

// main searching function based on GA
bool search_algo() {
	num_molecules = DEF_PROTEIN_NUMBER;

	// initialize random field, first index is the number of the input molecule, 
	// second index is to specify whether we are talking about a field representing a 0 or 1 logic bit
	double **inputField = new double*[num_in];
	for (int i = 0; i < num_in; i++) { inputField[i] = new double[2]; }

	// GA params
	int generation = 1;
	bool compare = false;

	int n = ITER_NUMBER;		// number of children
	int stuff = 1;				// counter for while cycle
	int for_i;					// for for cycles
	double x, y;				// these will be random numbers

	// initialize original candidates
	for (int i = 0; i < num_in; i++) {
		for (int j = 0; j < 2; j++) {
			inputField[i][j] = START_POINT;
		}
	}

	// best candidate	
	double **bestField = new double*[num_in];
	for (int i = 0; i < num_in; i++) { bestField[i] = new double[2]; }

	for (int i = 0; i < num_in; i++) {
		for (int j = 0; j < 2; j++) {
			bestField[i][j] = inputField[i][j];
		}
	}

	// more GA params
	DNA population[DEF_POP_SIZE];
	double mutationRate = DEF_MUT_RATE;
	std::vector<DNA> matingPool;
	double last_best = 0;

	// search loop
	while (!compare && n>0) {
		// delete previous stuff, except first run
		if (n == ITER_NUMBER) {
			for (int k = 0; k < num_molecules; k++) {
				protein[k];
				protein[k].initialize_molecule(18, 18, 18, false, 2, false);
			}
		}

		double best_fitness = 0;
		for (int i = 0; i < DEF_POP_SIZE; i++) {
			population[i].calcFitness();

			if (population[i].fitness > best_fitness) best_fitness = population[i].fitness;
			// also check logikai_hasonlitas
			compare = checkLogic(population[i].mol_num);

			if (compare) {
				num_molecules = population[i].mol_num;
				bestField = population[i].getFields();
				i = DEF_POP_SIZE;
			}
		}

		if (!compare) {
			// create the mutation pool
			for (int i = 0; i < DEF_POP_SIZE; i++) {
				int nn = int(population[i].fitness * DEF_MATING_POOL_COEFF);
				for (int j = 0; j < nn; j++) {
					matingPool.push_back(population[i]);
				}
			}

			// do mating and mutation
			for (int i = 0; i < DEF_POP_SIZE; i++) {
				int a = int(fRand(0, matingPool.size() - 0.00000000001));
				int b = int(fRand(0, matingPool.size() - 0.00000000001));

				DNA partnerA = matingPool[a];
				DNA partnerB = matingPool[b];
				DNA child = partnerA.crossover(partnerB);

				child.mutate(mutationRate);
				population[i] = child;
			}
		}

		// check the params and write them
		if (n % 50 == 0) {
			if (last_best - 0.01 < best_fitness && best_fitness < last_best + 0.01)  mutationRate = 0.3;
			else mutationRate = DEF_MUT_RATE;

			last_best = best_fitness;
		}

		cout << best_fitness << "         " << mutationRate << endl;
		matingPool.clear();

		// iteration counters
		generation++;
		n--;
	}

	/* Write the data out */
	if (compare) {
		for (int i = 0; i < num_in; i++) {
			for (int j = 0; j < 2; j++) {
				cout << "Input molecule "<<i  << ", magnitude of logic field "<<j <<" : "<< bestField[i][j] << "   ";
			}
			cout << endl;
		}
	}
	cout << "number of simulations: " << DEF_POP_SIZE*(ITER_NUMBER - n) << endl;

	// delete everything
	for (int i = 0; i <num_in; i++) { delete[] bestField[i]; }
	delete[] bestField;
	for (int i = 0; i <num_in; i++) { delete[] inputField[i]; }
	delete[] inputField;

	return compare;
}

// this takes the longest time, inside it is the search_algo function with the main searching algorithm
void main_func() {
	bool success = false;

	while (!success) {
		// returns a bool, whether the structure needed was found
		success = search_algo();

		// write out the data
		if (success) {
			cout << "Number of molecules: " << num_molecules << endl;
			for (int j = 0; j < num_molecules; j++) {
				if (protein[j].hasField) {
					cout << protein[j].input_num << "<-applied this field number to: "
						<< protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
				}

				if (protein[j].isOutput) {
					cout << "output: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << "    ";
					cout << "dipole: ";
					for (int r = 0; r < pow(2, num_in); r++) {
						cout << protein[j].actual[r] << " ";
					}
					cout << endl;
				}

				if (!protein[j].hasField && !protein[j].isOutput) {
					cout << "other molecules: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
				}
			}
		}
	}
}

// save a protein structure to file
void save() {
	ofstream strukt;
	strukt.open("strukt.csv");
	int i, j, k;

	for (i = 1; i <= 36; i++) {
		for (j = 1; j <= 36; j++) {
			for (k = 1; k <= 36; k++) {
				if (dronpa[i][j][k].is) {

					strukt << ";" << i - 1 << ";" << j - 1 << ";" << k - 1 << ";";

					if (dronpa[i][j][k].field) strukt << 1 << ";";
					else strukt << 0 << ";";
					strukt << endl;
				}
			}
		}
	}
	strukt.close();
	strukt.open("params.dat");

	strukt << th << endl;
	strukt << ph << endl;
	strukt << dim << endl;

	strukt.close();
}

// load a protein structure from file
void load() {
	ifstream fileIn;
	string line;
	fileIn.open("strukt.csv");

	int number = -1;
	int numbers_l[4] = { 0,0,0,0 };
	bool param = false;
	int counter_l = 2;

	while (getline(fileIn, line)) {
		int comma = 0;
		int position = 0;
		int linePos = 0;

		if (param) {
			int az = 0;
			int ss = 2;

			if (line.length() == 2) ss = 1;
			for (int i = 0; i < line.length(); i++) {
				switch (line.length()) {
				case 1: az = line[i] - 48;
					break;
				case 2: az += (line[i] - 48)*pow(10, ss);
					break;
				case 3: az += (line[i] - 48)*pow(10, ss);
					break;
				}
				ss--;
			}
			switch (counter_l) {
			case 2: // th = az;
				break;
			case 1: // ph = az;
				break;
			case 0: // dim = az;
				break;
			}
			counter_l--;
		}

		if (line != "." && !param) {
			for (int i = 0; i < line.length(); i++) {
				if (line[i] == ';') {
					if (i - comma == 2) {
						numbers_l[linePos] = (int)line[i - 1] - 48;
						linePos++;
					}

					if (i - comma == 3) {
						numbers_l[linePos] = ((int)line[i - 2] - 48) * 10 + (int)line[i - 1] - 48;
						linePos++;
					}
					comma = i;
				}
			}

			dronpa[numbers_l[0] + 1][numbers_l[1] + 1][numbers_l[2] + 1].is = true;
			if (numbers_l[3] == 1) dronpa[numbers_l[0] + 1][numbers_l[1] + 1][numbers_l[2] + 1].field = true;

			// switch to base 36
			number = numbers_l[0] * 1296 + numbers_l[1] * 36 + numbers_l[2];
			for (int k = -18; k < 18; k++) {
				for (int j = -18; j < 18; j++) {
					for (int i = -18; i < 18; i++) {
						// draw the cubes
						if (number == position) {
							cubes[number] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 } ,
							{ (float)numbers[0],(float)numbers[1],(float)numbers[2] } } };
						}
						position++;
					}
				}
			}
		}
		if (line == "." && !param) param = true;
	}
	fileIn.close();

	fileIn.open("params.dat");
	line = "";
	counter_l = 0;
	std::string::size_type sz;   // alias of size_t

	while (getline(fileIn, line)) {
		switch (counter_l) {
		case 0: th = stoi(line, &sz);
			break;
		case 1: ph = stoi(line, &sz);
			break;
		case 2: dim = stoi(line, &sz);
			break;
		}
		counter_l++;
	}
	fileIn.close();
}


/*  HANDLE KEYBOARD AND MOUSE INPUTS  */

// at the press of q,w,e handle movement and rotation
void windowPmotion(int x, int y) {
	mouseX = x;
	mouseY = y;

	// rotate the light
	if (Shift == "light") {
		if (mouseBtnPressed == "Left") {
			xcoord = mouseX;
			ycoord = mouseY;
			lightTh2 = lightTh;
			lightPh2 = lightPh;
			ecX2 = ecX;
			ecY2 = ecY;
			th2 = th;
			ph2 = ph;
		}

		if (mouseBtnPressed == "Right") {
			lightTh = (lightTh2 + (mouseX - xcoord));
			lightPh = (lightPh2 + (mouseX - xcoord));
		}
	}

	// rotate structure
	else if (Shift == "rotation") {
		if (mouseBtnPressed == "Left") {
			xcoord = mouseX;
			ycoord = mouseY;
			lightTh2 = lightTh;
			lightPh2 = lightPh;
			th2 = th;
			ph2 = ph;
			ecX2 = ecX;
			ecY2 = ecY;
		}

		if (mouseBtnPressed == "Right") {
			th = (th2 + (mouseX - xcoord));
			ph = (ph2 + (mouseY - ycoord));
		}
	}

	// move structure
	else if (Shift == "movement") {
		if (mouseBtnPressed == "Left") {
			xcoord = mouseX;
			ycoord = mouseY;
			lightTh2 = lightTh;
			lightPh2 = lightPh;
			ecX2 = ecX;
			ecY2 = ecY;
			th2 = th;
			ph2 = ph;
		}

		if (mouseBtnPressed == "Right") {
			ecX = ecX2 + (-mouseX + xcoord) / 10;
			ecY = ecY2 + (mouseY - ycoord) / 10;
		}
	}

	th %= 360;
	ph %= 360;
	redisplayAll();
}

// zoom
void mouseWheel(int scroll, int dir, int x, int y) {
	dim -= (double)dir;
	redisplayAll();
}

// handle mouse button presses
void windowMouse(int btn, int state, int x, int y) {
	if (btn == GLUT_LEFT_BUTTON) mouseBtnPressed = "Left";
	else if (btn == GLUT_RIGHT_BUTTON) mouseBtnPressed = "Right";

	if (state == GLUT_UP) mouseState = "up";
	if (state == GLUT_DOWN) mouseState = "down";
	
	redisplayAll();
}

// handle all keyboard presses (almost)
void windowKey(unsigned char key, int x, int y) {
	/* Exit on ESC */
	if (key == 27) exit(0);

	// space
	if (key == 32) {
		if (switcher == "parameters") switcher = "coordinates";
		else switcher = "parameters";
	}

	// define parameters
	if (switcher == "parameters") {
		// toggle axes and parameter display
		if (key == 'x' || key == 'X') toggleAxes = 1 - toggleAxes;
		else if (key == 'v' || key == 'V') toggleParams = 1 - toggleParams;

		/*  Toggle lighting */
		else if (key == 'l' || key == 'L') toggleLight = 1 - toggleLight;
		/*  Ambient level */
		else if (key == 'a' && ambient>0) ambient -= 5;
		else if (key == 'A' && ambient<100) ambient += 5;
		/*  Diffuse level */
		else if (key == 'd' && diffuse>0) diffuse -= 5;
		else if (key == 'D' && diffuse<100) diffuse += 5;

		// movement, rotation, light rotation
		else if (key == 'q') Shift = "movement";
		else if (key == 'w') Shift = "rotation";
		else if (key == 'e') Shift = "light";

		// set outputs
		else if (key == 'k') {
			enter = "output"; 
			set_outputs();
		}

		else if (key == 'r') run_old("graf.csv", false);// run simulation
		else if (key == 'f') main_func();				// run search algorithm
		else if (key == 't') enter = "field";			// switch enter mode to field
		else if (key == 's') save();					// save structure
		else if (key == 'o') load();					// load structure
	}

	/* Change field of view angle */
	if (key == '-' && key>1) fov--;
	if (key == '+' && key<179) fov++;

	/* Light elevation */
	if (key == '[') lightY -= 0.5;
	if (key == ']') lightY += 0.5;

	// enter and delete keys
	if (key == 13)	enter = "pressed";
	if (key == 8)  enter = "delete";

	// give alphanumeric coordinates
	if (switcher == "coordinates") {
		for (int i = 97; i < 123; i++) {
			if (key == i) { numbers[counter] = i - 87; counter++; key = ','; }
		}
		for (int i = 48; i < 58; i++) {
			if (key == i) { numbers[counter] = i - 48; counter++; key = ','; }
		}
	}

	// coordinates from 0 to 9
	for (int i = 48; i < 58; i++) {
		if (key == i) { numbers[counter] = i - 48; counter++; key = ','; }
	}

	/* Translate shininess power to value (-1 => 0) */
	shinyvec[0] = shininess<0 ? 0 : pow(2.0, shininess);
	redisplayAll();
}

/* Window menu is the same as the keyboard clicks */
void windowMenu(int value) {
	windowKey((unsigned char)value, 0, 0);
}

/* GLUT calls this routine when an arrow key is pressed */
void windowSpecial(int key, int x, int y) {
	int modifiers = glutGetModifiers();
	/* If holding shift, then rotate/elevate */
	if (modifiers == GLUT_ACTIVE_SHIFT) {
		/*  Right/Left - rotate */
		if (key == GLUT_KEY_RIGHT) th += 5;
		else if (key == GLUT_KEY_LEFT) th -= 5;
		/*  Up/Down - elevation */
		else if (key == GLUT_KEY_UP) ph += 5;
		else if (key == GLUT_KEY_DOWN) ph -= 5;

	}

	/* Otherwise, just shift the screen */
	else {
		/*  Shift */
		if (key == GLUT_KEY_RIGHT) ecX += .5;
		else if (key == GLUT_KEY_LEFT) ecX -= .5;
		else if (key == GLUT_KEY_DOWN) ecY += .5;
		else if (key == GLUT_KEY_UP) ecY -= .5;
	}

	/*  Keep angles to +/-360 degrees */
	th %= 360;
	ph %= 360;
	redisplayAll();
}