#include "screencasts.h"
using namespace std;


//specify outputs, runs if the k key is pressed
void kimenet()
{
	bemenetek_szama = log2(szamlalo - 1);
	kimenetek_szama = szamok[0];
	for (int i = 1; i < szamlalo; i++)
	{
		kimenetek[szamok[0] - 1][i - 1] = szamok[i];

		//set back to initial value
		szamok[i] = 36;
	}
	szamok[0] = 36;
	szamlalo = 0;
}

// this runs if the r key is pressed. it will run a simulation step for the specified structure
// can also save to file, defined in input arguments
void futasv(char* file_name, bool first_open)
{
	ofstream fileki;
	ifstream filebe;
	int i, j, k, l, n;
	string sor;
	n = 50 / dt;
	int szamlal = 0;
	if (first_open) fileki.open("tmp.csv");
	else fileki.open("tmp.csv", ios::app);


	if (first_open) fileki << ",";
	for (i = 1; i <= 36; i++)
	{
		for (j = 1; j <= 36; j++)
		{
			for (k = 1; k <= 36; k++)
			{
				if (dronpa[i][j][k].van)
				{
					if (first_open) fileki << i - 1 << " " << j - 1 << " " << k - 1 << ",";
					szamlal++;
				}
			}
		}
	}
	if (first_open) fileki << endl;

	int *itomb = new int[szamlal];
	int *jtomb = new int[szamlal];
	int *ktomb = new int[szamlal];
	szamlal = 0;
	for (i = 1; i <= 36; i++)
	{
		for (j = 1; j <= 36; j++)
		{
			for (k = 1; k <= 36; k++)
			{
				if (dronpa[i][j][k].van)
				{
					itomb[szamlal] = i;
					jtomb[szamlal] = j;
					ktomb[szamlal] = k;
					szamlal++;
				}
			}
		}
	}


	for (l = 0; l < n; l++)
	{
		t += dt;
		fileki << t << ",";
		for (i = 0; i < szamlal; i++)
		{
			//case A
			if (l % 2 == 0)
			{
				U = K / tav*(dronpa[itomb[i] - 1][jtomb[i]][ktomb[i]].dipB + dronpa[itomb[i] + 1][jtomb[i]][ktomb[i]].dipB +
					dronpa[itomb[i]][jtomb[i] - 1][ktomb[i]].dipB + dronpa[itomb[i]][jtomb[i] + 1][ktomb[i]].dipB +
					dronpa[itomb[i]][jtomb[i]][ktomb[i] - 1].dipB + dronpa[itomb[i]][jtomb[i]][ktomb[i] + 1].dipB);

				if (dronpa[itomb[i]][jtomb[i]][ktomb[i]].ter) U += dronpa[itomb[i]][jtomb[i]][ktomb[i]].terMag;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB +
					(U / Ce1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB / (Ce1*Ce2))*dt;

				if (U > (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B;

				if (U < (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].dip = dronpa[itomb[i]][jtomb[i]][ktomb[i]].dipA = -100 +
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A;
			}

			//case B
			else
			{
				U = K / tav*(dronpa[itomb[i] - 1][jtomb[i]][ktomb[i]].dipA + dronpa[itomb[i] + 1][jtomb[i]][ktomb[i]].dipA +
					dronpa[itomb[i]][jtomb[i] - 1][ktomb[i]].dipA + dronpa[itomb[i]][jtomb[i] + 1][ktomb[i]].dipA +
					dronpa[itomb[i]][jtomb[i]][ktomb[i] - 1].dipA + dronpa[itomb[i]][jtomb[i]][ktomb[i] + 1].dipA);

				if (dronpa[itomb[i]][jtomb[i]][ktomb[i]].ter) U += dronpa[itomb[i]][jtomb[i]][ktomb[i]].terMag;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA +
					(U / Ce1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA / (Ce1*Ce2))*dt;

				if (U > (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A;

				if (U < (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].dip = dronpa[itomb[i]][jtomb[i]][ktomb[i]].dipB = -100 +
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B;
			}
			fileki << dronpa[itomb[i]][jtomb[i]][ktomb[i]].dip << ",";
		}
		fileki << endl;

	}
	fileki.close();

	filebe.open("tmp.csv");
	fileki.open(file_name);
	while (getline(filebe, sor))
	{
		for (i = 0; i < sor.length(); i++)
		{
			if (sor[i] == '.') fileki << ".";
			else fileki << sor[i];
		}
		fileki << endl;
	}
	filebe.close();
	fileki.close();

	delete[] itomb;
	delete[] jtomb;
	delete[] ktomb;
}

//same as the function before, but more fields can be applied at the same time
void futas()
{
	int i, j, k, l, n;
	string sor;
	n = 50 / dt;
	int szamlal = 0;



	for (i = 1; i <= 36; i++) {
		for (j = 1; j <= 36; j++) {
			for (k = 1; k <= 36; k++) {
				if (dronpa[i][j][k].van) {
					szamlal++;
				}
			}
		}
	}

	int *itomb = new int[szamlal];
	int *jtomb = new int[szamlal];
	int *ktomb = new int[szamlal];
	szamlal = 0;

	for (i = 1; i <= 36; i++) {
		for (j = 1; j <= 36; j++) {
			for (k = 1; k <= 36; k++) {
				if (dronpa[i][j][k].van) {
					itomb[szamlal] = i;
					jtomb[szamlal] = j;
					ktomb[szamlal] = k;
					szamlal++;
				}
			}
		}
	}

	for (l = 0; l < n; l++) {
		for (i = 0; i < szamlal; i++) {
			//case A
			if (l % 2 == 0) {
				U = K / tav*(dronpa[itomb[i] - 1][jtomb[i]][ktomb[i]].dipB + dronpa[itomb[i] + 1][jtomb[i]][ktomb[i]].dipB +
					dronpa[itomb[i]][jtomb[i] - 1][ktomb[i]].dipB + dronpa[itomb[i]][jtomb[i] + 1][ktomb[i]].dipB +
					dronpa[itomb[i]][jtomb[i]][ktomb[i] - 1].dipB + dronpa[itomb[i]][jtomb[i]][ktomb[i] + 1].dipB);

				if (dronpa[itomb[i]][jtomb[i]][ktomb[i]].ter) U += dronpa[itomb[i]][jtomb[i]][ktomb[i]].terMag;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB +
					(U / Ce1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB / (Ce1*Ce2))*dt;

				if (U > (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B;

				if (U < (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].dip = dronpa[itomb[i]][jtomb[i]][ktomb[i]].dipA = -100 +
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A;
			}

			//case B
			else {
				U = K / tav*(dronpa[itomb[i] - 1][jtomb[i]][ktomb[i]].dipA + dronpa[itomb[i] + 1][jtomb[i]][ktomb[i]].dipA +
					dronpa[itomb[i]][jtomb[i] - 1][ktomb[i]].dipA + dronpa[itomb[i]][jtomb[i] + 1][ktomb[i]].dipA +
					dronpa[itomb[i]][jtomb[i]][ktomb[i] - 1].dipA + dronpa[itomb[i]][jtomb[i]][ktomb[i] + 1].dipA);

				if (dronpa[itomb[i]][jtomb[i]][ktomb[i]].ter) U += dronpa[itomb[i]][jtomb[i]][ktomb[i]].terMag;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA +
					(U / Ce1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeA / (Ce1*Ce2))*dt;

				if (U > (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1A;

				if (U < (dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A / Cp2)) dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B =
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A + (U / Cp1 - dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A /
					(Cp1*Cp2))*dt;
				else dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B = dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2A;

				dronpa[itomb[i]][jtomb[i]][ktomb[i]].dip = dronpa[itomb[i]][jtomb[i]][ktomb[i]].dipB = -100 +
					dronpa[itomb[i]][jtomb[i]][ktomb[i]].qeB + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp1B + dronpa[itomb[i]][jtomb[i]][ktomb[i]].qp2B;
			}
		}
	}

	delete[] itomb;
	delete[] jtomb;
	delete[] ktomb;
}

//it calculates and sorts the neighbour graph numbers of a structure's molecules 
int *grafszam(int molekulaSzam)
{
	int *molekulak = new int[molekulaSzam];
	for (int i = 0; i < molekulaSzam; i++) {
		molekulak[i] = 0;
	}

	int *itomb = new int[molekulaSzam];
	int *jtomb = new int[molekulaSzam];
	int *ktomb = new int[molekulaSzam];
	int l = 0;
	int m = 0;

	for (int i = 18 - (molekulaSzam - 1); i <= 18 + molekulaSzam - 1; i++) {
		for (int j = 18 - (molekulaSzam - 1); j <= 18 + molekulaSzam - 1; j++) {
			for (int k = 18 - (molekulaSzam - 1); k <= 18 + molekulaSzam - 1; k++) {
				if (dronpa[i][j][k].van) {
					itomb[l] = i;
					jtomb[l] = j;
					ktomb[l] = k;
					l++;
				}
			}
		}
	}

	for (l = 0; l < molekulaSzam; l++) {
		for (m = 0; m < molekulaSzam; m++) {
			if (abs(itomb[m] - itomb[l]) + abs(jtomb[m] - jtomb[l]) + abs(ktomb[m] - ktomb[l]) == 1)
				molekulak[m]++;
		}
	}


	//bubble sort
	bool swapped = true;
	while (swapped) {
		swapped = false;
		for (int j = 0; j < molekulaSzam - 1; j++) {
			if (molekulak[j] < molekulak[j + 1]) {
				int k = molekulak[j];
				molekulak[j] = molekulak[j + 1];
				molekulak[j + 1] = k;

				swapped = true;
			}
		}
	}

	return molekulak;

	delete[] itomb;
	delete[] jtomb;
	delete[] ktomb;
	delete[] molekulak;
}

//compare be1 and be2
bool hasonlitas(int be1, int be2, int kacsacsor) {
	bool hasonlit = false;
	if ((be1) < be2 && kacsacsor == 0) hasonlit = true;

	if (be1>(be2) && kacsacsor == 1) hasonlit = true;

	if (kacsacsor == 2) hasonlit = true;

	return hasonlit;
}

//return the factorial of input number
int factorial(int f)
{
	if (f == 0) return 1;
	return(f * factorial(f - 1));
}

//random number generator
double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

//checks all the rows of a logical function, returns whether they are good
bool logikai_hasonlitas(int molekulaszam) {
	bool jo = true;
	int kimenetek_iteralo = 0;

	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			if (jo && protein[j].kimenet) {
				if (!hasonlitas(protein[j].actual[i], protein[j].desired[kimenetek[kimenetek_iteralo][i]], kimenetek[kimenetek_iteralo][i])) jo = false;
				kimenetek_iteralo++;
			}
		}
		kimenetek_iteralo = 0;
	}

	return jo;
}

//calculates fitness function
double fitness_func(int molekulaszam) {
	double fitness = 0;
	int kimenetek_iteralo = 0;
	bool same_dipole = true;

	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			if (protein[j].kimenet) {
				if (protein[j].actual[i] > protein[j].init_dipole + 0.01 || protein[j].actual[i] < protein[j].init_dipole - 0.01) same_dipole = false;
				if (kimenetek[kimenetek_iteralo][i]) {
					if (protein[j].actual[i] < protein[j].desired[kimenetek[kimenetek_iteralo][i]] + OVER_FIT) {
						fitness += pow(protein[j].desired[kimenetek[kimenetek_iteralo][i]] + OVER_FIT - protein[j].actual[i], 2);
					}
					kimenetek_iteralo++;
				}
				else if (!kimenetek[kimenetek_iteralo][i]) {
					if (protein[j].actual[i] > protein[j].desired[kimenetek[kimenetek_iteralo][i]] - OVER_FIT) {
						fitness += pow(protein[j].desired[kimenetek[kimenetek_iteralo][i]] - OVER_FIT - protein[j].actual[i], 2);
					}
					kimenetek_iteralo++;
				}
			}
		}
		kimenetek_iteralo = 0;
	}

	//filter the case when an output has no neighbours 
	if (same_dipole) fitness = 10000;
	fitness = fitness / 100;
	fitness = 1 / (1 + log(1 + fitness));
	return fitness;
}

//runs some simulation steps, if mentes is true then it saves the results to .csv file 
void SIMULATION(double **ter_vektor, bool mentes, int molekulaszam) {
	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			protein[j].reset_dipole(protein[j].init_dipole);
		}

		//apply field
		for (int j = 0; j < molekulaszam; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(ter_vektor[protein[j].bemenet_szam][bemenetek[i][protein[j].bemenet_szam]]);
			}
		}
		std::string ok = "graf" + std::to_string(i) + ".csv";
		char* c = &ok[0];
		if (mentes) futasv(c, true);
		else futas();


		for (int j = 0; j < molekulaszam; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(0);
			}
		}
		if (mentes) futasv(c, false);
		else futas();

		//save the dipole moments
		for (int j = 0; j <molekulaszam; j++) {
			if (protein[j].kimenet) {
				protein[j].update_actual(i);
			}
		}
	}
}

void protein_definialas() {
	

	for (int k = 0; k < molekulaSzam; k++) {
		protein[k];
	}

	// FOR TESTING ============ //
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				protein[k + j * 3 + i * 9].initialize_molekula(18 + i, 18 + j, 18 + k, false, 2, false);
			}
		}
	}

	//first run only once, and save the initial dipole values
	futas();
	for (int i = 0; i < molekulaSzam; i++) {
		protein[i].set_init_dipole();
		protein[i].set_desired();
		protein[i].set_actual();
	}
	// ========================= //
}


//main searching function based on GA
bool harmony_search() {
	molekulaSzam = DEF_PROTEIN_NUMBER;

	//initialize random field, first index is the number of the input molecule, 
	//second index is to specify whether we are talking about a field representing a 0 or 1 logic bit
	double **inputTer = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { inputTer[i] = new double[2]; }


	//GA params
	int generation = 1;
	bool hasonlit = false;

	int n = ITER_NUMBER;						//number of children
	int stuff = 1;								//counter for while cycle
	int fori;									//for for cycles

	double x, y;								//these will be random numbers

	//initialize original candidates
	for (int i = 0; i < bemenetek_szama; i++) {
		for (int j = 0; j < 2; j++) {
			inputTer[i][j] = START_POINT;
		}
	}

	//best candidate	
	double **best_ter = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { best_ter[i] = new double[2]; }

	for (int i = 0; i < bemenetek_szama; i++) {
		for (int j = 0; j < 2; j++) {
			best_ter[i][j] = inputTer[i][j];
		}
	}

	//more GA params
	DNA population[DEF_POP_SIZE];
	double mutationRate = DEF_MUT_RATE;
	std::vector<DNA> matingPool;
	double last_best = 0;

	//search loop
	while (!hasonlit && n>0) {
		//delete previous stuff, except first run
		if (n == ITER_NUMBER) {
			for (int k = 0; k < molekulaSzam; k++) {
				protein[k];
				protein[k].initialize_molekula(18, 18, 18, false, 2, false);
			}
		}

		//structure definition (not used in for GA algo)
		//protein_definialas();

		double best_fitness = 0;
		for (int i = 0; i < DEF_POP_SIZE; i++) {
			population[i].calcFitness();
			if (population[i].fitness > best_fitness) best_fitness = population[i].fitness;
			//also check logikai_hasonlitas
			hasonlit = logikai_hasonlitas(population[i].mol_szam);
			/*if (last_best == best_fitness) {
			hasonlit = true;
			}*/
			if (hasonlit) {
				molekulaSzam = population[i].mol_szam;
				best_ter = population[i].getFields();
				i = DEF_POP_SIZE;
			}
		}
		if (!hasonlit) {
			//create the mutation pool
			for (int i = 0; i < DEF_POP_SIZE; i++) {
				int nn = int(population[i].fitness * DEF_MATING_POOL_COEFF);
				for (int j = 0; j < nn; j++) {
					matingPool.push_back(population[i]);
				}
			}

			for (int i = 0; i < DEF_POP_SIZE; i++) {
				int a = int(fRand(0, matingPool.size() - 0.0000000001));
				int b = int(fRand(0, matingPool.size() - 0.0000000001));

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

		//iteration
		generation++;
		n--;
	}

	/* Write the data out */
	if (hasonlit) {
		for (int i = 0; i < bemenetek_szama; i++) {
			for (int j = 0; j < 2; j++) {
				cout << "Input molecule "<<i  << ", magnitude of logic field "<<j <<" : "<< best_ter[i][j] << "   ";
			}
			cout << endl;
		}
	}

	cout << "number of simulations: " << DEF_POP_SIZE*(ITER_NUMBER - n) << endl;
	//cout << "molekulak szama a strukturaban: " << molekulaSzam << endl;


	//delete everything
	for (int i = 0; i <bemenetek_szama; i++) { delete[] best_ter[i]; }
	delete[] best_ter;
	for (int i = 0; i <bemenetek_szama; i++) { delete[] inputTer[i]; }
	delete[] inputTer;

	return hasonlit;
}

//this takes the longest time, inside it is the harmony_search function with the main searching algorithm
void fofuggveny()
{
	bool sikerult = false;

	while (!sikerult) {

		//random inputs/outputs (not used for GA)
		//protein_definialas();

		//returns a bool, whether the structure needed was found
		sikerult = harmony_search();

		//write out the data
		if (sikerult) {
			cout << "Number of molecules: " << molekulaSzam << endl;
			for (int j = 0; j < molekulaSzam; j++) {
				if (protein[j].ter) {
					cout << protein[j].bemenet_szam << "<-applied this field number to: "
						<< protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
				}
				if (protein[j].kimenet) {
					cout << "output: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << "    ";
					cout << "dipole: ";
					for (int r = 0; r < pow(2, bemenetek_szama); r++) {
						cout << protein[j].actual[r] << " ";
					}
					cout << endl;
				}

				if (!protein[j].ter && !protein[j].kimenet) {
					cout << "other molecules: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
				}
			}
		}
	}
}




//save a structure to file
void save()
{
	ofstream strukt;
	strukt.open("strukt.csv");
	int i, j, k;
	for (i = 1; i <= 36; i++)
	{
		for (j = 1; j <= 36; j++)
		{
			for (k = 1; k <= 36; k++)
			{
				if (dronpa[i][j][k].van)
				{
					strukt << ";" << i - 1 << ";" << j - 1 << ";" << k - 1 << ";";
					if (dronpa[i][j][k].ter)
					{
						strukt << 1 << ";";
					}
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

//load a structure from file
void load()
{
	ifstream filebe;
	string sor;
	filebe.open("strukt.csv");

	int szam = -1;
	int szamok[4] = { 0,0,0,0 };
	bool param = false;
	int szamlalo = 2;

	while (getline(filebe, sor))
	{
		int vesszo = 0;
		int hely = 0;
		int sorhely = 0;
		if (param)
		{
			int az = 0;
			int ss = 2;
			if (sor.length() == 2) ss = 1;
			for (int i = 0; i < sor.length(); i++)
			{
				switch (sor.length())
				{
				case 1: az = sor[i] - 48;
					break;
				case 2: az += (sor[i] - 48)*pow(10, ss);
					break;
				case 3: az += (sor[i] - 48)*pow(10, ss);
					break;
				}
				ss--;
			}
			switch (szamlalo)
			{
			case 2: //th = az;
				break;
			case 1: //ph = az;
				break;
			case 0: //dim = az;
				break;
			}

			szamlalo--;
		}


		if (sor != "." && !param)
		{
			for (int i = 0; i < sor.length(); i++)
			{
				if (sor[i] == ';')
				{
					if (i - vesszo == 2)
					{
						szamok[sorhely] = (int)sor[i - 1] - 48;
						sorhely++;
					}
					if (i - vesszo == 3)
					{
						szamok[sorhely] = ((int)sor[i - 2] - 48) * 10 + (int)sor[i - 1] - 48;
						sorhely++;
					}
					vesszo = i;
				}

			}


			dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].van = true;
			if (szamok[3] == 1) dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].ter = true;

			szam = szamok[0] * 1296 + szamok[1] * 36 + szamok[2];

			for (int k = -18; k < 18; k++)
			{
				for (int j = -18; j < 18; j++)
				{
					for (int i = -18; i < 18; i++)
					{
						//draw the cubes
						if (szam == hely)
						{
							cubes[szam] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 } ,{ (float)szamok[0],(float)szamok[1],(float)szamok[2] } } };

						}
						hely++;
					}
				}
			}
		}
		if (sor == "." && !param) param = true;
	}
	filebe.close();

	filebe.open("params.dat");
	sor = "";
	szamlalo = 0;

	std::string::size_type sz;   // alias of size_t


	while (getline(filebe, sor))
	{
		switch (szamlalo)
		{
		case 0: th = stoi(sor, &sz);
			break;
		case 1: ph = stoi(sor, &sz);
			break;
		case 2: dim = stoi(sor, &sz);
			break;
		}
		szamlalo++;
	}
	filebe.close();
}



/*  HANDLE KEYBOARD AND MOUSE INPUTS  */

// at the press of q,w,e handle movement and rotation
void windowPmotion(int x, int y) {
	mouseX = x;
	mouseY = y;

	//rotate the light
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

	//rotate structure
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

	//move structure
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

//zoom
void mouseWheel(int scroll, int dir, int x, int y) {
	dim -= (double)dir;
	redisplayAll();
}

//handle mouse button presses
void windowMouse(int btn, int state, int x, int y) {
	if (btn == GLUT_LEFT_BUTTON) mouseBtnPressed = "Left";
	else if (btn == GLUT_RIGHT_BUTTON) mouseBtnPressed = "Right";

	if (state == GLUT_UP)
	{
		mouseState = "up";
	}

	if (state == GLUT_DOWN)
	{
		mouseState = "down";
	}

	redisplayAll();
}

//handle all keyboard presses (almost)
void windowKey(unsigned char key, int x, int y) {
	/*  Exit on ESC */
	if (key == 27) exit(0);

	//space
	if (key == 32) {
		if (valto == "parameters") valto = "coordinates";
		else valto = "parameters";
	}

	//define parameters
	if (valto == "parameters") {
		//toggle axes and parameter display
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


		else if (key == 't') enter = "field";

		//movement, rotation, light rotation
		else if (key == 'q') Shift = "movement";
		else if (key == 'w') Shift = "rotation";
		else if (key == 'e') Shift = "light";

		//set field magnitude
		else if (key == 'i') {
			//enter = "field magnitude";  ternagysag();
		}

		//set outputs
		else if (key == 'k') {
			enter = "output"; kimenet();
		}

		//run simulation
		else if (key == 'r') futasv("graf.csv", false);

		//run search algorithm
		else if (key == 'f') fofuggveny();


		//save structure
		else if (key == 's') save();
		//load structure
		else if (key == 'o') load();
	}


	/*  Spacebar */
	//else if (key == 32) toggleAnimation = 1 - toggleAnimation;
	/*  Change field of view angle */
	if (key == '-' && key>1) fov--;
	if (key == '+' && key<179) fov++;

	/*  Light elevation */
	if (key == '[') lightY -= 0.5;
	if (key == ']') lightY += 0.5;

	/*  Specular level */
	//else if (key == 's' && specular>0) specular -= 5;
	//else if (key == 'S' && specular<100) specular += 5;
	/*  Emission level */
	//else if (key == 'e' && emission>0) emission -= 5;
	//else if (key == 'E' && emission<100) emission += 5;
	/*  Shininess level */
	//else if (key == 'n' && shininess>-1) shininess -= 1;
	//else if (key == 'N' && shininess<7) shininess += 1;

	//enter and delete keys
	if (key == 13)	enter = "pressed";
	if (key == 8)  enter = "delete";

	//give alphanumeric coordinates
	if (valto == "coordinates") {
		for (int i = 97; i < 123; i++) {
			if (key == i) { szamok[szamlalo] = i - 87; szamlalo++; key = ','; }
		}

		for (int i = 48; i < 58; i++) {
			if (key == i) { szamok[szamlalo] = i - 48; szamlalo++; key = ','; }
		}
	}

	//coordinates from 0 to 9
	for (int i = 48; i < 58; i++) {
		if (key == i) { szamok[szamlalo] = i - 48; szamlalo++; key = ','; }
	}

	/*  Translate shininess power to value (-1 => 0) */
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
	/*  If holding shift, then rotate/elevate */
	if (modifiers == GLUT_ACTIVE_SHIFT) {
		/*  Right/Left - rotate */
		if (key == GLUT_KEY_RIGHT) th += 5;
		else if (key == GLUT_KEY_LEFT) th -= 5;
		/*  Up/Down - elevation */
		else if (key == GLUT_KEY_UP) ph += 5;
		else if (key == GLUT_KEY_DOWN) ph -= 5;

	}
	/*  Otherwise, just shift the screen */
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