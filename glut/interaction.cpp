#include "screencasts.h"
using namespace std;


//kimenetek megadása, k lenyomásakor fut le
void kimenet()
{
	bemenetek_szama = log2(szamlalo - 1);
	kimenetek_szama = szamok[0];
	for (int i = 1; i < szamlalo; i++)
	{
		kimenetek[szamok[0] - 1][i - 1] = szamok[i];

		//alaphelyzetre állítás
		szamok[i] = 36;
	}
	szamok[0] = 36;
	szamlalo = 0;
}

//r billentyû lenyomásakor ez fut le. A megadott struktúra szimulációját futtatja le
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
			//A eset
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

			//B eset
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

//ugyanaz mint a futas, csak egyszerre több tér is lehet
void futas()
{
	int i, j, k, l, n;
	string sor;
	n = 50 / dt;
	int szamlal = 0;



	for (i = 1; i <= 36; i++)
	{
		for (j = 1; j <= 36; j++)
		{
			for (k = 1; k <= 36; k++)
			{
				if (dronpa[i][j][k].van)
				{
					szamlal++;
				}
			}
		}
	}

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
		for (i = 0; i < szamlal; i++)
		{
			//A eset
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

			//B eset
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
		}
	}

	delete[] itomb;
	delete[] jtomb;
	delete[] ktomb;
}

//összehasonlítást végez be1 és be2 között, 0 bemenetre igazt ad ha be1<be2, 1 bemenetre igazt ad ha be1>be2
bool hasonlitas(int be1, int be2, int kacsacsor)
{
	bool hasonlit = false;
	if ((be1) < be2 && kacsacsor == 0) hasonlit = true;

	if (be1>(be2) && kacsacsor == 1) hasonlit = true;

	if (kacsacsor == 2) hasonlit = true;

	return hasonlit;
}

//faktoriális számítás
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

//logikai függvény összes sorát megnézi, hogy jó-e
bool logikai_hasonlitas() {
	bool jo = true;
	int kimenetek_iteralo = 0;

	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (jo && protein[j].kimenet) {
				if (!hasonlitas(protein[j].actual[i], protein[j].desired[kimenetek[kimenetek_iteralo][i]], kimenetek[kimenetek_iteralo][i])) jo = false;
				kimenetek_iteralo++;
			}
		}
		kimenetek_iteralo = 0;
	}

	return jo;
}

//fitness function számoló
double fitness_func() {
	double fitness = 0;
	int kimenetek_iteralo = 0;

	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (protein[j].kimenet) {
				if (kimenetek[kimenetek_iteralo][i]) {
					if (protein[j].actual[i] < protein[j].desired[kimenetek[kimenetek_iteralo][i]] + OVER_FIT) {
						fitness += sqrt(pow(protein[j].desired[kimenetek[kimenetek_iteralo][i]] + OVER_FIT - protein[j].actual[i], 2));
					}
					kimenetek_iteralo++;
				}
				else if (!kimenetek[kimenetek_iteralo][i]) {
					if (protein[j].actual[i] > protein[j].desired[kimenetek[kimenetek_iteralo][i]] - OVER_FIT) {
						fitness += sqrt(pow(protein[j].desired[kimenetek[kimenetek_iteralo][i]] - OVER_FIT - protein[j].actual[i], 2));
					}
					kimenetek_iteralo++;
				}
			}
		}
		kimenetek_iteralo = 0;
	}

	return fitness;
}

//szimuláció sorozatot lefuttat, ha mentes igaz, akkor elmenti .csv-be
void SIMULATION(double **ter_vektor, bool mentes) {
	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			protein[j].reset_dipole(protein[j].init_dipole);
		}

		//tér aplikálás
		int bemenet_iteralo = 0;
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(ter_vektor[bemenet_iteralo][bemenetek[i][bemenet_iteralo]]);
				bemenet_iteralo++;
			}
		}
		std::string ok = "graf" + std::to_string(i) + ".csv";
		char* c = &ok[0];
		if (mentes) futasv(c, true);
		else futas();

		bemenet_iteralo = 0;
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(0);
				bemenet_iteralo++;
			}
		}
		if (mentes) futasv(c, false);
		else futas();

		//dipól értékek elmentése
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (protein[j].kimenet) {
				protein[j].update_actual(i);
			}
		}
	}
}

//a fitnessek összehasonlításához függvényke
double *fitness_compare(double finalbest, double fitness, int stuff, int iteration) {
	finalbest += fitness;
	double finalbest1;
	int stuff1;
	if (stuff % 10 == 0) {
		finalbest = finalbest / 10; stuff1 = 1;
	}
	else stuff1 = 0;
	if (stuff > 1) finalbest1 = finalbest / ((stuff % 10) + 1);
	else if (iteration == 2) finalbest1 = finalbest / ((stuff % 10) + 1);
	else finalbest1 = finalbest;

	double fitnessbest;
	if (stuff % 10 != 0) fitnessbest = (finalbest - fitness) / ((stuff % 10) + stuff1);
	else fitnessbest = (finalbest * 10 - fitness) / (10);
	double fitnessfinal; 

	if (stuff % 10 == 0) {
		fitnessfinal = (finalbest * 10 - fitness + fitnessbest)/10;
	}
	else fitnessfinal = finalbest - fitness + fitnessbest;

	double zulul[4] = { finalbest1, fitnessbest,finalbest,fitnessfinal };
	return zulul;
}

//tér keresés
void harmony_search() {

	//az első futást csak egyszer kell, és elmentjük az alap dipól értékeket
	futas();
	for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
		protein[i].set_init_dipole();
		protein[i].set_desired();
		protein[i].set_actual();
	}

	//random tér inicializálás, első index az input molekula száma, második index, hogy a 0 logikai értékű térről, vagy az 1 logikai értékű térről van-e szó
	double **inputTer = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { inputTer[i] = new double[2]; }



	//simulated annealing paraméterek, ahol [2][2] van azt majd át kell rakni dinamikusra
	// /* SIMULATION */ algoritmust átrakni egy függvénybe
	int iteration = 1;
	bool hasonlit = false;


	int n = ITER_NUMBER;						//number of children

	int stuff = 1;								//a while számlálója
	int fori;									//for ciklusokhoz

	double x, y;								//ez lesz egy random szám

												//random tér inicializálás (original candidate)
	for (int i = 0; i < bemenetek_szama; i++) {
		for (int j = 0; j < 2; j++) {
			inputTer[i][j] = START_POINT;
		}
	}

	//egy lehetséges tér
	double **candidate_ter = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { candidate_ter[i] = new double[2]; }
	//új generáció
	double **child_ter = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { child_ter[i] = new double[2]; }
	//legjobb megoldás		
	double **best_ter = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { best_ter[i] = new double[2]; }

	for (int i = 0; i < bemenetek_szama; i++) {
		for (int j = 0; j < 2; j++) {
			best_ter[i][j] = inputTer[i][j];
		}
	}



	int t = DEF_TEMP;							//"temperature"
	double sigma = DEF_SIGMA;					//gaussian sigma-ja
	double nu = DEF_NU;							//gaussian nu-je


	double distro;								//amibe elmentjük a gaussian által létrehozott számot
	double z;									//a distrohoz kell
	double fitness;								//fitness számoláshoz
	double sugar = DEF_SUGAR;					//sugár a random generátorhoz
	double besto = 0;							//best fitness számoláshoz
	double bestoszam = 0;						//best fitness számoláshoz
	double finalbest = 0;						//best fitness számoláshoz
	double *compare=new double[4];				//ebbe tároljuk a fitness_compare return értékeit

												//keresés
	while (!hasonlit && t>0) {
		//random szám generálás 1
		for (int i = 0; i < bemenetek_szama; i++) {
			for (int j = 0; j < 2;) {
				z = 0;

				while (z <= 0 || z >= 1) {
					x = fRand(-sugar, sugar);
					y = fRand(-sugar, sugar);
					z = x*x + y*y;
				}
				//gaussian random szám
				distro = nu + x*sigma * sqrt(-2 * log(z) / z);

				candidate_ter[i][j] = inputTer[i][j] + distro;
				if (candidate_ter[i][j]>max_ter || candidate_ter[i][j]<-max_ter);
				else j++;
			}
		}

		//generáció szimulálása
		for (fori = 0; fori < n; fori++) {
			//random szám generálás 2
			for (int i = 0; i < bemenetek_szama; i++) {
				for (int j = 0; j < 2;) {
					z = 0;

					while (z <= 0 || z >= 1) {
						x = fRand(-sugar, sugar);
						y = fRand(-sugar, sugar);
						z = x*x + y*y;
					}
					//gaussian random szám
					distro = nu + x*sigma * sqrt(-2 * log(z) / z);

					child_ter[i][j] = inputTer[i][j] + distro;
					if (child_ter[i][j]>max_ter || child_ter[i][j]<-max_ter);
					else j++;
				}
			}

			/* SIMULATION */
			SIMULATION(child_ter, false);

			//összehasonlítás, fitness
			fitness = fitness_func();
			besto += fitness;

			//legyen-e csere?
			if (besto / (fori + (iteration - 1)*n) < (besto - fitness) / (fori + (iteration - 1)*n - 1)) {
				for (int i = 0; i < bemenetek_szama; i++) {
					for (int j = 0; j < 2; j++) {
						candidate_ter[i][j] = child_ter[i][j];
					}
				}
			}
			else if (fori>0) {
				besto -= fitness;
				besto = besto + besto / (fori + (iteration - 1)*n);
			}
		}

		/* SIMULATION */
		SIMULATION(candidate_ter, false);
		//összehasonlítás
		fitness = fitness_func();
		compare = fitness_compare(bestoszam, fitness, stuff, iteration);
		x = fRand(0, pow(e,-abs(1/(compare[0]-compare[1]))));
		bestoszam = compare[2];
		if (compare[0] <= compare[1] || (DEF_TEMP_BOOL &&
			x < pow(e, (1 / compare[0] - 1 / compare[1]) / (DEF_TEMP - t)*DEF_TEMP_CONST))) {
			for (int i = 0; i < bemenetek_szama; i++) {
				for (int j = 0; j < 2; j++) {
					inputTer[i][j] = candidate_ter[i][j];
				}
			}
		}
		else if (iteration>1) bestoszam = compare[3];
		//if (x < pow(e, (1 / compare[1] - 1 / compare[0]) / (t/1500))) cout << endl;
		/* SIMULATION */
		SIMULATION(inputTer, MENTES);
		//összehasonlítás 2
		fitness = fitness_func();
		compare = fitness_compare(finalbest,fitness,stuff,iteration);
		finalbest = compare[2];

		if ( compare[0]<compare[1] ) {
			for (int i = 0; i < bemenetek_szama; i++) {
				for (int j = 0; j < 2; j++) {
					best_ter[i][j] = inputTer[i][j];
				}
			}
		}
		else if (iteration > 1) finalbest = compare[3];
		cout << "legjobb: " << compare[0] << "    current fitness: " << fitness << endl;

		t--;
		iteration++;
		if (iteration>2) stuff++;

		//megfelelnek-e a logikai értékek
		hasonlit = logikai_hasonlitas();

		//cout << hasonlit << endl;
	}

	/* Adatok kiíratása */
	for (int i = 0; i < bemenetek_szama; i++) {
		for (int j = 0; j < 2; j++) {
			cout << i << ". input molekulara " << j << " logikai ter nagysaga: " << best_ter[i][j] << "   ";
		}
		cout << endl;
	}

	//adatok kiírása
	cout << "number of simulations: " << iteration * 2 + iteration*n << endl << endl;
	if (hasonlit) {
		cout << "sikerult" << endl;
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (protein[j].ter) {
				cout << "terrel terhelt: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
			}
			if (protein[j].kimenet) {
				cout << "kimenet: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << "    ";
				cout << "dipol: ";
				for (int r = 0; r < pow(2, bemenetek_szama); r++) {
					cout << protein[j].actual[r] << " ";
				}
				cout << endl;
			}

			if (!protein[j].kell) {
				cout << "tobbi molekula: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
			}
		}
	}

	//nullázó, hogy többször lehessen futtatni a keresést anélkül hogy újraindítnánk a programot
	for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
		protein[i].reset_dipole(DEF_DIPOL);
	}
}

//f-re lefutó struktúra kereső fõfüggvény
void fofuggveny()
{

	//XOR és XNOR struktúrát még nem talált
	protein[0].initialize_molekula(17, 18, 18, true, true, true);
	protein[1].initialize_molekula(19, 18, 18, true, true, true);
	protein[2].initialize_molekula(21, 18, 18, true, true, true);
	protein[3].initialize_molekula(23, 18, 18, true, true, true);

	protein[4].initialize_molekula(18, 18, 18, false, false, false);
	protein[5].initialize_molekula(20, 18, 18, false, false, false);
	protein[6].initialize_molekula(22, 18, 18, false, false, false);



	//tér keresés
	harmony_search();
}

//struktúra elmentése
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

//struktúra beolvasása
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
						//ez megrajzolja a kockát
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



/*  BILLENTYŰ ÉS EGÉR KEZELÉS  */

// q,w,e billentyûk lenyomásakor a mozgatások és forgások kezelése
void windowPmotion(int x, int y) {
	mouseX = x;
	mouseY = y;

	//fény forgatása
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

	//struktúra forgatása
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

	//struktúra mozgatása
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

//zoomolás
void mouseWheel(int scroll, int dir, int x, int y) {
	dim -= (double)dir;
	redisplayAll();
}

//egér gombok lenyomásának kezelése
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

//összes billentyû lenyomás kezelése
void windowKey(unsigned char key, int x, int y) {
	/*  Exit on ESC */
	if (key == 27) exit(0);

	//space
	if (key == 32) {
		if (valto == "parameters") valto = "coordinates";
		else valto = "parameters";
	}

	//paraméterek megadása
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

		//mozgás, forgás, fény forgás
		else if (key == 'q') Shift = "movement";
		else if (key == 'w') Shift = "rotation";
		else if (key == 'e') Shift = "light";

		//tér nagyságának megadása
		else if (key == 'i') {
			//enter = "field magnitude";  ternagysag();
		}

		//kimenet megadás
		else if (key == 'k') {
			enter = "output"; kimenet();
		}

		//szimuláció futtatása
		else if (key == 'r') futasv("graf.csv", false);

		//próba függvény
		else if (key == 'f') fofuggveny();


		//struktúra elmentése
		else if (key == 's') save();
		//struktúra betöltése
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

	//enter és delete gombok
	if (key == 13)	enter = "pressed";
	if (key == 8)  enter = "delete";

	//alfanumerikus koordináták megadása
	if (valto == "coordinates") {
		for (int i = 97; i < 123; i++) {
			if (key == i) { szamok[szamlalo] = i - 87; szamlalo++; key = ','; }
		}

		for (int i = 48; i < 58; i++) {
			if (key == i) { szamok[szamlalo] = i - 48; szamlalo++; key = ','; }
		}
	}

	//0-9ig koordináta megadás
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