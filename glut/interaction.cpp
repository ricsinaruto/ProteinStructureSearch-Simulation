#include "screencasts.h"
using namespace std;


//kimenetek megadása, k lenyomásakor fut le
void kimenet()
{
	bemenetek_szama = log2(szamlalo-1);
	kimenetek_szama = szamok[0];
	for (int i = 1; i < szamlalo; i++)
	{
		kimenetek[szamok[0]-1][i-1] = szamok[i];

		//alaphelyzetre állítás
		szamok[i] = 36;
	}
	szamok[0] = 36;
	szamlalo = 0;
}

//r billentyû lenyomásakor ez fut le. A megadott struktúra szimulációját futtatja le
void futasv(char* file_name,bool first_open)
{
	ofstream fileki;
	ifstream filebe;
	int i, j, k, l, n;
	string sor;
	n = 50/dt;
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
	n = 50/dt;
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

//kiszámolja és sorba rakja egy struktúra molekuláinak a szomszédsági gráfszámát
int *grafszam(int molekulaSzam)
{
	int *molekulak = new int[molekulaSzam];
	for (int i = 0; i < molekulaSzam; i++)
	{
		molekulak[i] = 0;
	}

	int *itomb = new int[molekulaSzam];
	int *jtomb = new int[molekulaSzam];
	int *ktomb = new int[molekulaSzam];
	int l = 0;
	int m = 0;

	for (int i = 18 - (molekulaSzam - 1); i <= 18 + molekulaSzam - 1; i++)
	{
		for (int j = 18 - (molekulaSzam - 1); j <= 18 + molekulaSzam - 1; j++)
		{
			for (int k = 18 - (molekulaSzam - 1); k <= 18 + molekulaSzam - 1; k++)
			{
				if (dronpa[i][j][k].van)
				{
					itomb[l] = i;
					jtomb[l] = j;
					ktomb[l] = k;
					l++;
				}
			}
		}
	}

	for (l = 0; l < molekulaSzam; l++)
	{
		for (m = 0; m < molekulaSzam; m++)
		{
			if (abs(itomb[m] - itomb[l]) + abs(jtomb[m] - jtomb[l]) + abs(ktomb[m] - ktomb[l]) == 1)
				molekulak[m]++;

		}
	}


	bool swapped = true;
	while (swapped)
	{
		swapped = false;
		for (int j = 0; j < molekulaSzam - 1; j++)
		{
			if (molekulak[j] < molekulak[j + 1])
			{
				int k = molekulak[j];
				molekulak[j] = molekulak[j + 1];
				molekulak[j + 1] = k;

				swapped = true;
			}
		}
	}

	//for (int i=0;i<molekulaSzam;i++) cout << molekulak[i]<<endl;
	return molekulak;

	delete[] itomb;
	delete[] jtomb;
	delete[] ktomb;
	delete[] molekulak;
}

//összehasonlítást végez be1 és be2 között, 0 bemenetre igazt ad ha be1<be2, 1 bemenetre igazt ad ha be1>be2
bool hasonlitas(int be1, int be2, int kacsacsor)
{
	bool hasonlit=false;
	if ((be1) < be2 && kacsacsor==0) hasonlit = true;
	
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
bool logikai_hasonlitas(int molekulaszam) {
	bool jo = true;
	int kimenetek_iteralo = 0;

	for (int i = 0; i < pow(2,bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			if (jo && protein[j].kimenet) {
				if (!hasonlitas(protein[j].actual[i], protein[j].desired[kimenetek[kimenetek_iteralo][i]],kimenetek[kimenetek_iteralo][i])) jo = false;
				kimenetek_iteralo++;
			}
		}
		kimenetek_iteralo = 0;
	}

	return jo;
}

//fitness function számoló
double fitness_func(int molekulaszam) {
	double fitness=0;
	int kimenetek_iteralo = 0;

	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			if (protein[j].kimenet) {
				if (kimenetek[kimenetek_iteralo][i]) {
					if (protein[j].actual[i] < protein[j].desired[kimenetek[kimenetek_iteralo][i]]+OVER_FIT) {
						fitness += sqrt(pow(protein[j].desired[kimenetek[kimenetek_iteralo][i]]+OVER_FIT - protein[j].actual[i], 2));
						kimenetek_iteralo++;
					}
				}
				else {
					if (protein[j].actual[i] > protein[j].desired[kimenetek[kimenetek_iteralo][i]] -OVER_FIT) {
						fitness += sqrt(pow(protein[j].desired[kimenetek[kimenetek_iteralo][i]] -OVER_FIT - protein[j].actual[i], 2));
						kimenetek_iteralo++;
					}
				}
			}
		}
		kimenetek_iteralo = 0;
	}

	return fitness;
}

//szimuláció sorozatot lefuttat, ha mentes igaz, akkor elmenti .csv-be
void SIMULATION(double **ter_vektor,bool mentes, int molekulaszam) {
	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			protein[j].reset_dipole(protein[j].init_dipole);
		}

		//tér aplikálás
		int bemenet_iteralo = 0;
		for (int j = 0; j < molekulaszam; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(ter_vektor[bemenet_iteralo][bemenetek[i][bemenet_iteralo]]);
				bemenet_iteralo++;
				
			}
		}
		std::string ok = "graf" + std::to_string(i) + ".csv";
		char* c = &ok[0];
		if (mentes) futasv(c,true);
		else futas();

		bemenet_iteralo = 0;
		for (int j = 0; j < molekulaszam; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(0);
				bemenet_iteralo++;
			}
		}
		if (mentes) futasv(c,false);
		else futas();

		//dipól értékek elmentése
		for (int j = 0; j <molekulaszam; j++) {
			if (protein[j].kimenet) {
				protein[j].update_actual(i);
			}
		}
	}
}

//tér keresés
bool harmony_search(int molekulaszam) {

	//az első futást csak egyszer kell, és elmentjük az alap dipól értékeket
	futas();
	for (int i = 0; i < molekulaszam; i++) {
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
		for (fori=0; fori < n; fori++) {
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
			SIMULATION(child_ter,false,molekulaszam);

			//összehasonlítás, fitness
			fitness=fitness_func(molekulaszam);
			besto += fitness;
			

			//legyen-e csere?
			if (besto / (fori + (iteration - 1)*n) < (besto - fitness) / (fori + (iteration - 1)*n - 1)) {
				for (int i = 0; i < bemenetek_szama; i++) {
					for (int j = 0; j < 2; j++) {
						candidate_ter[i][j] = child_ter[i][j];
					}
				}
			}
			else if (fori>4 || iteration>1) {
				besto -= fitness;
				besto = besto + besto / (fori + (iteration - 1)*n);
			}
		}

		/* SIMULATION */ 
		SIMULATION(candidate_ter,false,molekulaszam);
		//összehasonlítás
		x = fRand(0, 1);
		fitness = fitness_func(molekulaszam);
		bestoszam += fitness;
		if (bestoszam / stuff < (bestoszam-fitness) / (stuff - 1)/* || 
			x < pow(e, ((1 / (bestoszam / stuff) - 1 / (bestoszam - fitness) / (stuff - 1))) / t)*/) {
			for (int i = 0; i < bemenetek_szama; i++) {
				for (int j = 0; j < 2; j++) {
					inputTer[i][j] = candidate_ter[i][j];
				}
			}
		}
		else if (iteration>2) {
			bestoszam -= fitness;
			bestoszam = bestoszam + bestoszam / stuff;
		}

		/* SIMULATION */
		SIMULATION(inputTer,MENTES,molekulaszam);
		//összehasonlítás 2
		fitness = fitness_func(molekulaszam);
		finalbest += fitness;
		if (finalbest / stuff < (finalbest - fitness) / (stuff - 1)) {
			for (int i = 0; i < bemenetek_szama; i++) {
				for (int j = 0; j < 2; j++) {
					best_ter[i][j] = inputTer[i][j];
				}
			}
		}
		//cout << "legjobb: " << finalbest/stuff<< "    current fitness: "<<fitness << endl;

		t--;
		iteration++;
		stuff++;

		//megfelelnek-e a logikai értékek
		hasonlit=logikai_hasonlitas(molekulaszam);
		
		//cout << hasonlit << endl;
	}

	/* Adatok kiíratása */
	for (int i = 0; i < bemenetek_szama; i++) {
		for (int j = 0; j < 2; j++) {
			cout <<i<<". input molekulara "<<j<<" logikai ter nagysaga: "<< best_ter[i][j] << "   ";
		}
		cout << endl;
	}
	//megadja a próbálgatások számát
	cout <<"number of simulations: "<< iteration*2+iteration*n << endl;
	cout << "molekulak szama a strukturaban: " << molekulaszam << endl;


	//nullázó, hogy többször lehessen futtatni a keresést anélkül hogy újraindítnánk a programot
	for (int i = 0; i < molekulaszam; i++) {
		protein[i].reset_dipole(DEF_DIPOL);
	}

	for (int i = 0; i <bemenetek_szama; i++) { delete[] candidate_ter[i]; }
	delete[] candidate_ter;
	for (int i = 0; i <bemenetek_szama; i++) { delete[] child_ter[i]; }
	delete[] child_ter;
	for (int i = 0; i <bemenetek_szama; i++) { delete[] best_ter[i]; }
	delete[] best_ter;
	for (int i = 0; i <bemenetek_szama; i++) { delete[] inputTer[i]; }
	delete[] inputTer;

	return hasonlit;
}

//ez tart a legtöbb ideig, ebben van egy adott struktúrán belül a különbözõ terekkel való tesztelés
bool terteszt(int molekulaSzam, int bemenetek_szam, int kimenetek_szam)
{

	int *itomb = new int[molekulaSzam];
	int *jtomb = new int[molekulaSzam];
	int *ktomb = new int[molekulaSzam];
	int l = 0;
	int m = 0;
	int n = 0;
	int p = 0;
	bool sikerult = false;

	for (int i = 18 - (molekulaSzam - 1); i <= 18 + molekulaSzam - 1; i++)
	{
		for (int j = 18 - (molekulaSzam - 1); j <= 18 + molekulaSzam - 1; j++)
		{
			for (int k = 18 - (molekulaSzam - 1); k <= 18 + molekulaSzam - 1; k++)
			{
				if (dronpa[i][j][k].van)
				{
					itomb[l] = i;
					jtomb[l] = j;
					ktomb[l] = k;

					protein[l].initialize_molekula(i, j, k, false, false, false);
					l++;
					
				}
			}
		}
	}

	int szimulacioszam = 0;
	int lehetosegek = 0;
	int tr1 = 0, tr2 = 0, tr3 = 0, tr4 = 0;
	int l0, m0, n0, p0, tr10, tr20, tr30, tr40;

	ofstream fileki;
	string fileki_string = "";
	string fileki_vegso = "";
	int ter_max = 10;
	if (bemenetek_szam == 4) { l0 = 0, m0 = 1, n0 = 2, p0 = 3; }
	else if (bemenetek_szam == 3) { l0 = molekulaSzam - 4, m0 = 0, n0 = 1, p0 = 2;  }
	else if (bemenetek_szam == 2) { l0 = molekulaSzam - 4, m0 = molekulaSzam - 3, n0 = 0, p0 = 1; }
	else if (bemenetek_szam == 1) { l0 = molekulaSzam - 4, m0 = molekulaSzam - 3, n0 = molekulaSzam - 2, p0 = 0; }

	for (l = l0; l < molekulaSzam - 3; l++) {
		for (m = m0 + l - l0; m < molekulaSzam - 2; m++) {
			for (n = n0 + m - m0; n < molekulaSzam - 1; n++) {
				for (p = p0 + n - n0; p < molekulaSzam; p++) {
					for (int i = 0; i < molekulaSzam; i++) {
						protein[i].kell = true;
						protein[i].kimenet = true;

						//teszt
						protein[p].kell = true;
						protein[p].ter = true;
						protein[n].kell = true;
						protein[n].ter = true;
						

						sikerult = harmony_search(molekulaSzam);

						//adatok kiíratása
						if (sikerult) {
							for (int j = 0; j < molekulaSzam; j++) {
								if (protein[j].ter) {
									cout << "terrel terhelt: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
								}
								if (protein[j].kimenet) {
									cout<<"kimenet: "<< protein[j].x << " " << protein[j].y << " " << protein[j].z << "    ";
									cout << "dipol: ";
									for (int r = 0; r < pow(2, bemenetek_szama); r++) {
										cout  << protein[j].actual[r]<<" ";
									}
									cout << endl;
								}

								if (!protein[j].kell) {
									cout << "tobbi molekula: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
								}
							}
							i = molekulaSzam;
							l = molekulaSzam;
							m = molekulaSzam;
							n = molekulaSzam;
							p = molekulaSzam;
							//cout << "siker" << endl;
						}


						protein[p].kell = false;
						protein[p].ter = false;
						protein[n].kell = false;
						protein[n].ter = false;

						protein[i].kell = false;
						protein[i].kimenet = false;

						
					}
				}
			}
		}
	}

	delete[] itomb;
	delete[] jtomb;
	delete[] ktomb;

	return sikerult;
}

//f-re lefutó szimulációs függvény
bool proba(int molekulaSzam)
{

	int eddigiMolekulak = molekulaSzam - 1;
	int i = 18, j = 18, k = 18;
	bool sikerult = false;


	int hasonlit_int = 1;
	int **osszehasonlito = new int*[1000];
	for (int i = 0; i < 1000; i++) { osszehasonlito[i] = new int[molekulaSzam]; }
	osszehasonlito[0] = grafszam(eddigiMolekulak);


	for (int i = 18 - (molekulaSzam - 1); i <= 18 + molekulaSzam - 1; i++)
	{
		for (int j = 18 - (molekulaSzam - 1); j <= 18 + molekulaSzam - 1; j++)
		{
			for (int k = 18 - (molekulaSzam - 1); k <= 18 + molekulaSzam - 1; k++)
			{
				if (!dronpa[i][j][k].van)
				{
					for (int l = 0; l < eddigiMolekulak; l++)
					{
						if (abs(i - itomb_mol[l]) + abs(j - jtomb_mol[l]) + abs(k - ktomb_mol[l]) == 1)
						{
							dronpa[i][j][k].van = true;
							dronpa[i][j][k].dip = -100;
							dronpa[i][j][k].dipA = -100;
							dronpa[i][j][k].dipB = -100;
							dronpa[i][j][k].ter = false;
						}
					}

					if (dronpa[i][j][k].van)
					{
						eddigiMolekulak = molekulaSzam;
						bool hasonlit = false;
						osszehasonlito[hasonlit_int] = grafszam(eddigiMolekulak);
						int p;
						for (p = 0; p < hasonlit_int;)
						{
							bool anyad = false;
							for (int h = 0; h < molekulaSzam; h++)
							{
								if (osszehasonlito[p][h] != osszehasonlito[hasonlit_int][h])
									anyad = true;

							}
							if (anyad) p++;
							else p = hasonlit_int + 1;
						}

						if (p == hasonlit_int)
						{
							//run
							if (bemenetek_szama + kimenetek_szama <= molekulaSzam)	sikerult = terteszt(eddigiMolekulak, bemenetek_szama, kimenetek_szama);
							if (sikerult)
							{
								i = 18 + molekulaSzam;
								j = 18 + molekulaSzam;
								k = 18 + molekulaSzam;
							}

							osszehasonlito[hasonlit_int] = grafszam(eddigiMolekulak);
							
							if (!sikerult) {
								cout << endl<<endl;
								cout << endl<<endl;
								cout << endl;
								cout << endl;
								cout << "egyszer" << endl;
							}
							hasonlit_int++;

							itomb_mol[struktura_szamlal] = i;
							jtomb_mol[struktura_szamlal] = j;
							ktomb_mol[struktura_szamlal] = k;
							struktura_szamlal++;
						}

						dronpa[i][j][k].van = false;
						dronpa[i][j][k].dip = 0;
						dronpa[i][j][k].dipA = 0;
						dronpa[i][j][k].dipB = 0;
						dronpa[i][j][k].ter = false;

						eddigiMolekulak = molekulaSzam - 1;
					}
				}
			}
		}
	}

	for (int i = 0; i < 1000; i++) { delete[] osszehasonlito[i]; }
	delete[] osszehasonlito;



	return sikerult;
}

//f-re lefutó szimulációs fõfüggvény
void fofuggveny()
{
	ofstream fileki;
	
	bool sikerult = false;
	int i = 18, j = 18, k = 18;
	dronpa[i - 1][j][k].van = true;
	dronpa[i - 1][j][k].dip = -100;
	dronpa[i - 1][j][k].dipA = -100;
	dronpa[i - 1][j][k].dipB = -100;
	dronpa[i - 1][j][k].ter = false;

	dronpa[i][j][k].van = true;
	dronpa[i][j][k].dip = -100;
	dronpa[i][j][k].dipA = -100;
	dronpa[i][j][k].dipB = -100;
	dronpa[i][j][k].ter = false;

	dronpa[i][j - 1][k].van = true;
	dronpa[i][j - 1][k].dip = -100;
	dronpa[i][j - 1][k].dipA = -100;
	dronpa[i][j - 1][k].dipB = -100;
	dronpa[i][j - 1][k].ter = false;

	itomb_mol[0] = i - 1;
	jtomb_mol[0] = j;
	ktomb_mol[0] = k;

	itomb_mol[1] = i;
	jtomb_mol[1] = j;
	ktomb_mol[1] = k;

	itomb_mol[2] = i;
	jtomb_mol[2] = j - 1;
	ktomb_mol[2] = k;


	if (bemenetek_szama + kimenetek_szama < 4)
	{
		struktura_szamlal = 3;
		int molekulaszam = 4;


		while (!sikerult)
		{
			if (molekulaszam == 4)
			{
				sikerult = proba(molekulaszam);
				molekulaszam++;
				struktura_szamlal = struktura_szamlal - 2;
			}
			else
			{
				int lul = struktura_szamlal;
				int i = 0;
				int j = 0;
				while (i<struktura_szamlal && j<molekulaszam - 1)
				{
					if (!dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].van)
					{
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].van = true;
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dip = -100;
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipA = -100;
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipB = -100;
						i++;
						j++;
					}
					else i++;
				}


				sikerult = proba(molekulaszam);
				cout << molekulaszam << endl;

				for (int i = 0; i<struktura_szamlal; i++)
				{
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].van = false;
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dip = 0;
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipA = 0;
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipB = 0;
				}


				molekulaszam++;
				//struktura_szamlal = molekulaszam;
			}

		}
	}

	else
	{
		dronpa[i][j][k - 1].van = true;
		dronpa[i][j][k - 1].dip = -100;
		dronpa[i][j][k - 1].dipA = -100;
		dronpa[i][j][k - 1].dipB = -100;
		dronpa[i][j][k - 1].ter = false;

		itomb_mol[3] = i;
		jtomb_mol[3] = j;
		ktomb_mol[3] = k - 1;

		struktura_szamlal = 4;
		int molekulaszam = 5;


		while (!sikerult)
		{
			if (molekulaszam == 5)
			{
				sikerult = proba(molekulaszam);
				molekulaszam++;
				struktura_szamlal = struktura_szamlal - 2;
			}
			else
			{
				int lul = struktura_szamlal;
				int i = 0;
				int j = 0;
				while (i<struktura_szamlal && j<molekulaszam - 1)
				{
					if (!dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].van)
					{
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].van = true;
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dip = -100;
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipA = -100;
						dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipB = -100;
						i++;
						j++;
					}
					else i++;
				}


				sikerult = proba(molekulaszam);
				cout << molekulaszam << endl;

				for (int i = 0; i<struktura_szamlal; i++)
				{
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].van = false;
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dip = 0;
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipA = 0;
					dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dipB = 0;
				}


				molekulaszam++;
				//struktura_szamlal = molekulaszam;
			}

		}
	}
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
					strukt << ";"<< i - 1 << ";" << j - 1 << ";" << k - 1 << ";";
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
	else if (Shift=="rotation") {
		if (mouseBtnPressed == "Left")	{
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
	else if (Shift=="movement") {
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
			ecX = ecX2 + (-mouseX + xcoord)/10;
			ecY = ecY2 + (mouseY - ycoord)/10;
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
		else if (key == 'r') futasv("graf.csv",false);

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
			if (key == i) { szamok[szamlalo] = i-87; szamlalo++; key = ','; }
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