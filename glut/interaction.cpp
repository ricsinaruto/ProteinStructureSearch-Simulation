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

//fitness function számoló
double fitness_func(int molekulaszam) {
	double fitness = 0;
	int kimenetek_iteralo = 0;

	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			if (protein[j].kimenet) {
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
	fitness = fitness / 1000;
	fitness = 1 / (1 + log(1 + fitness));
	return fitness;
}

//szimuláció sorozatot lefuttat, ha mentes igaz, akkor elmenti .csv-be
void SIMULATION(double **ter_vektor, bool mentes, int molekulaszam) {
	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < molekulaszam; j++) {
			protein[j].reset_dipole(protein[j].init_dipole);
		}

		//tér aplikálás
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

		//dipól értékek elmentése
		for (int j = 0; j <molekulaszam; j++) {
			if (protein[j].kimenet) {
				protein[j].update_actual(i);
			}
		}
	}
}

void protein_definialas() {
	///*
	bool ter_lok = false;
	bool kimenet_lok = false;
	bool redo = false;
	bool redo_structure = true;
	std::vector<int> placeholder_koordok;


	//előzőeket törölni, kivéve első futás
	if (molekulaSzam <= DEF_PROTEIN_NUMBER) {
		for (int k = 0; k < molekulaSzam; k++) {
			protein[k].delete_molekula();
		}
	}

	//random mennyiségű molekula
	molekulaSzam = fRand(3, DEF_PROTEIN_NUMBER + 0.9999999999999);
	cout << endl;
	cout << "molekulak szama: " << molekulaSzam << endl;

	for (int k = 0; k < molekulaSzam; k++) {
		protein[k];
	}

	//struktúra felépítés
	int elemek = 0;
	for (elemek = 0; elemek < molekulaSzam; elemek++) {
		if (elemek == 0) {
			protein[elemek].initialize_molekula(18, 18, 18, false, 2, false);
			protein[elemek].set_szomszedok();
		}

		else {
			std::vector<std::vector<int>> lista;
			std::set<std::vector<int>> elem_koord;
			for (int iter = 0; iter < elemek; iter++) {
				for (int i = 0; i < 6; i++) {
					lista.push_back(protein[iter].szomszedok[i]);
				}
				placeholder_koordok.push_back(protein[iter].x);
				placeholder_koordok.push_back(protein[iter].y);
				placeholder_koordok.push_back(protein[iter].z);
				elem_koord.insert(placeholder_koordok);
				for (int j = 0; j < 3; j++) placeholder_koordok.pop_back();
			}
			std::set<std::vector<int>> koord_set(lista.begin(), lista.end());
			lista.clear();

			//csak a különbözőket megtartani
			std::set<std::vector<int>> vegso_koordok;
			std::set_difference(koord_set.begin(), koord_set.end(), elem_koord.begin(), elem_koord.end(), std::inserter(vegso_koordok, vegso_koordok.end()));

			//random elem kiválasztás
			int k = fRand(0, vegso_koordok.size() - 0.0000000001);
			int i = 0;
			for (auto f : vegso_koordok) {
				if (i == k) {
					placeholder_koordok = f;
				}
				i++;
			}
			protein[elemek].initialize_molekula(placeholder_koordok[0], placeholder_koordok[1], placeholder_koordok[2], false, 2, false);
			protein[elemek].set_szomszedok();

			placeholder_koordok.clear();
			elem_koord.clear();
			koord_set.clear();
			vegso_koordok.clear();
		}
	}


	//1 darab kimenet
	int kimen = fRand(0, molekulaSzam - 0.000000000001);
	protein[kimen].kimenet = true;

	//random bemenetek
	for (int i = 0; i < molekulaSzam; i++) {
		int ter_legyen = fRand(0, 1.9999999999);
		if (ter_legyen) {
			int melyik_ter = fRand(0, bemenetek_szama - 0.00000000001);
			protein[i].ter = true;
			protein[i].set_ter_mol();
			protein[i].bemenet_szam = melyik_ter;
		}
	}



	// FOR TESTING ============ //
	/*
	molekulaSzam = 2;
	protein[0].initialize_molekula(18, 18, 18, false, 2, false);
	protein[1].initialize_molekula(18, 18, 19, false, 2, false);
	//1 darab kimenet
	protein[1].kimenet = true;
	protein[0].ter = true;
	protein[0].set_ter_mol();
	protein[0].bemenet_szam = 0;
	protein[1].ter = true;
	protein[1].set_ter_mol();
	protein[1].bemenet_szam = 1;
	// ========================= //

	/* ELŐZŐ ALGORITMUS */
	//manuális definiálás
	/*for (int i = 0; i < 3; i++) {
	for (int j = 0; j < 3; j++) {
	for (int k = 0; k < 3; k++) {
	protein[k + j * 3 + i * 9].initialize_molekula(18 + i, 18 + j, 18 + k, false, 2, false);
	}
	}
	}

	//töröljünk egyet vagy kettőt random
	int torles = fRand(0, 3 - 0.00000000001);
	if (torles > 0) {
	int egyik = fRand(0, molekulaSzam - 0.00000000001);
	protein[egyik].delete_molekula();
	if (torles > 1) {
	int masik = fRand(0, molekulaSzam - 0.00000000001);
	protein[masik].delete_molekula();
	while (egyik == masik) {
	masik = fRand(0, molekulaSzam - 0.00000000001);
	protein[masik].delete_molekula();
	}
	}
	}
	molekulaSzam -= torles;
	int shift = 0;
	for (int i = 0; i < molekulaSzam; i++) {
	if (protein[i].torolt) {
	shift++;
	}
	if (!protein[i+shift].torolt) {
	protein[i] = protein[i + shift];
	}
	else protein[i] = protein[i + shift + 1];

	}*/
}


//tér keresés
bool harmony_search() {

	//az első futást csak egyszer kell, és elmentjük az alap dipól értékeket
	futas();
	for (int i = 0; i < molekulaSzam; i++) {
		protein[i].set_init_dipole();
		protein[i].set_desired();
		protein[i].set_actual();
	}

	//random tér inicializálás, első index az input molekula száma, második index, hogy a 0 logikai értékű térről, vagy az 1 logikai értékű térről van-e szó
	double **inputTer = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { inputTer[i] = new double[2]; }



	//simulated annealing paraméterek, ahol [2][2] van azt majd át kell rakni dinamikusra
	int generation = 1;
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



	double nu = DEF_NU;							//gaussian nu-je


	double distro;								//amibe elmentjük a gaussian által létrehozott számot
	double z;									//a distrohoz kell
	double fitness;								//fitness számoláshoz
	double sugar = DEF_SUGAR;					//sugár a random generátorhoz
	double besto = 0;							//best fitness számoláshoz
	double bestoszam = 0;						//best fitness számoláshoz
	double finalbest = 0;						//best fitness számoláshoz

	DNA population[DEF_POP_SIZE];
	double mutationRate = DEF_MUT_RATE;
	std::vector<DNA> matingPool;

	//keresés
	while (!hasonlit && n>0) {
		double best_fitness = 0;
		for (int i = 0; i < DEF_POP_SIZE; i++) {
			population[i].calcFitness();
			if (population[i].fitness > best_fitness) best_fitness = population[i].fitness;
			//also check logikai_hasonlitas
			hasonlit = logikai_hasonlitas(molekulaSzam);
			if (hasonlit) {
				best_ter = population[i].getFields();
				i = DEF_POP_SIZE;
			}
		}

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


		cout << best_fitness << endl;;
		matingPool.clear();
		//iterálók
		generation++;
		n--;
	}

	/* Adatok kiíratása */
	if (hasonlit) {
		for (int i = 0; i < bemenetek_szama; i++) {
			for (int j = 0; j < 2; j++) {
				cout << i << ". input molekulara " << j << " logikai ter nagysaga: " << best_ter[i][j] << "   ";
			}
			cout << endl;
		}
	}

	//megadja a próbálgatások számát
	cout << "number of simulations: " << DEF_POP_SIZE*(ITER_NUMBER - n) << endl;
	//cout << "molekulak szama a strukturaban: " << molekulaSzam << endl;


	//nullázó, hogy többször lehessen futtatni a keresést anélkül hogy újraindítnánk a programot
	for (int i = 0; i < molekulaSzam; i++) {
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
void fofuggveny()
{
	bool sikerult = false;

	while (!sikerult) {

		//random input/outputok
		protein_definialas();

		sikerult = harmony_search();

		//adatok kiíratása
		if (sikerult) {
			for (int j = 0; j < molekulaSzam; j++) {
				if (protein[j].ter) {
					cout << protein[j].bemenet_szam << "-dik terrel terhelt: "
						<< protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
				}
				if (protein[j].kimenet) {
					cout << "kimenet: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << "    ";
					cout << "dipol: ";
					for (int r = 0; r < pow(2, bemenetek_szama); r++) {
						cout << protein[j].actual[r] << " ";
					}
					cout << endl;
				}

				if (!protein[j].ter && !protein[j].kimenet) {
					cout << "tobbi molekula: " << protein[j].x << " " << protein[j].y << " " << protein[j].z << endl;
				}
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