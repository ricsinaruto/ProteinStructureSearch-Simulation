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
			if (protein[j].kimenet) {
				if (!hasonlitas(protein[j].actual[i], protein[j].desired[kimenetek[kimenetek_iteralo][i]], kimenetek[kimenetek_iteralo][i])) {
					jo = false;
					szamok[i] = 0;
				}
				else szamok[i] = 1;
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
		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(ter_vektor[protein[j].bemenet_szam][bemenetek[i][protein[j].bemenet_szam]]);
			}
		}
		std::string ok = "graf" + std::to_string(i) + ".csv";
		char* c = &ok[0];
		if (mentes) futasv(c, true);
		else futas();


		for (int j = 0; j < DEF_PROTEIN_NUMBER; j++) {
			if (protein[j].ter) {
				protein[j].set_ter(0);
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

//gráf és proteinek definiálása
void protein_definialas() {
	double *szomszedsag = new double[DEF_PROTEIN_NUMBER];
	
	bool ter_lok=false;
	bool kimenet_lok=false;
	bool redo = false;
	bool redo_structure = true;
	
	for (int k = 0; k < DEF_PROTEIN_NUMBER; k++) {
		protein[k];
	}

	/*while (redo_structure) {
		redo_structure = false;
		for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
			redo = false;
			for (int k = 0; k < DEF_PROTEIN_NUMBER; k++) {
				szomszedsag[k] = 0;
			}
			int elso_szamlal = 0;

			//szomszédok feltöltése
			for (int j = 1 + i; j < DEF_PROTEIN_NUMBER; j++) {
				if (0.5 > fRand(0, 1)) szomszedsag[j] = 0;
				else szomszedsag[j] = 1;
				if (i == 0 && szomszedsag[j]==1) elso_szamlal++;
			}
			if (i==0 && elso_szamlal == 0) redo = true;
			if (i == 0) protein[i].initialize_molekula(18, 18, 18, ter_lok, kimenet_lok, szomszedsag);
			else {
				int hany_szomszed = 0;
				//felső háromszög feltöltése
				for (int it = 0; it < i; it++) {
					szomszedsag[it] = protein[it].szomszedok[i];
					if (szomszedsag[it] == 1) hany_szomszed++;
				}

				//nézni, hogy szomszédokhoz próbálunk csatolni-e
				for (int szomszed_it = 1 + i; szomszed_it < DEF_PROTEIN_NUMBER; szomszed_it++) {
					if (szomszedsag[szomszed_it] == 1) {
						hany_szomszed++;
						for (int elozo = 0; elozo < i; elozo++) {
							if (protein[elozo].szomszedok[szomszed_it] == 1 && protein[elozo].szomszedok[szomszed_it - i] == 1) redo = true;
						}
					}
				}

				if (hany_szomszed == 0) redo = true; 
				if (hany_szomszed==0 && i==DEF_PROTEIN_NUMBER-1) { redo_structure = true; i = DEF_PROTEIN_NUMBER+2; }
			
				//nézni hogy mikor lesz 3 azonos érték, Ha felső háromszögre jön ki a 3 azonos szomszéd akkor redo structure
				if (hany_szomszed >= 3) {
					for (int elozo = 0; elozo < i; elozo++) {
						int szamlalo = 0;
						bool felso = true;
						for (int it = 0; it < DEF_PROTEIN_NUMBER; it++) {
							if (protein[elozo].szomszedok[it] - szomszedsag[it] == 0 && szomszedsag[it] == 1) {
								szamlalo++;
								if (it > i) felso = false;
							}
						}
						if (szamlalo >= 3) redo = true;
						if (szamlalo>=3 && felso) { redo_structure = true; i = DEF_PROTEIN_NUMBER + 2; }
					}
				}


			}

			if (redo) i--;
			//struktúra rakás
			else if(i!=0) {
				bool redo_this = true;
				int stuff = 0;
				while (redo_this && stuff<100) {
					redo_this = false;
					int x = fRand(18 - i, 18 + i+0.9999999);
					int y = fRand(18 - i, 18 + i+0.9999999);
					int z = fRand(18 - i, 18 + i+0.9999999);

					for (int it = 0; it < i; it++) {
						if (szomszedsag[it] == 1 && abs(protein[it].x - x)+ abs(protein[it].y - y) + abs(protein[it].z - z) != 1) redo_this = true;
						if (szomszedsag[it] == 0 && abs(protein[it].x - x) + abs(protein[it].y - y) + abs(protein[it].z - z) == 1) redo_this = true;
						if (abs(protein[it].x - x) + abs(protein[it].y - y) + abs(protein[it].z - z) == 0) redo_this = true;
					}
					if (redo_this == false) {
						protein[i].initialize_molekula(x, y, z, ter_lok, kimenet_lok, szomszedsag);
					}
					stuff++;
					//cout << stuff << endl;
				}
				if (stuff >= 100) { redo_structure = true; i = DEF_PROTEIN_NUMBER; }
			}
			//DEBUGGING, 4 molekulára még nem mükszik
			//cout << i << endl;
		}
	}*/

	//manuális definiálás
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				protein[k + j * 3 + i * 9].initialize_molekula(18 + i, 18 + j, 18 + k, false, 2, false);
			}
		}
	}

	//1 darab kimenet
	int kimen = fRand(0, DEF_PROTEIN_NUMBER - 0.000000000001);
	protein[kimen].kimenet = true;

	for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
		int ter_legyen= fRand(0, 1.9999999999);
		if (ter_legyen) {
			int melyik_ter= fRand(0, bemenetek_szama-0.00000000001);
			protein[i].ter = true;
			protein[i].set_ter_mol();
			protein[i].bemenet_szam = melyik_ter;
		}
	}




	//két bemenetre, 1 kimenetre
	/*bool ter_figyel = true;
	int egyik, masik;
	while (ter_figyel) {
		ter_figyel = false;
		egyik = fRand(0, DEF_PROTEIN_NUMBER - 0.000000000001);
		masik = fRand(0, DEF_PROTEIN_NUMBER - 0.000000000001);
		if (egyik == masik) ter_figyel = true;
	}
	
	int kimen = fRand(0, DEF_PROTEIN_NUMBER - 0.000000000001);
	protein[egyik].ter = true; protein[egyik].set_ter_mol();
	protein[masik].ter = true; protein[masik].set_ter_mol();
	protein[kimen].kimenet = true;

	for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
		if (i != egyik && i != masik) {
			protein[i].ter = false;
		}
	}*/

	//initial dipol futás
	futas();
	for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
		protein[i].set_init_dipole();
		protein[i].set_desired();
		protein[i].set_actual();
	}

}

//tér keresés
void harmony_search() {
	ofstream weka_file,fileki,fileki_input,fileki_bin_input,fileki_bool,fileki_dipol;
	fileki.open("output_data.dat");
	fileki_bool.open("output_data_4_bool.dat");
	fileki_dipol.open("output_data_dipol.dat");
	fileki_input.open("input_data.dat");
	fileki_bin_input.open("bin_input_data.dat");
	weka_file.open("weka_file.dat");

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


	int distro;								//amibe elmentjük a gaussian által létrehozott számot
	double z;									//a distrohoz kell
	double fitness;								//fitness számoláshoz
	double sugar = DEF_SUGAR;					//sugár a random generátorhoz
	double besto = 0;							//best fitness számoláshoz
	double bestoszam = 0;						//best fitness számoláshoz
	double finalbest = 0;						//best fitness számoláshoz
	double *compare=new double[4];				//ebbe tároljuk a fitness_compare return értékeit
	bool jo = true;
	int megjegyez = 0;
	int logikai_fv = 0;
	int ter_megj = 0;
	int egyesek = 0;
	int nullasok = 0;

	/* KERESÉS */
	while (t==DEF_TEMP) {
		//több dolog is optimalizálva van 3 molekulára jelenleg

		logikai_fv = 0;
		/*for (int i = 0; i < bemenetek_szama; i++) {
			for (int j = 0; j < 2;j++) {
				 child_ter[i][j] = fRand(-max_ter,max_ter);

				 distro = fRand(0, 1.99999999999999999);
				 kimenetek[0][2*i+j] = distro;
			}
		}*/
		

		//protein tér, struktúra definiálás
		//protein_definialas();
		

		/*ELŐZŐ STRUKTÚRÁS MEGOLDÁS*/
		//if (protein[0].szomszedok[1] == protein[0].szomszedok[2]) fileki_input << 1 << " ";
		//else if (protein[0].szomszedok[1] ==1) fileki_input << 2 << " ";
		//else fileki_input << 3 << " ";
		
		/*megjegyez = 0;
		ter_megj = 0;
		for (int szamlal = 0; szamlal < DEF_PROTEIN_NUMBER; szamlal++) {
			if (protein[szamlal].ter) ter_megj += (szamlal + 1);
			if (protein[szamlal].kimenet) {
				megjegyez = szamlal;
			}
		}
		fileki_input << ter_megj-2 << " ";
		fileki_input << megjegyez + 1 << " ";*/
		

		//tesztelés
		child_ter[0][0] = 100;
		child_ter[0][1] = 75.013;
		child_ter[1][0] = -7.777;
		child_ter[1][1] = -0.681;

		/* SIMULATION */
		SIMULATION(child_ter, false);
		jo = logikai_hasonlitas();
		if (jo) egyesek++;
		else nullasok++;

		//kiíratás
		if (nullasok < DEF_TEMP || jo) {
			logikai_fv = 0;
			for (int i = 0; i < bemenetek_szama; i++) {
				for (int j = 0; j < 2; j++) {
					fileki_input << child_ter[i][j] << " ";
					fileki_bin_input << child_ter[i][j] << " ";
					weka_file << child_ter[i][j] << ",";

					logikai_fv += kimenetek[0][2 * i + j] *pow(2, 2 * i + j);
					fileki_bin_input << kimenetek[0][2 * i + j] << " ";
					weka_file << kimenetek[0][2 * i + j] << ",";
				}
			}
			fileki_input << logikai_fv << " ";
			weka_file<< logikai_fv << ",";

			megjegyez = 0;
			for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
				fileki_input << protein[i].bemenet_szam + 1 << " ";
				fileki_bin_input << protein[i].bemenet_szam << " ";
				weka_file << protein[i].bemenet_szam+1 << ",";

				if (protein[i].kimenet) megjegyez = i;
			}
			fileki_input << megjegyez + 1 << " ";
			fileki_bin_input << megjegyez << " ";
			weka_file << megjegyez+1 << ",";

			
			
			for (int i = 0; i < pow(2, bemenetek_szama); i++) {
				fileki_dipol << protein[megjegyez].actual[i] << " ";
				cout << protein[megjegyez].actual[i] << " ";
				weka_file << protein[megjegyez].actual[i] << ",";
				fileki_bool << szamok[i] << " ";
				cout << szamok[i] << " ";
				weka_file << szamok[i] << ",";
			}
			weka_file << jo;
		
			weka_file << endl;
			fileki << endl;
			fileki_input << endl;
			fileki_bin_input << endl;
			fileki_bool << endl;
			fileki_dipol << endl;
		}
		

		//cout << jo<<endl;
		
		t--;
		if (t%100==0) cout << egyesek << endl;

		//nullázó, hogy többször lehessen futtatni a keresést anélkül hogy újraindítnánk a programot
		for (int i = 0; i < DEF_PROTEIN_NUMBER; i++) {
			protein[i].reset_dipole(DEF_DIPOL);
		}
	}
	fileki.close();
	fileki_bool.close();
	fileki_dipol.close();
	fileki_input.close();
	fileki_bin_input.close();

	
}

//f-re lefutó struktúra kereső fõfüggvény
void fofuggveny()
{


	protein[0].initialize_molekula(18, 18, 18, false, 2, false);
	protein[1].initialize_molekula(18, 18, 19, true, 1, false);
	protein[2].initialize_molekula(18, 19, 18, true, 0, false);
	protein[3].initialize_molekula(18, 19, 19, false, 2, false);
	protein[4].initialize_molekula(19, 18, 18, true, 0, false);
	protein[5].initialize_molekula(19, 18, 19, true, 1, false);
	protein[6].initialize_molekula(19, 19, 18, false, 2, false);
	protein[7].initialize_molekula(19, 19, 19, true, 0, true);





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