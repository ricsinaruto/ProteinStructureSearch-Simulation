#include "screencasts.h"
using namespace std;

//tér nagyság megadása, i lenyomásakor fut le
/*
void ternagysag() {
	if (szamlalo == 3)
	{
		if (szamok[0]==0) ter = -(szamok[1] * 10 + szamok[2]);
		else ter = szamok[1] * 10 + szamok[2];
		
		//alaphelyzetre állítás
		szamlalo = 0;
		szamok[0] = 36;
		szamok[1] = 36;
		szamok[2] = 36;
	}
}*/

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




//összehasonlítást végez be1 és be2 között, 0 bemenetre igazt ad ha be1<be2, 1 bemenetre igazt ad ha be1>be2
bool hasonlitas(int be1, int be2, int kacsacsor)
{
	bool hasonlit=false;
	if ((be1) < be2 && kacsacsor==0) hasonlit = true;
	
	if (be1>(be2) && kacsacsor == 1) hasonlit = true;

	if (kacsacsor == 2) hasonlit = true;

	return hasonlit;
}

int factorial(int f)
{
	if (f == 0) return 1;
	return(f * factorial(f - 1));
}



//random number gen
double fRand(double fMin, double fMax)
{
	
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}

//logikai függvény összes sorát megnézi, hogy jó-e
bool logikai_hasonlitas(double **actual, double **desired) {

	bool jo = true;
	int bemenetek_iteralo = 0;
	int kimenetek_iteralo = 0;

	for (int i = 0; i < pow(2,bemenetek_szama); i++) {
		for (int j = 0; j < bemenetek_szama+kimenetek_szama; j++) {
			if (jo && dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].ter) {
				if (!hasonlitas(actual[i][j], desired[j][bemenetek[i][bemenetek_iteralo]], bemenetek[i][bemenetek_iteralo])) jo = false;
				bemenetek_iteralo++;
			}
			if (jo && !dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].ter && dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].kell) {
				if (!hasonlitas(actual[i][j], desired[j][kimenetek[kimenetek_iteralo][i]],kimenetek[kimenetek_iteralo][i])) jo = false;
				kimenetek_iteralo++;
			}
		}
		bemenetek_iteralo = kimenetek_iteralo = 0;
	}

	return jo;
}

//reset dipole moment
void reset_dipole(int i, int j, int k, double dipole) {
	dronpa[i][j][k].dip = dipole;
	dronpa[i][j][k].dipA = dipole;
	dronpa[i][j][k].dipB = dipole;
	dronpa[i][j][k].qeA = 0;
	dronpa[i][j][k].qeB = 0;
	dronpa[i][j][k].qp1A = 0;
	dronpa[i][j][k].qp1B = 0;
	dronpa[i][j][k].qp2A = 0;
	dronpa[i][j][k].qp2B = 0;
}

//tér keresés
void harmony_search() {

	//az első futást csak egyszer kell, és elmentjük az alap dipól értékeket
	futas();
	double *dipol = new double[struktura_szamlal];
	for (int i = 0; i < struktura_szamlal; i++) {
		dipol[i] = dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].dip;
	}

	//kívánt érték, alaptól 10 dipol eltérés 
	double **desired = new double*[bemenetek_szama+kimenetek_szama];
	for (int i = 0; i <bemenetek_szama+kimenetek_szama; i++) { desired[i] = new double[2]; }

	int kell_iterator = 0;
	for (int i = 0; i < struktura_szamlal; i++) {
		if (dronpa[itomb_mol[i]][jtomb_mol[i]][ktomb_mol[i]].kell) {
			for (int j = 0; j < 2; j++) {
				if (j == 0) desired[kell_iterator][j] = dipol[i] - tolerance;
				else desired[kell_iterator][j] = dipol[i] + tolerance;
			}
			kell_iterator++;
		}
	}
	

	//jelenlegi érték, első index a logikai sor, második a molekula száma
	double **actual = new double*[pow(2,bemenetek_szama)];
	for (int i = 0; i < pow(2, bemenetek_szama); i++) { actual[i] = new double[bemenetek_szama+kimenetek_szama]; }

	kell_iterator = 0;
	for (int i = 0; i < pow(2, bemenetek_szama); i++) {
		for (int j = 0; j < 3; j++) {
			if (dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].kell) {
				actual[i][kell_iterator] = dipol[j];
				kell_iterator++;
			}
		}
		kell_iterator = 0;
	}
	
	//random tér inicializálás
	double **inputTer = new double*[bemenetek_szama];
	for (int i = 0; i < bemenetek_szama; i++) { inputTer[i] = new double[2]; } //első index az input molekula száma, második index, hogy a 0 logikai értékű térről, vagy az 1 logikai értékű térről van-e szó
	
	

	//fitness
	double fitness = 0;
	double bestFitness = 100000000;
	int iteration = 0;
	bool hasonlit = false;

	//keresés
	while (!hasonlit) {
		
		//random tér inicializálás
		for (int i = 0; i < bemenetek_szama; i++) {
			for (int j = 0; j < 2; j++) {
				inputTer[i][j] = fRand(-max_ter,max_ter);
			}
		}


		/* SIMULATION */
		for (int i = 0; i < pow(2,bemenetek_szama); i++) {
			for (int j = 0; j < struktura_szamlal; j++) {
				reset_dipole(itomb_mol[j], jtomb_mol[j], ktomb_mol[j], dipol[j]);
			}

			//tér aplikálás
			int bemenet_iteralo = 0;
			for (int j = 0; j < struktura_szamlal; j++) {
				if (dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].ter) {
					dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].terMag = inputTer[bemenet_iteralo][bemenetek[i][bemenet_iteralo]];
						bemenet_iteralo++;
				}
			}
			std::string ok = "graf"+std::to_string(i) + ".csv";
			char* c = &ok[0];
			//futasv(c,true);
			futas();

			bemenet_iteralo = 0;
			for (int j = 0; j < struktura_szamlal; j++) {
				if (dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].ter) {
					dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].terMag = 0;
					bemenet_iteralo++;
				}
			}
			//futasv(c,false);
			futas();

		
			//dipól értékek elmentése
			int kell_iteralo = 0;
			for (int j = 0; j < struktura_szamlal; j++) {
				if (dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].kell) {
					actual[i][kell_iteralo] = dronpa[itomb_mol[j]][jtomb_mol[j]][ktomb_mol[j]].dip;
					kell_iteralo++;
				}
			}
		}
		

		//kiértékelés
		fitness = 0;
		//cout << endl;
		/*if (actual[0] > desired[0]) fitness += (actual[0] - desired[0])*(actual[0] - desired[0]);
		if (actual[2] > desired[2]) fitness += (actual[2] - desired[2])*(actual[2] - desired[2]);
		if (actual[1] < desired[1]) fitness += (actual[1] - desired[1])*(actual[1] - desired[1]);*/

		fitness = sqrt(fitness);
		
		if (fitness < bestFitness) bestFitness = fitness;
		//cout << fitness << endl;

		iteration++;

		//megfelelnek-e a logikai értékek
		hasonlit=logikai_hasonlitas(actual,desired);
		
		//cout << hasonlit << endl;
	}

	/* Adatok kiíratása */
	for (int i = 0; i < bemenetek_szama; i++) {
		for (int j = 0; j < 2; j++) {
			cout <<i<<". input molekulara "<<j<<" logikai ter nagysaga: "<< inputTer[i][j] << "   ";
		}
		cout << endl;
	}
	//megadja a próbálgatások számát
	cout <<"number of iterations: "<< iteration << endl<<endl;


	//nullázó, hogy többször lehessen futtatni a keresést anélkül hogy újraindítnánk a programot
	for (int i = 0; i < struktura_szamlal; i++) {
		reset_dipole(itomb_mol[i], jtomb_mol[i], ktomb_mol[i], -100);
	}

	/* tömb és mátrixok törlése */
	//for (int i = 0; i < bemenetek_szama + kimenetek_szama; i++) { delete[] desired[i]; }
	//delete[] desired;
	//desired = nullptr;

	//for (int i = 0; i < pow(2, bemenetek_szama); i++) { delete[] actual[i]; }
	//delete[] actual;
	//actual = nullptr;

	//for (int i = 0; i < bemenetek_szama; i++) { delete[] inputTer[i]; }
	//delete[] inputTer;
	//inputTer = nullptr;

	//delete[] dipol;
	//dipol = nullptr;
}

//f-re lefutó szimulációs fõfüggvény
void fofuggveny()
{
	int i = 18, j = 18, k = 18;
	struktura_szamlal = 4;	//ennyi molekula inicializálva

	//XOR és XNOR struktúrát még nem talált
	itomb_mol[0] = 17;
	jtomb_mol[0] = 17;
	ktomb_mol[0] = 18;
	dronpa[itomb_mol[0]][jtomb_mol[0]][ktomb_mol[0]].van = true;
	dronpa[itomb_mol[0]][jtomb_mol[0]][ktomb_mol[0]].kell = true;
	dronpa[itomb_mol[0]][jtomb_mol[0]][ktomb_mol[0]].dip = -100;
	dronpa[itomb_mol[0]][jtomb_mol[0]][ktomb_mol[0]].dipA = -100;
	dronpa[itomb_mol[0]][jtomb_mol[0]][ktomb_mol[0]].dipB = -100;
	dronpa[itomb_mol[0]][jtomb_mol[0]][ktomb_mol[0]].ter = true;
	
	
	itomb_mol[1] = 16;
	jtomb_mol[1] = 18;
	ktomb_mol[1] = 18;
	dronpa[itomb_mol[1]][jtomb_mol[1]][ktomb_mol[1]].van = true;
	dronpa[itomb_mol[1]][jtomb_mol[1]][ktomb_mol[1]].kell = true;
	dronpa[itomb_mol[1]][jtomb_mol[1]][ktomb_mol[1]].dip = -100;
	dronpa[itomb_mol[1]][jtomb_mol[1]][ktomb_mol[1]].dipA = -100;
	dronpa[itomb_mol[1]][jtomb_mol[1]][ktomb_mol[1]].dipB = -100;
	dronpa[itomb_mol[1]][jtomb_mol[1]][ktomb_mol[1]].ter = true;
	

	itomb_mol[2] = 18;
	jtomb_mol[2] = 18;
	ktomb_mol[2] = 18;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].van = true;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].kell = true;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].dip = -100;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].dipA = -100;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].dipB = -100;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].ter = false;


	itomb_mol[3] = 17;
	jtomb_mol[3] = 18;
	ktomb_mol[3] = 18;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].van = true;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].kell = false;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].dip = -100;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].dipA = -100;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].dipB = -100;
	dronpa[itomb_mol[2]][jtomb_mol[2]][ktomb_mol[2]].ter = false;
	
	
	

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

//beolvasása a struktúrának
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

// q,w,e billentyûk lenyomásakor a mozgatások és forgások kezelése
void windowPmotion(int x, int y)
{
	mouseX = x;
	mouseY = y;

	//fény forgatása
	if (Shift == "light")
	{
		if (mouseBtnPressed == "Left")
		{
			xcoord = mouseX;
			ycoord = mouseY;
			lightTh2 = lightTh;
			lightPh2 = lightPh;
			ecX2 = ecX;
			ecY2 = ecY;
			th2 = th;
			ph2 = ph;
		}

		if (mouseBtnPressed == "Right")
		{
			lightTh = (lightTh2 + (mouseX - xcoord));
			lightPh = (lightPh2 + (mouseX - xcoord));
		}
	}

	//struktúra forgatása
	else if (Shift=="rotation")
	{
		if (mouseBtnPressed == "Left")
		{
			xcoord = mouseX;
			ycoord = mouseY;
			lightTh2 = lightTh;
			lightPh2 = lightPh;
			th2 = th;
			ph2 = ph;
			ecX2 = ecX;
			ecY2 = ecY;
		}

		if (mouseBtnPressed == "Right")
		{
			th = (th2 + (mouseX - xcoord));
			ph = (ph2 + (mouseY - ycoord));
		}
	}

	//struktúra mozgatása
	else if (Shift=="movement")
	{
		if (mouseBtnPressed == "Left")
		{
			xcoord = mouseX;
			ycoord = mouseY;
			lightTh2 = lightTh;
			lightPh2 = lightPh;
			ecX2 = ecX;
			ecY2 = ecY;
			th2 = th;
			ph2 = ph;
		}

		if (mouseBtnPressed == "Right")
		{
			ecX = ecX2 + (-mouseX + xcoord)/10;
			ecY = ecY2 + (mouseY - ycoord)/10;
		}
	}

	th %= 360;
	ph %= 360;
	redisplayAll();
}

//zoomolás
void mouseWheel(int scroll, int dir, int x, int y)
{
	dim -= (double)dir;
	redisplayAll();
}

//egér gombok lenyomásának kezelése
void windowMouse(int btn, int state, int x, int y)
{
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
void windowKey(unsigned char key, int x, int y)
{
	/*  Exit on ESC */
	if (key == 27) exit(0);

	//space
	if (key == 32)
	{
		if (valto == "parameters") valto = "coordinates";
		else valto = "parameters";
	}

	//paraméterek megadása
	if (valto == "parameters")
	{
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
	if (valto == "coordinates")
	{
		for (int i = 97; i < 123; i++)
		{
			if (key == i) { szamok[szamlalo] = i-87; szamlalo++; key = ','; }
		}

		for (int i = 48; i < 58; i++)
		{
			if (key == i) { szamok[szamlalo] = i - 48; szamlalo++; key = ','; }
		}
	}

	//0-9ig koordináta megadás
	for (int i = 48; i < 58; i++)
	{
		if (key == i) { szamok[szamlalo] = i - 48; szamlalo++; key = ','; }
	}

	/*  Translate shininess power to value (-1 => 0) */
	shinyvec[0] = shininess<0 ? 0 : pow(2.0, shininess);

	redisplayAll();
}

/*
*  windowMenu
*  ------
*  Window menu is the same as the keyboard clicks
*/
void windowMenu(int value)
{
	windowKey((unsigned char)value, 0, 0);
}

/*
*  windowSpecial()
*  ------
*  GLUT calls this routine when an arrow key is pressed
*/
void windowSpecial(int key, int x, int y)
{
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