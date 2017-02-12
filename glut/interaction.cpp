﻿#include "screencasts.h"
using namespace std;


//egyeneletkhez paraméterek
double dt = 0.1;
double t = 0;
double K = -0.07;
double Ce = 0.18;
double tav = 343;
double taue = 1.51;
double Cp = 0.153;
double taup = 5.7;

double U = 0;
double Ce1 = 0.008;
double Ce2 = 180;
double Cp1 = 0.037;
double Cp2 = 153;

int *itomb_mol = new int[1000];
int *jtomb_mol = new int[1000];
int *ktomb_mol = new int[1000];


int struktura_szamlal = 0;





//bemenetek a proba függvényhez
int bemenetek[16][4] = { { 0,0,0,0 },{ 1,0,0,0 },{ 0,1,0,0 },{ 1,1,0,0 },{ 0,0,1,0 },{ 1,0,1,0 },{ 0,1,1,0 },{ 1,1,1,0 },
{ 0,0,0,1 },{ 1,0,0,1 },{ 0,1,0 ,1 },{ 1,1,0,1 },{ 0,0,1,1 },{ 1,0,1,1 },{ 0,1,1,1 },{ 1,1,1,1 } };



//tér nagyság megadása, i lenyomásakor fut le
void ternagysag()
{
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
}

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
void futasv()
{
	ofstream fileki;
	ifstream filebe;
	int i, j, k, l, n;
	string sor;
	n = 500;
	int szamlal = 0;
	if (!t) fileki.open("tmp.csv");
	else fileki.open("tmp.csv", ios::app);

	
	
	
	if (!t) fileki << ",";
	for (i = 1; i <= 36; i++)
	{
		for (j = 1; j <= 36; j++)
		{
			for (k = 1; k <= 36; k++)
			{
				if (dronpa[i][j][k].van)
				{
					if (!t) fileki << i - 1 << " " << j - 1 << " " << k - 1 << ",";
					szamlal++;
				}
			}
		}
	}
	if (!t) fileki << endl;

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
	fileki.open("graf.csv");
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
	n = 500;
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

//-100-ra állítja dipól értékeket
void wipe()
{
	int i, j, k;
	for (i = 0; i < 36; i++)
	{
		for (j = 0; j < 36; j++)
		{
			for (k = 0; k < 36; k++)
			{

				dronpa[i][j][k].dipA = -100;
				dronpa[i][j][k].dipB = -100;
				dronpa[i][j][k].dip = -100;
			}
		}
	}
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

//összehasonlítást végez be1 és be2 között, 0 bemenetre igazt ad ha be1<be2-10, 1 bemenetre igazt ad ha be1>be2+10
bool hasonlitas(int be1, int be2, int kacsacsor)
{
	bool hasonlit=false;
	if ((be1+1) < be2 && kacsacsor==0) hasonlit = true;
	
	if (be1>(be2+1) && kacsacsor == 1) hasonlit = true;

	if (kacsacsor == 2) hasonlit = true;

	return hasonlit;
}

int factorial(int f)
{
	if (f == 0) return 1;
	return(f * factorial(f - 1));
}

//ez tart a legtöbb ideig, ebben van egy adott struktúrán belül a különbözõ terekkel való tesztelés
bool terteszt(int molekulaSzam,int bemenetek_szam, int kimenetek_szam)
{
	//cout << "terteszt kezd" << endl;
	int **kimenet = new int*[kimenetek_szama];
	for (int i = 0; i < kimenetek_szama; i++) { kimenet[i] = new int[pow(2, bemenetek_szama)]; }

	int **kimenetek_t = new int*[kimenetek_szam];
	for (int i = 0; i < kimenetek_szam; i++) { kimenetek_t[i] = new int[pow(2, bemenetek_szam)]; }
	for (int i = 0; i < kimenetek_szam; i++)
	{
		for (int j = 0; j < pow(2, bemenetek_szam); j++) 
		{ 
			kimenetek_t[i][j] = kimenetek[i][j];
		}
	}

	int **bemenetek_t = new int*[pow(2, bemenetek_szam)];
	for (int i = 0; i < pow(2, bemenetek_szam); i++) { bemenetek_t[i] = new int[bemenetek_szam]; }
	for (int i = 0; i < pow(2, bemenetek_szam); i++)
	{
		for (int j = 0; j < bemenetek_szam; j++) { bemenetek_t[i][j] = bemenetek[i][j]; }
	}
	

	int *itomb = new int[molekulaSzam];
	int *jtomb = new int[molekulaSzam];
	int *ktomb = new int[molekulaSzam];
	int l = 0;
	int m = 0;
	int n = 0;
	int p = 0;
	bool sikerult=false;
	double *dipolmoment = new double[molekulaSzam];

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
	int *j_elmentes=new int[kimenetek_szam];
	int szimulacioszam = 0;
	int lehetosegek = 0;
	int tr1 = 0, tr2 = 0, tr3 = 0, tr4 = 0;
	int l0, m0, n0, p0, tr10, tr20, tr30, tr40;

	ofstream fileki;
	string fileki_string = "";
	string fileki_vegso = "";
	int ter_max = 10;
	if (bemenetek_szam == 4) { l0 = 0, m0 = 1, n0 = 2, p0 = 3; tr10 = -ter_max, tr20 = -ter_max, tr30 = -ter_max, tr40 = -ter_max; }
	else if (bemenetek_szam==3) { l0 = molekulaSzam-4, m0 = 0, n0 = 1, p0 = 2; tr10 = -ter_max, tr20 = -ter_max, tr30 = -ter_max, tr40 = ter_max;}
	else if (bemenetek_szam == 2) { l0 = molekulaSzam - 4, m0 = molekulaSzam-3, n0 = 0, p0 = 1;
	tr10 = -ter_max, tr20 = -ter_max, tr30 = ter_max, tr40 = ter_max;}
	else if (bemenetek_szam == 1) { l0 = molekulaSzam - 4, m0 = molekulaSzam-3, n0 = molekulaSzam-2, p0 = 0;
	tr10 = -ter_max, tr20 = ter_max, tr30 = ter_max, tr40 = ter_max;}

	for (l = l0; l < molekulaSzam - 3; l++)
	{
		for (m = m0+l-l0; m < molekulaSzam - 2; m++)
		{
			for (n = n0+m-m0; n < molekulaSzam - 1; n++)
			{
				for (p = p0 + n - n0; p < molekulaSzam; p++)
				{
					
					
					
					int **tesztelt_kimenetek = new int*[kimenetek_szam];
					for (int lll = 0; lll < kimenetek_szam; lll++) 
					{  
						tesztelt_kimenetek[lll] = new int[molekulaSzam - bemenetek_szam]; 
						for (int cuck = 0; cuck < molekulaSzam - bemenetek_szam; cuck++) { tesztelt_kimenetek[lll][cuck] = -1; }
					}
					
					for (int megegy = 0; megegy < molekulaSzam - bemenetek_szam; megegy++)
					{
						for (int megegy_kimenet=0; megegy_kimenet < kimenetek_szam; megegy_kimenet++)
						{
							fileki_string = "";
							fileki.open("talalatok.txt");
							lehetosegek = 0;
							int lul = 0;
							for (int z = 0; z < pow(2, bemenetek_szam); z++)
							{
								for (tr1 = tr10; tr1 <= ter_max; tr1++)
								{
									for (tr2 = tr20; tr2 <= ter_max; tr2++)
									{
										for (tr3 = tr30; tr3 <= ter_max; tr3++)
										{
											for (tr4 = tr40; tr4 <= ter_max; tr4++)
											{
												for (int eszter = 0; eszter < factorial(bemenetek_szam); eszter++)
												{
													int l1, l3, l5, l7;
													switch (eszter)
													{
													case 0: l1 = 1, l3 = 3, l5 = 5, l7 = 7;
														break;
													case 1: l1 = 3, l3 = 1, l5 = 5, l7 = 7;
														break;
													case 2: l1 = 5, l3 = 1, l5 = 3, l7 = 7;
														break;
													case 3: l1 = 5, l3 = 3, l5 = 1, l7 = 7;
														break;
													case 4: l1 = 3, l3 = 5, l5 = 1, l7 = 7;
														break;
													case 5: l1 = 1, l3 = 5, l5 = 3, l7 = 7;
														break;

													case 6: l1 = 7, l3 = 3, l5 = 5, l7 = 1;
														break;
													case 7: l1 = 7, l3 = 1, l5 = 5, l7 = 3;
														break;
													case 8: l1 = 7, l3 = 1, l5 = 3, l7 = 5;
														break;
													case 9: l1 = 7, l3 = 3, l5 = 1, l7 = 5;
														break;
													case 10: l1 = 7, l3 = 5, l5 = 1, l7 = 3;
														break;
													case 11: l1 = 7, l3 = 5, l5 = 3, l7 = 1;
														break;

													case 12: l1 = 3, l3 = 7, l5 = 5, l7 = 1;
														break;
													case 13: l1 = 1, l3 = 7, l5 = 5, l7 = 3;
														break;
													case 14: l1 = 1, l3 = 7, l5 = 3, l7 = 5;
														break;
													case 15: l1 = 3, l3 = 7, l5 = 1, l7 = 5;
														break;
													case 16: l1 = 5, l3 = 7, l5 = 1, l7 = 3;
														break;
													case 17: l1 = 5, l3 = 7, l5 = 3, l7 = 1;
														break;

													case 18: l1 = 3, l3 = 5, l5 = 7, l7 = 1;
														break;
													case 19: l1 = 1, l3 = 5, l5 = 7, l7 = 3;
														break;
													case 20: l1 = 1, l3 = 3, l5 = 7, l7 = 5;
														break;
													case 21: l1 = 3, l3 = 1, l5 = 7, l7 = 5;
														break;
													case 22: l1 = 5, l3 = 1, l5 = 7, l7 = 3;
														break;
													case 23: l1 = 5, l3 = 3, l5 = 7, l7 = 1;
														break;

													}

													int szamlal = 0;
													for (int i = 0; i < bemenetek_szam * 2 + 1; i++)
													{
														if (i == 0)
														{
															futas();
															for (int j = 0; j < molekulaSzam; j++)
															{
																dipolmoment[j] = dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip;
															}
														}
														if (i == l1)
														{
															ter = tr1;
															dronpa[itomb[p]][jtomb[p]][ktomb[p]].ter = true;
															futas();
														}
														if (i == l1 + 1)
														{
															dronpa[itomb[p]][jtomb[p]][ktomb[p]].ter = false;
															futas();
														}
														if (i == l3)
														{
															ter = tr2;
															dronpa[itomb[n]][jtomb[n]][ktomb[n]].ter = true;
															futas();
														}
														if (i == l3 + 1)
														{
															dronpa[itomb[n]][jtomb[n]][ktomb[n]].ter = false;
															futas();
														}
														if (i == l5)
														{
															ter = tr3;
															dronpa[itomb[m]][jtomb[m]][ktomb[m]].ter = true;
															futas();
														}
														if (i == l5 + 1)
														{
															dronpa[itomb[m]][jtomb[m]][ktomb[m]].ter = false;
															futas();
														}
														if (i == l7)
														{
															ter = tr4;
															dronpa[itomb[l]][jtomb[l]][ktomb[l]].ter = true;
															futas();
														}
														if (i == l7 + 1)
														{
															dronpa[itomb[l]][jtomb[l]][ktomb[l]].ter = false;
															futas();
														}
													}

													bool *figyel = new bool[kimenetek_szam];
													for (int i = 0; i < kimenetek_szam; i++) { figyel[i] = false; }

													for (int j = 0; j < molekulaSzam; j++)
													{
														bool kimenet_e = true;
														for (int i = 0; i <= z; i++)
														{
															if (z == i)
															{
																if (bemenetek_szam >= 1) {
																	if (itomb[p] == itomb[j] && jtomb[p] == jtomb[j] && ktomb[p] == ktomb[j]) {
																		if (hasonlitas(dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip,
																			dipolmoment[j], bemenetek_t[z][0]))
																			szamlal++;
																		kimenet_e = false;
																	}

																	if (bemenetek_szam >= 2) {
																		if (itomb[n] == itomb[j] && jtomb[n] == jtomb[j] && ktomb[n] == ktomb[j]) {
																			if (hasonlitas(dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip,
																				dipolmoment[j], bemenetek_t[z][1]))
																				szamlal++;
																			kimenet_e = false;
																		}

																		if (bemenetek_szam >= 3) {
																			if (itomb[m] == itomb[j] && jtomb[m] == jtomb[j] && ktomb[m] == ktomb[j]) {
																				if (hasonlitas(dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip,
																					dipolmoment[j], bemenetek_t[z][2]))
																					szamlal++;
																				kimenet_e = false;
																			}

																			if (bemenetek_szam >= 4) {
																				if (itomb[l] == itomb[j] && jtomb[l] == jtomb[j] && ktomb[l] == ktomb[j]) {
																					if (hasonlitas(dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip,
																						dipolmoment[j], bemenetek_t[z][3]))
																						szamlal++;
																					kimenet_e = false;
																				}
																			}
																		}
																	}
																}

																//ha nem bemenet
																if (kimenet_e)
																{
																	for (int b = 0; b < kimenetek_szam; b++)
																	{
																		bool teszteles = true;
																		for (int najo = 0; najo < molekulaSzam - bemenetek_szam; najo++)
																		{
																			if (tesztelt_kimenetek[b][najo] == j && z==0) teszteles = false;
																		}
																		if (hasonlitas(dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip, dipolmoment[j], kimenetek_t[b][z])) {
																			bool nemjo = false;
																			for (int a = 0; a < z; a++)
																			{
																				if (kimenet[b][a] != j) nemjo = true;
																			}

																			if (!figyel[b] && !nemjo) {
																				tesztelt_kimenetek[megegy_kimenet][megegy] = j;

																				if (teszteles)
																				{
																					figyel[b] = true;
																					j_elmentes[b] = j;
																					kimenet[b][z] = j;
																				}
																			}
																		}
																	}
																}
															}
														}


														if (lul == 0) {

															fileki_string += (to_string(itomb[j]) + " ");
															fileki_string += (to_string(jtomb[j]) + " " + to_string(ktomb[j]));
															fileki_string += "     ";

															fileki << itomb[j] << " " <<
																jtomb[j] << " " << ktomb[j] << endl;
														}
														//cout << dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip << endl;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].dip = -100;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].dipA = -100;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].dipB = -100;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].qeA = 0;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].qeB = 0;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].qp1A = 0;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].qp1B = 0;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].qp2A = 0;
														dronpa[itomb[j]][jtomb[j]][ktomb[j]].qp2B = 0;

													}
													lul++;
													bool jo = true;
													for (int i = 0; i < kimenetek_szam; i++)
													{
														if (!figyel[i]) jo = false;
													}
													if (szamlal >= bemenetek_szam && jo)
													{
														int faszombamar = 0;
														if (tr1 == 0)faszombamar++;
														if (tr2 == 0)faszombamar++;
														if (tr3 == 0)faszombamar++;
														if (tr4 == 0)faszombamar++;
														if (faszombamar < bemenetek_szam) {

															for (int i = 0; i < bemenetek_szam; i++)
															{
																int lol = 0;
																int terr = 0;
																if (i == 0) {
																	lol = p; terr = tr1;
																}
																if (i == 1) {
																	lol = n; terr = tr2;
																}
																if (i == 2) {
																	lol = m; terr = tr3;
																}
																if (i == 3) {
																	lol = l; terr = tr4;
																}

																fileki_string += (string("field ") + to_string(i + 1));
																fileki_string += (string(": ") + to_string(terr)+"     ");
																fileki_string += (string("field ") + to_string(i + 1) + string(" applied to protein: ") + to_string(itomb[lol]) + " " +
																	to_string(jtomb[lol]) + " " + to_string(ktomb[lol]) + "     ");

																fileki << "field " << i + 1 << ": " << terr << endl;
																fileki << "field " << i + 1 << " applied to protein: " << itomb[lol] << " " <<
																	jtomb[lol] << " " << ktomb[lol] << endl;
																//cout << "ter " << i + 1 << ": " << terr << endl;
																//cout << "ter " << i + 1 << "-el terhelt molekula: " << itomb[lol] << " " <<
																	//jtomb[lol] << " " << ktomb[lol] << endl;

															}

															for (int i = 0; i < kimenetek_szam; i++)
															{

																fileki_string+= (string("output protein: ") + to_string(itomb[j_elmentes[i]]) + string(" ") +
																	to_string(jtomb[j_elmentes[i]]) + " " + to_string(ktomb[j_elmentes[i]]) + "     ");

																fileki << "output protein: " << itomb[j_elmentes[i]] << " " <<
																	jtomb[j_elmentes[i]] << " " << ktomb[j_elmentes[i]] << endl;
																//cout << "kimenet molekula: " << itomb[j_elmentes[i]] << " " <<
																	//jtomb[j_elmentes[i]] << " " << ktomb[j_elmentes[i]] << endl;
															}

															//cout << endl;
															lehetosegek++;
															//cout << lehetosegek << endl;
															tr1 = 10;
															tr2 = 10;
															tr3 = 10;
															tr4 = 10;
															eszter = 30;
														}
													}
													if (lehetosegek < z)
													{
														z = pow(2, bemenetek_szam);
														tr1 = 10;
														tr2 = 10;
														tr3 = 10;
														tr4 = 10;
														eszter = 30;
													}
													delete[] figyel;
												}
												szimulacioszam++;
												//cout << szimulacioszam << endl;
											}
										}
									}
								}


							}
							bool idonteven = true;
							for (int i = 0; i < kimenetek_szam; i++)
							{
								for (int v = 0; v < pow(2, bemenetek_szam) - 1; v++)
								{
									if (kimenet[i][v] != kimenet[i][v + 1]) idonteven = false;
									if (kimenetek_szam >= 2 && i < kimenetek_szam - 1)
									{
										if (kimenet[i][v] == kimenet[i + 1][v]) idonteven = false;
									}
								}
							}
							if (lehetosegek == pow(2, bemenetek_szam) && idonteven)
							{

								cout << "sikerult" << endl;
								m = molekulaSzam + 1;
								l = molekulaSzam + 1;
								p = molekulaSzam + 1;
								n = molekulaSzam + 1;
								megegy = 1000;
								megegy_kimenet = 1000;
								//fileki.open("talalatok.txt");
								//fileki << fileki_string;
								//fileki.close();
								sikerult = true;
							}
							//fileki << fileki_string;
							fileki.close();
						}
					}
					for (int llll = 0; llll < kimenetek_szam; llll++)
					{
						delete[] tesztelt_kimenetek[llll];
					}
					delete[] tesztelt_kimenetek;
				}
			}
		}
	}

	
	for (int i = 0; i < kimenetek_szama; i++) { delete[] kimenet[i]; }
	delete[] kimenet;

	
	for (int i = 0; i < kimenetek_szam; i++) { delete[] kimenetek_t[i]; }
	delete[] kimenetek_t;
	

	for (int i = 0; i < pow(2, bemenetek_szam); i++) { delete[] bemenetek_t[i]; }
	delete[] bemenetek_t;

	delete[] itomb;
	delete[] jtomb;
	delete[] ktomb;
	delete[] dipolmoment;
	delete[] j_elmentes;
	
	
	//cout << "terteszt veg" << endl;
	return sikerult;
}

//f-re lefutó szimulációs függvény
bool proba(int molekulaSzam)
{
	
	int eddigiMolekulak = molekulaSzam-1;
	int i = 18, j = 18, k = 18;
	bool sikerult = false;

	
	int hasonlit_int = 1;
	int **osszehasonlito = new int*[1000];
	for (int i = 0; i < 1000; i++) { osszehasonlito[i] = new int[molekulaSzam]; }
	osszehasonlito[0]=grafszam(eddigiMolekulak);
	

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

						if (p==hasonlit_int)
						{
							//run
							if (bemenetek_szama+kimenetek_szama<=molekulaSzam)	sikerult=terteszt(eddigiMolekulak,bemenetek_szama,kimenetek_szama);
							if (sikerult)
							{
								i = 18 + molekulaSzam;
								j = 18 + molekulaSzam;
								k = 18 + molekulaSzam;
							}

							osszehasonlito[hasonlit_int] = grafszam(eddigiMolekulak);
							if (!sikerult) cout << "egyszer" << endl;

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
						
						eddigiMolekulak = molekulaSzam-1;
					}
				}
			}
		}
	}
	
	for (int i = 0; i < 1000; i++) { delete[] osszehasonlito[i]; }
	delete[] osszehasonlito;
	
	

	return sikerult;
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

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			if (jo && j==0) {
				if (!hasonlitas(actual[i][j], desired[j][bemenetek[i][0]], bemenetek[i][0])) jo = false;
			}
			if (jo && j == 1) {
				if (!hasonlitas(actual[i][j], desired[j][kimenetek[0][i]],kimenetek[0][i])) jo = false;
			}
			if (jo && j == 2) {
				if (hasonlitas(actual[i][j], desired[j][bemenetek[i][1]], bemenetek[i][1])) jo = false;
			}
		}
	}

	return jo;
}

//tér keresés
void harmony_search() {

	//az első futást csak egyszer kell, és elmentjük az alap dipól értékeket
	futas();
	double *dipol = new double[3];
	for (int i = 0; i < 3; i++) {
		dipol[i] = dronpa[i+17][18][18].dip;
	}

	//kívánt érték
	double **desired = new double*[3];
	for (int i = 0; i < 3; i++) { desired[i] = new double[2]; }
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			if (j==0) desired[i][j] = dipol[i]-1;
			else desired[i][j] = dipol[i] + 1;
		}
	}
	

	//jelenlegi érték, első index a logikai sor, második a molekula száma
	double **actual = new double*[4];
	for (int i = 0; i < 4; i++) { actual[i] = new double[3]; }
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 3; j++) {
			actual[i][j] = dipol[j];
		}
	}
	
	//random tér inicializálás
	double **inputTer = new double*[2];
	for (int i = 0; i < 2; i++) { inputTer[i] = new double[2]; } //első index az input molekula száma, második index, hogy a 0 logikai értékű térről, vagy az 1 logikai értékű térről van-e szó
	
	
	

	//fitness
	double fitness = 0;
	double bestFitness = 100000000;
	int iteration = 0;
	bool hasonlit = false;

	//keresés
	while (!hasonlit) {
		
		//random tér inicializálás
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				inputTer[i][j] = fRand(-10,10);
			}
		}



		//szimuláció
		for (int j = 0; j < 4; j++) {
			for (int i = 0; i < 3; i++) {
				dronpa[i + 17][18][18].dip = dipol[i];
				dronpa[i + 17][18][18].dipA = dipol[i];
				dronpa[i + 17][18][18].dipB = dipol[i];
				dronpa[i + 17][18][18].qeA = 0;
				dronpa[i + 17][18][18].qeB = 0;
				dronpa[i + 17][18][18].qp1A = 0;
				dronpa[i + 17][18][18].qp1B = 0;
				dronpa[i + 17][18][18].qp2A = 0;
				dronpa[i + 17][18][18].qp2B = 0;
			}

			
			dronpa[17][18][18].terMag = inputTer[0][bemenetek[j][0]];
			dronpa[19][18][18].terMag = inputTer[1][bemenetek[j][1]];
			//cout << ter0 << "   " << ter2 << endl;
			futas();

			dronpa[17][18][18].terMag = 0;
			dronpa[19][18][18].terMag = 0;
			futas();

			for (int i = 0; i < 3; i++) {
				actual[j][i] = dronpa[i + 17][18][18].dip;
				//cout << actual[i] << "  ";
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
		cout << hasonlit << endl;

	}

	//megadja a próbálgatások számát
	cout << iteration << endl;


	//nullázó, hogy többször lehessen futtatni a keresést anélkül hogy újraindítnánk a programot
	for (int i = 0; i < 3; i++) {
		dronpa[i + 17][18][18].dip = -100;
		dronpa[i + 17][18][18].dipA = -100;
		dronpa[i + 17][18][18].dipB = -100;
		dronpa[i + 17][18][18].qeA = 0;
		dronpa[i + 17][18][18].qeB = 0;
		dronpa[i + 17][18][18].qp1A = 0;
		dronpa[i + 17][18][18].qp1B = 0;
		dronpa[i + 17][18][18].qp2A = 0;
		dronpa[i + 17][18][18].qp2B = 0;
	}



}

//f-re lefutó szimulációs fõfüggvény
void fofuggveny()
{
	ofstream fileki;

	//NAND struktúra, korlátok bemeneti térre: -10,+10..... kimeneti treshholdok: -10,+10 dipól eltérés az alaptól
	bool sikerult = false;
	int i = 18, j = 18, k = 18;
	dronpa[i-1][j][k].van = true;
	dronpa[i-1][j][k].dip = -100;
	dronpa[i-1][j][k].dipA = -100;
	dronpa[i-1][j][k].dipB = -100;
	dronpa[i-1][j][k].ter = true;
	

	dronpa[i][j][k].van = true;
	dronpa[i][j][k].dip = -100;
	dronpa[i][j][k].dipA = -100;
	dronpa[i][j][k].dipB = -100;
	dronpa[i][j][k].ter = false;

	dronpa[i+1][j][k].van = true;
	dronpa[i+1][j][k].dip = -100;
	dronpa[i+1][j][k].dipA = -100;
	dronpa[i+1][j][k].dipB = -100;
	dronpa[i+1][j][k].ter = true;
	

	itomb_mol[0] = i-1;
	jtomb_mol[0] = j;
	ktomb_mol[0] = k;

	itomb_mol[1] = i;
	jtomb_mol[1] = j;
	ktomb_mol[1] = k;

	itomb_mol[2] = i+1;
	jtomb_mol[2] = j;
	ktomb_mol[2] = k;


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
			enter = "field magnitude";  ternagysag();
		}

		//kimenet megadás
		else if (key == 'k') {
			enter = "output"; kimenet();
		}

		//szimuláció futtatása
		else if (key == 'r') futasv();

		//próba függvény
		else if (key == 'f') fofuggveny();

		//dipólok törlése
		else if (key == 'p') wipe();
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