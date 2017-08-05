#include "screencasts.h"


/* Initializes all of our objects */
void initializeObjs(void) {
	int hely=0;
	int szam = szamok[0] * 1296 + szamok[1] * 36 + szamok[2];	//switch to base 36

	/* If we specify 3 coordinates this will run and a red or green (field applied) cube will be created */
	for (int k = -18; k < 18; k++) {
		for (int j = -18; j < 18; j++) {
			for (int i = -18; i < 18; i++) {

				//this draws the red cube
				if (szamlalo == 3&&szam==hely&&enter=="pressed") {
					cubes[szam] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 } ,
										{(float)szamok[0],(float)szamok[1],(float)szamok[2]}} };

					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].van = true;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dip = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipA = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipB = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].ter = false;

					szamlalo = 0;
					szamok[0] = szamok[1] = szamok[2] = 36;
				}

				//this sets it back to black cube
				if (szamlalo == 3 && szam == hely&&enter == "delete") {
					cubes[szam] = { { { (float)i, (float)j, (float)k },{ 0,0,0 },{ 90,0,0 },
									{ (float)szamok[0],(float)szamok[1],(float)szamok[2] } } };

					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].van = false;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].ter = false;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dip = 0;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipA = 0;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipB = 0;

					szamlalo = 0;
					szamok[0] = szamok[1] = szamok[2] = 36;
				}

				//this draws the green cube which has a field applied to it
				if (szamlalo == 3 && szam == hely&&enter == "field") {
					cubes[szam] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 },
									{ (float)szamok[0],(float)szamok[1],(float)szamok[2] } } };

					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].van = true;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dip = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipA = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipB = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].ter = true;

					szamlalo = 0;
					szamok[0] = szamok[1] = szamok[2] = 36;
				}

				//3 array fot the 3 coordinates (will be needed for textures)
				h1[hely] = hely / 1296;
				h2[hely] = (hely - h1[hely] * 1296) / 36;
				h3[hely] = hely - h1[hely] * 1296 - h2[hely] * 36;
				hely++;
			}
		}
	}
}

//initialize molecule grid, with 0 dipole moment, and false existing
void initializeProteins(void) {
	for (int i = 0; i < 38; i++) {
		for (int j = 0; j < 38; j++) {
			for (int k = 0; k < 38; k++) {
				dronpa[i][j][k].van = false;
				dronpa[i][j][k].ter = 0;
				dronpa[i][j][k].dipA = 0;
				dronpa[i][j][k].dipB = 0;
				dronpa[i][j][k].dip = 0;
				dronpa[i][j][k].qeA = 0;
				dronpa[i][j][k].qeB = 0;
				dronpa[i][j][k].qp1A = 0;
				dronpa[i][j][k].qp1B = 0;
				dronpa[i][j][k].qp2A = 0;
				dronpa[i][j][k].qp2B = 0;
			}
		}
	}
}

//initialize textures
void initializeTextures(void) {
	//numbers
	for (int i = 0; i < 10; i++) {
		std::string ok = std::to_string(i)+".bmp";
		char* c = &ok[0];
		textures[i] = loadTexBMP(c);
	}

	//characters
	for (int i = 10; i < 36; i++) {
		std::stringstream ss;
		std::string ok;
		char v = i+87;

		ss << v;
		ss >> ok;
		ok += ".bmp";
		char* c = &ok[0];
		textures[i] = loadTexBMP(c);
	}
}