#include "main.h"


/* Initializes all of our objects */
void initializeObjs(void) {
	int position=0;
	int index = numbers[0] * 1296 + numbers[1] * 36 + numbers[2];	//switch to base 36

	/* If we specify 3 coordinates this will run and a red or green (field applied) cube will be created */
	for (int k = -18; k < 18; k++) {
		for (int j = -18; j < 18; j++) {
			for (int i = -18; i < 18; i++) {

				// this draws the red cube
				if (counter == 3 && index==position && enter=="pressed") {
					cubes[index] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 } ,
										{(float)numbers[0],(float)numbers[1],(float)numbers[2]}} };

					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].is = true;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dip = -100;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dipA = -100;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dipB = -100;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].field = false;

					counter = 0;
					numbers[0] = numbers[1] = numbers[2] = 36;
				}

				// this sets it back to black cube
				if (counter == 3 && index == position && enter == "delete") {
					cubes[index] = { { { (float)i, (float)j, (float)k },{ 0,0,0 },{ 90,0,0 },
									{ (float)numbers[0],(float)numbers[1],(float)numbers[2] } } };

					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].is = false;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].field = false;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dip = 0;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dipA = 0;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dipB = 0;

					counter = 0;
					numbers[0] = numbers[1] = numbers[2] = 36;
				}

				// this draws the green cube which has a field applied to it
				if (counter == 3 && index == position && enter == "field") {
					cubes[index] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 },
									{ (float)numbers[0],(float)numbers[1],(float)numbers[2] } } };

					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].is = true;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dip = -100;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dipA = -100;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].dipB = -100;
					dronpa[numbers[0] + 1][numbers[1] + 1][numbers[2] + 1].field = true;

					counter = 0;
					numbers[0] = numbers[1] = numbers[2] = 36;
				}

				// 3 arrays fot the 3 digits of base 36 coordinates (will be needed for textures)
				base36_1[position] = position / 1296;
				base36_2[position] = (position - base36_1[position] * 1296) / 36;
				base36_3[position] = position - base36_1[position] * 1296 - base36_2[position] * 36;
				position++;
			}
		}
	}
}

// initialize molecule grid, with 0 dipole moment, and false existing
void initializeProteins(void) {
	for (int i = 0; i < 38; i++) {
		for (int j = 0; j < 38; j++) {
			for (int k = 0; k < 38; k++) {
				dronpa[i][j][k].is = false;
				dronpa[i][j][k].field = 0;
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