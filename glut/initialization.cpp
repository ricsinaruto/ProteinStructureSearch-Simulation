#include "screencasts.h"

//globális változók inicializálása
void initializeGlobals(void)

{
	/* WINDOW */
	windowHeight = DEF_WINDOW_HEIGHT;
	windowWidth = DEF_WINDOW_WIDTH;

	/* TOGGLE */
	toggleAxes = DEF_AXES;
	toggleParams = DEF_PARAMS;

	/* PROJECTION */
	dim = DEF_DIM;
	th = DEF_TH;
	ph = DEF_PH;
	fov = DEF_FOV;
	asp = DEF_ASP;
	ecX = DEF_ECX;
	ecY = DEF_ECY;
	ecZ = DEF_ECZ;

	/* LIGHTING */
	toggleLight = DEF_LIGHT;
	distance = DEF_DISTANCE;
	ambient = DEF_AMBIENT;
	diffuse = DEF_DIFFUSE;
	emission = DEF_EMISSION;
	specular = DEF_SPECULAR;
	shininess = DEF_SHININESS;
	lightY = DEF_L_Y;
	lightPh = DEF_L_PH;
	lightTh = DEF_L_TH;

	/* TEXTURES */
	currentTexture = TEX_DEFAULT;

	//ugyanazok mint screencasts.h-ban
	
	mouseBtnPressed = "Press the left mouse button before pressing the right one when using the q,w,e keys!";
	mouseState = "";
	Shift = "rotation";
	jobbx = 0;
	mouseX = 0, mouseY = 0;
	xcoord = 0, ycoord = 0;
	th2 = 0, ph2 = 0;
	lightTh2 = 0, lightPh2 = 0;
	fely = 0, timed = 0;
	ecX2 = 0, ecY2 = 0;
	main_window = 0;
	szamlalo = 0;
	ter = 0;
	enter = "pressed";
	valto = "parameters";
	kimenetek_szama = 6;
	bemenetek_szama = 3;

}

/*
* initializeObjs(void)
* ------
* Initializes all of our objects
*/
void initializeObjs(void)
{
	int hely=0;
	int szam = szamok[0] * 1296 + szamok[1] * 36 + szamok[2];

	//ha beírunk 3 koordinátát ez lefut és egy piros vagy zöld (térrel terhelt) kockát létrehozz az adott koordinátára
	for (int k = -18; k < 18; k++)
	{
		for (int j = -18; j < 18; j++)
		{
			for (int i = -18; i < 18; i++)
			{
				//ez megrajzolja a piros kockát
				if (szamlalo == 3&&szam==hely&&enter=="pressed")
				{
					cubes[szam] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 } ,{(float)szamok[0],(float)szamok[1],(float)szamok[2]}} };
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].van = true;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dip = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipA = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipB = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].ter = false;

					szamlalo = 0;
					szamok[0] = 36;
					szamok[1] = 36;
					szamok[2] = 36;
				}

				//ez visszaállítja az alap fekete kockát
				if (szamlalo == 3 && szam == hely&&enter == "delete")
				{
					cubes[szam] = { { { (float)i, (float)j, (float)k },{ 0,0,0 },{ 90,0,0 },{ (float)szamok[0],(float)szamok[1],(float)szamok[2] } } };
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].van = false;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].ter = false;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dip = dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipA = dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipB = 0;

					szamlalo = 0;
					szamok[0] = 36;
					szamok[1] = 36;
					szamok[2] = 36;
				}

				//ez megrajzolja a térrel terhelt zöld kockát
				if (szamlalo == 3 && szam == hely&&enter == "field")
				{
					cubes[szam] = { { { (float)i, (float)j, (float)k },{ 1,1,1 },{ 90,0,0 },{ (float)szamok[0],(float)szamok[1],(float)szamok[2] } } };
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].van = true;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dip = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipA = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].dipB = -100;
					dronpa[szamok[0] + 1][szamok[1] + 1][szamok[2] + 1].ter = true;

					szamlalo = 0;
					szamok[0] = 36;
					szamok[1] = 36;
					szamok[2] = 36;
				}
				hely++;
			}
		}
	}

	hely = 0;
	for (int k = -18; k < 18; k++)
	{
		for (int j = -18; j < 18; j++)
		{
			for (int i = -18; i < 18; i++)
			{
				h1[hely] = hely / 1296;
				h2[hely] = (hely - h1[hely] * 1296) / 36;
				h3[hely] = hely - h1[hely] * 1296 - h2[hely] * 36;

				hely++;
			}
		}
	}
}

//molekula rács létrehozása, alapból 0 minden molekula dipólja, és false az értéke
void initializeProteins(void) {
	for (int i = 0; i < 38; i++)
	{
		for (int j = 0; j < 38; j++)
		{
			for (int k = 0; k < 38; k++)
			{
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

//textúrák inicializálása
void initializeTextures(void) {
	//számok
	for (int i = 0; i < 10; i++) {
		std::string ok = std::to_string(i)+".bmp";
		char* c = &ok[0];
		textures[i] = loadTexBMP(c);
	}

	//karakterek
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