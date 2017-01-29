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
	toggleAnimation = DEF_ANIMATE;
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
}


//textúrák inicializálása
void initializeTextures(void)
{
	/*
	TEX_DEFAULT 0
	TEX_BRICK 1
	TEX_CRATE 2
	TEX_ICE 3
	TEX_FIRE 4
	TEX_EARTH 5
	TEX_WOOD 6
	TEX_VENUS 7
	*/
	textures[TEX_BRICK] = loadTexBMP("txBrick14.bmp");
	textures[TEX_CRATE] = loadTexBMP("txCrate.bmp");
	textures[TEX_ICE] = loadTexBMP("txIce7.bmp");
	textures[TEX_FIRE] = loadTexBMP("txLava1.bmp");
	textures[TEX_EARTH] = loadTexBMP("txRock5.bmp");
	textures[TEX_WOOD] = loadTexBMP("txWood3.bmp");
	textures[TEX_VENUS] = loadTexBMP("txVenus.bmp");

	textures[0] = loadTexBMP("0.bmp");
	textures[1] = loadTexBMP("1.bmp");
	textures[2] = loadTexBMP("2.bmp");
	textures[3] = loadTexBMP("3.bmp");
	textures[4] = loadTexBMP("4.bmp");
	textures[5] = loadTexBMP("5.bmp");
	textures[6] = loadTexBMP("6.bmp");
	textures[7] = loadTexBMP("7.bmp");
	textures[8] = loadTexBMP("8.bmp");
	textures[9] = loadTexBMP("9.bmp");

	textures[10] = loadTexBMP("a.bmp");
	textures[11] = loadTexBMP("b.bmp");
	textures[12] = loadTexBMP("c.bmp");
	textures[13] = loadTexBMP("d.bmp");
	textures[14] = loadTexBMP("e.bmp");
	textures[15] = loadTexBMP("f.bmp");
	textures[16] = loadTexBMP("g.bmp");
	textures[17] = loadTexBMP("h.bmp");
	textures[18] = loadTexBMP("i.bmp");
	textures[19] = loadTexBMP("j.bmp");
	textures[20] = loadTexBMP("k.bmp");
	textures[21] = loadTexBMP("l.bmp");
	textures[22] = loadTexBMP("m.bmp");
	textures[23] = loadTexBMP("n.bmp");
	textures[24] = loadTexBMP("o.bmp");
	textures[25] = loadTexBMP("p.bmp");
	textures[26] = loadTexBMP("q.bmp");
	textures[27] = loadTexBMP("r.bmp");
	textures[28] = loadTexBMP("s.bmp");
	textures[29] = loadTexBMP("t.bmp");
	textures[30] = loadTexBMP("u.bmp");
	textures[31] = loadTexBMP("v.bmp");
	textures[32] = loadTexBMP("w.bmp");
	textures[33] = loadTexBMP("x.bmp");
	textures[34] = loadTexBMP("y.bmp");
	textures[35] = loadTexBMP("z.bmp");
}