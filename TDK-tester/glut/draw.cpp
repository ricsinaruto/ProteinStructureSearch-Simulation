#include "screencasts.h"


//alapkockák megrajzolása
void drawAxes(void) {
	if (toggleAxes)	negyzetracs(0, 0, 0, 1, 1, 1, 0, 0, 0, 0);
}


//Draws a cube
void drawCube(cube_s cb) {
	//A 0,0,0 molekula nem használható!!! (nem fog kirajzolódni)
	if (cb.tsr.hely.x >0 || cb.tsr.hely.y>0 || cb.tsr.hely.z>0)
		cube(cb.tsr.t.x, cb.tsr.t.y, cb.tsr.t.z, cb.tsr.s.x, cb.tsr.s.y, cb.tsr.s.z,
			 cb.tsr.r.y, cb.tsr.hely.x, cb.tsr.hely.y, cb.tsr.hely.z);
}


//Draws the light
void drawLight(void)
{
	/*  Light switch */
	if (toggleLight) {
		/*  Translate intensity to color vectors */
		GLfloat Ambient[] = { 0.01*ambient ,0.01*ambient ,0.01*ambient ,1.0 };
		GLfloat Diffuse[] = { 0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,1.0 };
		GLfloat Specular[] = { 0.01*specular,0.01*specular,0.01*specular,1.0 };
		GLfloat Position[] = { distance*Sin(lightPh),lightY,distance*Cos(lightPh),1.0 };

		/*  Draw light position as sphere (still no lighting here) */
		glColor3fv(white);
		glDisable(GL_LIGHTING);
	

		/*  Set ambient, diffuse, specular components and position of light 0 */
		glLightfv(GL_LIGHT0, GL_AMBIENT, Ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, Diffuse);
		glLightfv(GL_LIGHT0, GL_POSITION, Position);
		glLightfv(GL_LIGHT0, GL_SPECULAR, Specular);

		/* some other params */
		glEnable(GL_COLOR_MATERIAL);
		glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
		glEnable(GL_NORMALIZE);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
	}
	else glDisable(GL_LIGHTING);
}


//Draw the parameters in the lower left corner
void drawParameters(void) {
	if (toggleParams) {
		glColor3fv(white);

		/*  Print on the screen some parameters */
		printAt(5, 5, "Angle=%d,%d  Dim=%.1f FOV=%d Coordinates=%i %i %i Light=%s",th, ph, dim, 
				fov, szamok[0],szamok[1],szamok[2], toggleLight ? "On" : "Off");
		
		//output and input params
		printAt(5, 340, "Output: ");
		int i, j = 0;
		for (i = 0; i < kimenetek_szama; i++) {
			for (j = 0; j < pow(2, bemenetek_szama); j++) 
				printAt(50+j*10, 340 - i * 20, "%i", kimenetek[i][j]);
		}
		printAt(5, 340 - (i + 1) * 20, "Number of inputs: %i", bemenetek_szama);
		printAt(5, 340 - (i + 2) * 20, "Number of outputs: %i", kimenetek_szama);

		
		//some other params
		printAt(5, 45, "Diffuse=%d  //  field=%i", diffuse,0);
		printAt(5, 25, "enter=%s // movement type=%s // mousebtn=%s // switch=%s", 
				enter, Shift, mouseBtnPressed, valto);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glColor4f(1, 1, 1, 0.4);
	}
}


//Draws the scene with everything we need in it
void drawScene(void) {
	//all draw functions
	drawLight();
	drawAxes();
	drawParameters();

	if (szamlalo == 3)	initializeObjs();				//megadott kockák inicializálása
	for (int i = 0; i < 46656; i++)	drawCube(cubes[i]);	//inicializált kockák rajzolása
	
	
	//default textúra betöltése
	currentTexture = textures[TEX_DEFAULT];
}