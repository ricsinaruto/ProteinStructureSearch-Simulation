#include "main.h"


// Draw the basic cubes
void drawAxes(void) {
	if (toggleAxes)	grid(0, 0, 0, 1, 1, 1, 0, 0, 0, 0);
}

// Draws a cube
void drawCube(cube_s cb) {
	// the molecule at coordinates 0;0;0 can not be used!!! (it won't be drawn)
	if (cb.tsr.pos.x >0 || cb.tsr.pos.y>0 || cb.tsr.pos.z>0)
		cube(cb.tsr.t.x, cb.tsr.t.y, cb.tsr.t.z, cb.tsr.s.x, cb.tsr.s.y, cb.tsr.s.z,
			 cb.tsr.r.y, cb.tsr.pos.x, cb.tsr.pos.y, cb.tsr.pos.z);
}

// Draws the light
void drawLight(void)
{
	/* Light switch */
	if (toggleLight) {
		/* Translate intensity to color vectors */
		GLfloat Ambient[] = { 0.01*ambient ,0.01*ambient ,0.01*ambient ,1.0 };
		GLfloat Diffuse[] = { 0.01*diffuse ,0.01*diffuse ,0.01*diffuse ,1.0 };
		GLfloat Specular[] = { 0.01*specular,0.01*specular,0.01*specular,1.0 };
		GLfloat Position[] = { distance*Sin(lightPh),lightY,distance*Cos(lightPh),1.0 };

		/* Draw light position as sphere (still no lighting here) */
		glColor3fv(white);
		glDisable(GL_LIGHTING);
	
		/* Set ambient, diffuse, specular components and position of light 0 */
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

// Draw the parameters in the lower left corner
void drawParameters(void) {
	// Parameter printing switch
	if (toggleParams) {
		glColor3fv(white);

		/* Print on the screen some parameters */
		printAt(5, 5, "Angle=%d,%d  Dim=%.1f FOV=%d Coordinates=%i %i %i Light=%s",th, ph, dim, 
				fov, numbers[0],numbers[1],numbers[2], toggleLight ? "On" : "Off");
		
		// output and input params
		printAt(5, 340, "Output: ");
		int i, j = 0;
		for (i = 0; i < num_out; i++) {
			for (j = 0; j < pow(2, num_in); j++) 
				printAt(50+j*10, 340 - i * 20, "%i", outputs[i][j]);
		}
		printAt(5, 340 - (i + 1) * 20, "Number of inputs: %i", num_in);
		printAt(5, 340 - (i + 2) * 20, "Number of outputs: %i", num_out);

		
		// some other params
		printAt(5, 45, "Diffuse=%d  //  field=%i", diffuse,0);
		printAt(5, 25, "enter=%s // movement type=%s // mousebtn=%s // switch=%s", 
				enter, Shift, mouseBtnPressed,switcher);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glColor4f(1, 1, 1, 0.4);
	}
}

// Draw the scene with everything we need in it
void drawScene(void) {
	//all draw functions
	drawLight();
	drawAxes();
	drawParameters();

	if (counter == 3)	initializeObjs();				// initialize the created cubes
	for (int i = 0; i < 46656; i++)	drawCube(cubes[i]);	// draw the initialized cubes
	
	//load the default texture
	currentTexture = textures[TEX_DEFAULT];
}