#ifndef SCREENCASTS
#define SCREENCASTS

/* standard headers */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <string>


/* OpenGL and friends */
#ifdef USEGLEW
#include <gl/glew.h>
#endif
#define GL_GLEXT_PROTOTYPES
#ifdef __APPLE__
#include <GLUT/freeglut.h>
#else
#include <GL/freeglut.h>
#endif


/* includes */
#include "common.h"			/* common is just defines */
#include "display.h"		/* display -> setup scene to draw */
#include "draw.h"			/* draw -> draw whatever objects in the scene */
#include "error.h"			/* error convenience */
#include "fatal.h"			/* fatal convenience */
#include "interaction.h"    /* user interactions (keyboard, mouse, etc) */
#include "initialization.h" /* initialization */
#include "print.h"			/* printing functions */
#include "shapes.h"			/* basic shapes (cube, cone, etc) */
#include "textures.h"		/* texture functionality */
#include "glWindowPos.h"	/* needed for printinc characters */
#include "structs.h"		/* common structs */
	

/* GLOBALS */
extern int screencastID;
extern char *windowName;
extern int windowHeight;
extern int windowWidth;

/* toggle views */
extern int toggleAxes;
extern int toggleParams;

/* view */
extern double asp;			/* aspect ratio */
extern double dim;			/* dimension of orthogonal box */
extern int th;				/* azimuth of view angle */
extern int ph;				/* elevation of view angle */
extern int fov;				/* field of view for perspective */
extern double ecX;			/* eye center position x */
extern double ecY;			/* eye center position y */
extern double ecZ;			/* eye center position z */

/* lighting */
extern int toggleLight;		/* toggle light */
extern double distance;		/* light distance */
extern int ambient;			/* ambient intensity % */
extern int diffuse;			/* diffuse intensity % */
extern int emission;		/* emission intensity % */
extern int specular;		/* specular intensity % */
extern int shininess;		/* shininess (power of two) */
extern float shinyvec[1];	/* shininess (value) */
extern float lightY;		/* elevation of light */
extern float white[];		/* the color white */
extern int lightPh;			/* light movement */
extern int lightTh;			/* light rotation */

/* Textures */
extern unsigned int textures[44];	/* holds our textures */
extern int currentTexture;			/* current texture */

/*  Animation  */
extern int toggleAnimation;			/* toggle animation */

/*  Objects  */
extern cube_s cubes[60000];			/* cube objects */
extern molekula dronpa[38][38][38]; /* molecule class forming 36*36*36 grid */

// my own stuff
extern char *mouseBtnPressed;	// left or right
extern char *mouseState;		// up or down
extern char *Shift;				// type of movement or rotation
extern double jobbx;			// obsolete
extern int mouseX, mouseY;		// needed for q,w,e keys
extern int xcoord, ycoord;		// needed for q,w,e keys
extern int th2, ph2;			// needed for rotation
extern int lightTh2, lightPh2;	// needed for lighting

extern double fely, timed;		// maybe obsolete
extern double ecX2, ecY2;		// needed for movement along x and y axis
extern int main_window;			// obsolete
extern int szamok[17];			// contatins the coordinate of the specified molecule, I also use this for the user to specify output bits
extern int szamlalo;			// counts how many numbers have been pressed, if it reaches 3 it goes to 0
extern char *enter;				// is enter pressed?
extern char *valto;				// switching the enter key
extern int kimenetek_szama;		// number of outputs
extern int bemenetek_szama;		// number of inputs
extern int kimenetek[7][16];	// contains the bits for the specified outputs
extern int ter;					// field
extern int h1[46656], h2[46656], h3[46656];
			

#endif