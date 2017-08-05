#ifndef MAIN
#define MAIN

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
#include <sstream>
#include <string>
#include <random>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>


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


/* project headers */
#include "defines.h"		/* define constants */
#include "display.h"		/* setup scene to draw */
#include "draw.h"			/* draw whatever objects in the scene */
#include "error.h"			/* error convenience */
#include "fatal.h"			/* fatal convenience */
#include "interaction.h"    /* user interactions and simulation functions */
#include "initialization.h" /* initialization of variables and objects */
#include "print.h"			/* printing functions */
#include "shapes.h"			/* basic shapes (cubes) */
#include "textures.h"		/* texture functionality */
#include "classes.h"		/* all classes */
#include "structs.h"		/* all structs */
#include "glWindowPos.h"	/* stuff to handle writing characters on the window */


/* GLOBALS */
extern char *windowName;
extern int windowHeight;
extern int windowWidth;

/* toggle views */
extern int toggleAxes;
extern int toggleParams;

/* view */
extern double asp;		/* aspect ratio */
extern double dim;		/* dimension of orthogonal box */
extern int th;			/* azimuth of view angle */
extern int ph;			/* elevation of view angle */
extern int fov;			/* field of view for perspective */
extern double ecX;      /* eye center position x */
extern double ecY;      /* eye center position y */
extern double ecZ;      /* eye center position z */

/* lighting */
extern int toggleLight;   /* toggle light */
extern double distance;   /* light distance */
extern int ambient;       /* ambient intensity % */
extern int diffuse;       /* diffuse intensity % */
extern int emission;      /* emission intensity % */
extern int specular;      /* specular intensity % */
extern int shininess;     /* shininess (power of two) */
extern float shinyvec[1]; /* shininess (value) */
extern float lightY;      /* elevation of light */
extern float white[];     /* the color white */
extern int lightPh;       /* light movement */
extern int lightTh;		  /* light rotation */

/* textures */
extern unsigned int textures[44];	 /* holds our textures */
extern int currentTexture;			 /* current texture */

/*  Animation  */
extern int toggleAnimation;			 /* toggle animation */

/*  Objects  */
extern cube_s cubes[CUBE_COUNT_BOUNDARY];					/* cube objects */
extern molecule dronpa[MAX_COORD][MAX_COORD][MAX_COORD];	/* molecule class forming 36*36*36 grid */
extern init_molecule protein[DEF_PROTEIN_NUMBER];			/* protein class, used for searching algo */


/* All other variables */
extern char *mouseBtnPressed;	// left or right
extern char *mouseState;		// up or down
extern char *Shift;				// type of movement or rotation
extern int mouseX, mouseY;		// needed for q,w,e keys
extern int xcoord, ycoord;		// needed for q,w,e keys
extern int th2, ph2;			// needed for rotation (angles)
extern int lightTh2, lightPh2;	// needed for lighting (angles)

extern double ecX2, ecY2;		// needed for movement along x and y axis
extern int main_window;			// holds the OpenGL window
extern int numbers[17];			// contatins the coordinate of the specified molecule, I also use this for the user to specify output bits
extern int counter;				// counts how many numbers have been pressed, if it reaches 3 it goes to 0
extern char *enter;				// is enter pressed?
extern char *switcher;			// switching the enter key
extern int num_out;				// number of outputs
extern int num_in;				// number of inputs
extern int outputs[7][16];		// contains the bits for the specified outputs
extern int inputs[16][4];		// contains the bits for the specified inputs

// variables for the searching algo
extern double tolerance;
extern double max_field;
extern int num_molecules;

// parameters for the equations
extern double dt;
extern double t;
extern double K;
extern double dist;
extern double U;
extern double Ce1, Ce2;
extern double Cp1, Cp2;

// store all coordinates in 3 base 36 digits
extern int base36_1[46656], base36_2[46656], base36_3[46656];

// store the coordinates of created molecules
extern int *i_array_mol, *j_array_mol, *k_array_mol;

#endif