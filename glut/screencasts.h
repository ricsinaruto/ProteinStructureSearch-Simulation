#ifndef SCREENCASTS
#define SCREENCASTS

/*standard headers*/
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


/*OpenGL and friends*/

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
#include "common.h" /* common is just defines */
#include "display.h"/* display -> setup scene to draw */
#include "draw.h"   /* draw -> draw whatever objects in the scene */
#include "error.h"  /* error convenience */
#include "fatal.h"  /* fatal convenience */
#include "interaction.h"    /* user interactions (keyboard, mouse, etc) */
#include "initialization.h" /* initialization */
#include "print.h"   /* printing functions */
#include "shapes.h"  /* basic shapes (cube, cone, etc) */
#include "textures.h"/* texture functionality */


/*  Structs  */
#include "structs.h" /* common structs */

#include "glWindowPos.h"	//vmi plusz cucc ami kell nem tom már mihez


/*GLOBALS*/
extern int screencastID;
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
extern int lightTh;
extern int ter;

/* textures */
extern unsigned int textures[44]; /* holds our textures */
extern int currentTexture;       /* current texture */

								 /*  Animation  */
extern int toggleAnimation;  /* toggle animation */


							 /*  Objects  */
extern cube_s cubes[60000]; /* cube objects */
extern molekula dronpa[38][38][38]; //molekula class 36*36*36os elrendezés

//my own shit
extern char *mouseBtnPressed;	//left vagy right
extern char *mouseState;		//up vagy down
extern char *Shift;				//forgás/mozgás fajtája 
extern double jobbx;			//asszem obsolete
extern int mouseX, mouseY;		//a q,w,e-hez kell
extern int xcoord, ycoord;		//a q,w,e-hez kell
extern int th2, ph2;			//forgatáshoz kell
extern int lightTh2, lightPh2;	//fény forgatáshoz kell
extern double fely, timed;		//asszem obsolete
extern double ecX2, ecY2;		//x és y menti mozgatáshoz kell
extern int main_window;			//obsolete
extern int szamok[17];			//tartalmazza a beírt molekula koordinátáit, ezt használom a kimenetek megadására is
extern int szamlalo;			//számolja, hogy hány szám billentyû lett lenyomva, ha értéke 3 nullázódik
extern char *enter;				//enter lenyomva
extern char *valto;				//ez váltja az entert
extern int kimenetek_szama;
extern int bemenetek_szama;
extern int kimenetek[7][16];
extern int h1[46656], h2[46656], h3[46656];


#endif