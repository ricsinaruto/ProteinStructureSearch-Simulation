#include "screencasts.h"
/*  Global descriptions on screencasts.h  */

/*  ID-used to keep screencasts separate  */
int screencastID = 0;

/*  WINDOW  */
char *windowName = "OpenGL screenscasts XX: Placeholder";
int windowHeight = DEF_WINDOW_HEIGHT;
int windowWidth = DEF_WINDOW_WIDTH;

/*  TOGGLE DRAW DISPLAYS  */
int toggleAxes = DEF_AXES;
int toggleParams = DEF_PARAMS;

/*  PROJECTION  */
double asp = DEF_ASP;
double dim = DEF_DIM;
int th = DEF_TH;
int ph = DEF_PH;
int fov = DEF_FOV;
double ecX = DEF_ECX;
double ecY = DEF_ECY;
double ecZ = DEF_ECZ;

/*  LIGHTING  */
int toggleLight = DEF_LIGHT;
double distance = DEF_DISTANCE;
int ambient = DEF_AMBIENT;
int diffuse = DEF_DIFFUSE;
int emission = DEF_EMISSION;
int specular = DEF_SPECULAR;
int shininess = DEF_SHININESS;
float shinyvec[1] = { 1 };
float lightY = DEF_L_Y;
float white[] = { 1,1,1,1 };
int lightPh = DEF_L_PH;
int lightTh = DEF_L_TH;

/*  TEXTURES  */
unsigned int textures[44];
int currentTexture = TEX_DEFAULT;

// same as in screencasts.h
char *mouseBtnPressed = "";
char *mouseState = "";
char *Shift = "";
char *fel="";
double jobbx = 0;
int mouseX = 0, mouseY = 0;
int xcoord = 0, ycoord = 0;
int th2 = 0, ph2 = 0;
double ecX2 = 0, ecY2 = 0;
int lightTh2 = 0, lightPh2 = 0;
double fely = 0, timed = 0;
int main_window;
int szamok[17] = { 36,36,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
int szamlalo=0;
int ter = 0;
char *enter="";
char *valto = "";
int kimenetek_szama = 6;
int bemenetek_szama = 3;
int kimenetek[7][16] = { { 0,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0 },{ 0,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1 },
{ 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },{ 0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,1 },{ 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },{ 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },{ 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 }};
int h1[46656], h2[46656], h3[46656];


/*  OBJECTS  */
cube_s cubes[60000];
molekula dronpa[38][38][38];
