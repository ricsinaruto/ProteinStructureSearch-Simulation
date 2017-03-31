#include "screencasts.h"
/*  Global descriptions on screencasts.h  */

/*  WINDOW  */
char *windowName =	"Placeholder";
int windowHeight =	DEF_WINDOW_HEIGHT;
int windowWidth =	DEF_WINDOW_WIDTH;

/*  TOGGLE DRAW DISPLAYS  */
int toggleAxes =	DEF_AXES;
int toggleParams =	DEF_PARAMS;

/*  PROJECTION  */
double asp =	DEF_ASP;
double dim =	DEF_DIM;
int	   th  =	DEF_TH;
int    ph  =	DEF_PH;
int    fov =	DEF_FOV;
double ecX =	DEF_ECX;
double ecY =	DEF_ECY;
double ecZ =	DEF_ECZ;

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
unsigned int textures[TEX_NUMBER];
int currentTexture = TEX_DEFAULT;

/* EQUATION PARAMETERS */
double t = 0;
double dt = TIME_STEP;
double K =	K_CONST;
double tav = DIST_CONST;

double U = 0;
double Ce1 = CE1_CONST;
double Ce2 = CE2_CONST;
double Cp1 = CP1_CONST;
double Cp2 = CP2_CONST;

/* INTERACTION VARIABLES */
char *mouseBtnPressed = INIT_TEXT;
char *enter = DEF_ENTER;
char *valto = DEF_VALTO;
char *Shift = DEF_SHIFT;
char *mouseState = "";
char *fel= "";

int mouseX = 0, mouseY = 0;
int xcoord = 0, ycoord = 0;
int th2 = 0, ph2 = 0;
int lightTh2 = 0, lightPh2 = 0;
double ecX2 = 0, ecY2 = 0;
double jobbx = 0;

int main_window=0;
int szamlalo=0;
int struktura_szamlal = 0;

/* SEARCHING ALGORITHM */
int kimenetek_szama = DEF_KIMENETEK_SZAMA;
int bemenetek_szama = DEF_BEMENETEK_SZAMA;
int molekulaSzam = DEF_PROTEIN_NUMBER;
double tolerance	= DEF_TOLERANCE;
double max_ter = DEF_MAX_TER;
int szamok[17] = { 36,36,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

int *itomb_mol = new int[COORD_ARRAY_LENGTH];
int *jtomb_mol = new int[COORD_ARRAY_LENGTH];
int *ktomb_mol = new int[COORD_ARRAY_LENGTH];

//nem tudom mi ez a logikai függvény lol
int kimenetek[7][16] = { { 0,0,0,1,1,0,0,1,1,1,0,0,1,0,1,0 },{ 0,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1 },
						 { 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },{ 0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,1 },
						 { 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },{ 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },
						 { 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 }};

//bemenetek a proba függvényhez
int bemenetek[16][4] = { { 0,0,0,0 },{ 1,0,0,0 },{ 0,1,0,0 },{ 1,1,0,0 },{ 0,0,1,0 },{ 1,0,1,0 },
						{ 0,1,1,0 },{ 1,1,1,0 },{ 0,0,0,1 },{ 1,0,0,1 },{ 0,1,0 ,1 },{ 1,1,0,1 },
						{ 0,0,1,1 },{ 1,0,1,1 },{ 0,1,1,1 },{ 1,1,1,1 } };

/* VISUALIZATION */
cube_s cubes[CUBE_COUNT_BOUNDARY];
molekula dronpa[MAX_COORD][MAX_COORD][MAX_COORD];
init_molekula protein[DEF_PROTEIN_NUMBER];
int h1[CUBE_COUNT], h2[CUBE_COUNT], h3[CUBE_COUNT];