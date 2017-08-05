#include "main.h"
/*  Global descriptions on main.h  */

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
double dist = DIST_CONST;

double U = 0;
double Ce1 = CE1_CONST;
double Ce2 = CE2_CONST;
double Cp1 = CP1_CONST;
double Cp2 = CP2_CONST;

/* INTERACTION VARIABLES */
char *mouseBtnPressed = INIT_TEXT;
char *enter = DEF_ENTER;
char *switcher = DEF_SWITCH;
char *Shift = DEF_SHIFT;
char *mouseState = "";

int mouseX = 0, mouseY = 0;
int xcoord = 0, ycoord = 0;
int th2 = 0, ph2 = 0;
int lightTh2 = 0, lightPh2 = 0;
double ecX2 = 0, ecY2 = 0;

int main_window=0;
int counter=0;

/* SEARCHING ALGORITHM */
int num_out = DEF_NUM_OUT;
int num_in = DEF_NUM_IN;
int num_molecules = DEF_PROTEIN_NUMBER+1;
double tolerance	= DEF_TOLERANCE;
double max_field = DEF_MAX_FIELD;
int numbers[17] = { 36,36,36,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };

int *i_array_mol = new int[COORD_ARRAY_LENGTH];
int *j_array_mol = new int[COORD_ARRAY_LENGTH];
int *k_array_mol = new int[COORD_ARRAY_LENGTH];

// a logical function
int outputs[7][16] = { { 0,0,0,1,1,0,0,1,1,1,0,0,1,0,1,0 },{ 0,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1 },
						 { 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },{ 0,0,0,1,1,1,0,1,0,0,1,0,1,0,0,1 },
						 { 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },{ 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 },
						 { 0,0,1,0,1,0,1,1,0,0,1,1,1,1,1,0 }};

// inputs of logic function; these are always the same
int inputs[16][4] = { { 0,0,0,0 },{ 1,0,0,0 },{ 0,1,0,0 },{ 1,1,0,0 },{ 0,0,1,0 },{ 1,0,1,0 },
						{ 0,1,1,0 },{ 1,1,1,0 },{ 0,0,0,1 },{ 1,0,0,1 },{ 0,1,0 ,1 },{ 1,1,0,1 },
						{ 0,0,1,1 },{ 1,0,1,1 },{ 0,1,1,1 },{ 1,1,1,1 } };

/* OBJECTS */
cube_s cubes[CUBE_COUNT_BOUNDARY];
molecule dronpa[MAX_COORD][MAX_COORD][MAX_COORD];
init_molecule protein[DEF_PROTEIN_NUMBER];
int base36_1[CUBE_COUNT], base36_2[CUBE_COUNT], base36_3[CUBE_COUNT];