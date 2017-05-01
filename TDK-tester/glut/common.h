/*  Poor man's approximation of PI */
#define PI 3.1415926535898
#define e 2.71828182845904523536
/*  Macro for sin & cos in degrees */
#define Cos(th) cos(PI/180*(th))
#define Sin(th) sin(PI/180*(th))

/*  Common #defines */
/*  Defaults for window sizing  */
#define DEF_WINDOW_HEIGHT	720
#define DEF_WINDOW_WIDTH	1280

/*  Projection  */
#define DEF_ASP		1
#define DEF_DIM		5
#define DEF_TH		340
#define DEF_PH		30
#define DEF_FOV		55
#define DEF_ECX		2
#define DEF_ECY		0
#define DEF_ECZ		4

/*  Draw defaults  */
#define DEF_AXES	1
#define DEF_PARAMS	1

/*  Shape degrees  */
#define DEF_D 5

/*  Lighting  */
#define DEF_LIGHT		1
#define DEF_DISTANCE	5
#define DEF_AMBIENT		35
#define DEF_DIFFUSE		100
#define DEF_EMISSION	0
#define DEF_SPECULAR	0
#define DEF_SHININESS	0
#define DEF_L_Y			0
#define DEF_L_PH		90
#define DEF_L_TH		90

/*  Textures  */
#define TEX_DEFAULT 36
#define TEX_NUMBER	44


/* EQUATION PARAMETERS */
#define TIME_STEP	1
#define K_CONST		-0.07
#define DIST_CONST	343
#define CE1_CONST	0.008
#define CE2_CONST	180
#define CP1_CONST	0.037
#define CP2_CONST	153
#define DEF_DIPOL	-100

/* INTERACTION CONSTANTS */
#define INIT_TEXT "Press the left mouse button before pressing the right one when using the q,w,e keys!";
#define DEF_SHIFT "rotation"
#define DEF_ENTER "pressed"
#define DEF_VALTO "parameters"

/* SEARCHING ALGORITHM */
#define DEF_KIMENETEK_SZAMA 1
#define DEF_BEMENETEK_SZAMA 2
#define COORD_ARRAY_LENGTH	1000
#define DEF_TOLERANCE		10
#define DEF_MAX_TER			100
#define DEF_PROTEIN_NUMBER	8
#define MENTES				false

/* SIMULATED ANNEALING PARAMS */
#define ITER_NUMBER			50			///number of iterations per each while cycle
#define DEF_SIGMA			0.01		///sigma a normál eloszláshoz amiből a random bemenetek választódnak
#define DEF_TEMP			40000		///valamilyen választási paraméter, ennyiszer futhat max a while ciklus
#define DEF_TEMP_CONST		DEF_TEMP/(100/n)		/// mennyivel legyen elosztva a temp
#define DEF_TEMP_BOOL		false		///legyen-e temp randomizálás
#define OVER_FIT			5			///a fitness fuctiont a desired-hez képest mennyire toljuk el
#define START_POINT			0			///a bemeneti terek kezdő értéke
#define DEF_NU				0			///normál eloszlás várható értéke
#define DEF_SUGAR			0.0008		///random generáláshoz paraméter
/*külön sigma és sugár kéne mindegyik bemenethez*/
/*adaptívan kéne változtatni a paramétereket az algoritmus futása közben*/
/*temp-el kezdeni vmit*/
/*ez a stuff-os finalbest reset is shady*/



/* VISUALIZATION */
#define CUBE_COUNT			46656
#define CUBE_COUNT_BOUNDARY	60000
#define MAX_COORD			38