/*  Poor man's approximation of PI and e */
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
#define TIME_STEP	0.1
#define K_CONST		-0.07
#define DIST_CONST	343
#define CE1_CONST	0.008
#define CE2_CONST	180
#define CP1_CONST	0.037
#define CP2_CONST	153
#define DEF_DIPOL	-100

/* INTERACTION STRINGS */
#define INIT_TEXT	"Press the left mouse button before pressing the right one when using the q,w,e keys!";
#define DEF_SHIFT	"rotation"
#define DEF_ENTER	"pressed"
#define DEF_SWITCH	"parameters"

/* SEARCHING ALGORITHM PARAMS */
#define DEF_NUM_OUT				1		/// number of outputs
#define DEF_NUM_IN				2		/// number of inputs
#define COORD_ARRAY_LENGTH		1000	/// length of the coordinate arrays
#define DEF_TOLERANCE			0.1		/// tolerance of dipole values
#define DEF_MAX_FIELD			100		/// maximum magnitude of electric field
#define DEF_PROTEIN_NUMBER		1000	/// maximum number of proteins 
#define DEF_SAVE				false	/// whether to save simulation steps to .csv file

/* GENETIC ALGORITHM PARAMS */
#define ITER_NUMBER				100		///	number of iterations per each while cycle
#define DEF_POP_SIZE			1500	///	population size for the GA
#define DEF_MUT_RATE			0.001	///	mutation rate for the GA
#define DEF_MATING_POOL_COEFF	200		///	used to multiply fitness for mating pool sizes

#define OVER_FIT				5		///	how much should we adjust objective function compared to the desired dipole value
#define START_POINT				0		///	starting value of input fields
#define DEF_NU					0		///	expected value of normal distribution
#define DEF_R					0.008	///	parameter for random generation

/* VISUALIZATION */
#define CUBE_COUNT				46656	/// number of cubes: 36^3
#define CUBE_COUNT_BOUNDARY		60000	/// an upper boundary of number of cubes: more than 38^3
#define MAX_COORD				38		/// maximum number of coordinates in one dimension: 36+2 because of padding at edges