typedef struct point {
	float x;
	float y;
	float z;
} point;

typedef struct tsr {
	point t;	/* translation */
	point s;	/* scale */
	point r;	/* rotation */
	point pos;	/* coordinate position */
} tsr;

typedef struct cube_s {
	tsr tsr;
} cube_s;


// the main molecule class, containing all protein objects
class molecule {
public:
	bool is;			// is there a molecule at that place?
	bool needed;		// is it used for the logic function?
	bool field;			// is there a field applied to it?
	double fieldMag;	// what is the magnitude of the field applied?
	double dipA, dipB, dip;						// dipole values for equations
	double qeA, qeB, qp1A, qp1B, qp2A, qp2B;	// for the new equations (equations.docx)
};