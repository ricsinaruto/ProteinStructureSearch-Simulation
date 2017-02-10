typedef struct point {
	float x;
	float y;
	float z;
} point;

typedef struct tsr {
	point t;	/* translation */
	point s;	/* scale */
	point r;	/* rotation */
	point hely;	/*coordinate position*/
} tsr;

typedef struct cube_s {
	tsr tsr;
} cube_s;


//the main protein class, containing all protein objects
class molekula
{
public:
	bool van;	//van-e molekula az adott helyen
	bool ter;	//t�rrel terhelt-e a molekula az adott helyen
	double terMag;	//mekkora t�rrel terhelt a molekula
	double dipA, dipB, dip;		//dip�lok
	double qeA, qeB, qp1A, qp1B, qp2A, qp2B;	//�j egyenletekhez (equations.docx)
};
