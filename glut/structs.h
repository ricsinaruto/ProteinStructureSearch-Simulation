typedef struct point {
	float x;
	float y;
	float z;
} point;

typedef struct tsr {
	point t; /* translation */
	point s; /* scale */
	point r; /* rotation */
	point hely;
} tsr;

typedef struct cube_s {
	tsr tsr;
} cube_s;
class molekula
{
public:
	bool van;	//van-e molekula az adott helyen
	bool ter;	//t�rrel terhelt-e a molekula az adott helyen
	double dipA, dipB, dip;		//dip�lok
	double qeA, qeB, qp1A, qp1B, qp2A, qp2B;	//�j egyenletekhez (equations.docx)
};
