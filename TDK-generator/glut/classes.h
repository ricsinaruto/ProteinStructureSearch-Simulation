//molecule class for searching algorithm
class init_molekula {
public:
	int x, y, z; //coordinates
	
	bool ter; //does the molecule have a field applied
	bool kimenet; //kimenet-e az adott molekula

	double init_dipole;		//dipole moment after first run
	double desired[2];		//desired vector
	double *actual;			//actual dipole value
	double *szomszedok;		//neighbours of the molecule
	int	   bemenet_szam;	//which input is the molecule

	//basic constructor		
	init_molekula();

	//constructor
	void initialize_molekula(int _x, int  _y, int  _z, bool _ter, int _bemenet_szam, bool _kimenet,double *_szomszedok);

	//overload
	void initialize_molekula(int _x, int  _y, int  _z, bool _ter, int _bemenet_szam, bool _kimenet);

	//apply field to molecule
	void set_ter_mol();

	//get the dipole value
	double get_dipole();

	//set initial dipole
	void set_init_dipole();

	//set desired vector
	void set_desired();

	//set magnitude of field
	void set_ter(double terMag);

	//set actual dipole value
	void set_actual();

	//update actual dipole value
	void update_actual(int sor);

	//reset dipole moment
	void reset_dipole(double dipole);
};
