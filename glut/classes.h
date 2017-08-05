/* molecule class for searching algorithm */
class init_molecule {
public:
	int x, y, z;			// coordinates
				
	bool hasField;			// does the molecule have a field applied?
	bool isOutput;			// is the molecule an output?
	bool deleted;			// is the molecule deleted?

	double init_dipole;		// dipole moment after first run
	double desired[2];		// desired vector
	double *actual;			// actual dipole value
	int	   input_num;		// which input is the molecule

	// vector to store the molecule's neighbours
	std::vector<std::vector<int>> neighbours; 


	// basic constructor		
	init_molecule(); 

	// constructor
	void initialize_molecule(int _x, int  _y, int  _z, bool _hasField, int _input_num, bool _isOutput);

	// set neighbours
	void set_neighbours();

	// delete molecule
	void delete_molecule();

	// apply field to molecule
	void set_field_mol();

	// unset field on molecule
	void unset_field_mol();

	// get the dipole value
	double get_dipole();

	// set initial dipole
	void set_init_dipole();

	// set desired vector
	void set_desired();

	// set magnitude of field
	void set_field(double fieldMag);

	// set actual dipole value
	void set_actual();

	// update actual dipole value
	void update_actual(int row);

	// reset dipole moment
	void reset_dipole(double dipole);
};
  

/* genetic algorithm class */
class DNA {
public:
	double *genes;		// genes represented by a double array
	double fitness;		// fitness score
	bool similar;		// is it similar?
	int vec_len;		// length of genes array
	int mol_num;		// number of molecules

	// constructor
	DNA();

	// get the fields from genes
	double **getFields();

	// calculate fitness score
	void calcFitness();

	// crossover function
	DNA crossover(DNA partner);

	// mutate function
	void mutate(float mutationRate);
};