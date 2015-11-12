#ifndef DATATYPES_H
#define DATATYPES_H

#define STR_SIZE			1024 				// Maximum number of characters for a string
#define MAX_LAYERS			16 					// Maximum number of layers allowed
#define PI					3.141592653589793
#define C_0					0.299792458 		// Speed of light in vacuum (mm/ps)
#define TOL					1.0E-10 			// Some small non-zero number

typedef struct
{
	double x, y, z;
} REAL3;

typedef struct
{
	/* Photon Data */
	REAL3 r;				// Position vector of photon
	REAL3 v;				// Unit vector giving the direction of the photon
	double w;				// Photon weight
	double t;				// Time
	double s;				// Demensionless step size (see Wang1995)
	TAUS_SEED seed;			// Seed for the random number generator
	/* Focusing Parameters */
	double z_f;				// Focal depth for the photon
	double z_R;				// Rayleigh length for the photon
	/* Flags */
	unsigned int r_id;		// Region photon is in
	bool scat;				// Flag for whether or not the photon has scattered
	int det;				// Flag indicating if the photon has been detected ( 0-no; 1-backwards; 2-forwards; 3-specular; 4-killed by roulette )
	
} PHOTON;

typedef struct
{
	/* Properties to define a layer */
	double z_0;			// Starting depth of layer
	double z_f;			// Ending depth of layer
	double n;			// Index of refraction of layer (Note: nonindex matched samples are not included in this version)
	double g;			// Anisotropy Parameter
	double u_s;			// Scattering coefficient (1/mm)
	double u_a;			// Absorption coefficient (1/mm)
} LAYER_PROP;

typedef struct
{
	char data_dir[STR_SIZE];		// Directory where data files are dumped to ("here" dumps to the same directory as config)
	double x_i;						// Minimum x-value for data arrays
	double x_f;						// Maximum x-value for data arrays
	int n_x;						// Number of x bins for data arrays
	double y_i;						// Minimum y-value for data arrays
	double y_f;						// Maximum y-value for data arrays
	int n_y;						// Number of y bins for data arrays
	double t_i;						// Minimum t-value for data arrays
	double t_f;						// Maximum t-value for data arrays
	int n_t;						// Number of t bins for data arrays
	double t_max;					// Maximum time to simulate photons before they are terminated
	bool seed_from_clock;			// Seed the simulation from the clock
	int N_photons;					// Number of photons to simulate
	int N_layers;					// Number of layers that make up the sample
	double beam_d_rho;				// 1/e^2 beam diameter incident on the lens in mm
	double beam_d_tau;				// FWHM temporal width of the pulse simulated in ps
	double beam_delay;				// Delay of pulse in ps
	double focal_depth;				// Depth of the focus from the surface of the sample
	double focal_length;			// Focal length of the lens used
	double beam_wavelength;			// Wavelength of the beam in microns
	LAYER_PROP layer[MAX_LAYERS];	// Properties of the layers
	double max_step;				// Maximum single step a photon can take (in mm)
	double epsilon;					// Parameter effects the step size taken for curved paths
	double w_thresh;				// Weight when Russian Roulette procedure kicks in
	int global_seed;				// Global seed used in the run
	TAUS_SEED global_taus_seed;		// Seeds used for the hybrid Taus RNG
} PROP;

typedef struct
{
	double x;		// center x-coordinate for pixel
	double y;		// center y-coordinate for pixel
	double w;		// weight for pixel
} DATA_ENTRY_2D;

typedef struct
{
	double t;		// center t-coordinate for value
	double w;		// weight of value
} DATA_ENTRY_1D;

typedef struct
{
	DATA_ENTRY_2D *r_pos; 		// Spatial distribution of reflected photons
	DATA_ENTRY_2D *t_pos; 		// Spatial distribution of transmitted photons
	DATA_ENTRY_1D *r_time; 		// Temporal distribution of reflected photons
	DATA_ENTRY_1D *t_time; 		// Temporal distribution of transmitted photons
	double R; 					// Total reflection coefficient
	double T; 					// Total transmission coefficient
} DATA;

/* Function definitions */
void MC_Focusing(PROP &prop, DATA &data);
double cos_theta( double xi, double g );
double sin_theta( double xi, double g );
REAL3 evaluate_velocity(REAL3 pos, PHOTON &photon, PROP &prop);
double traj_curvature(PHOTON &photon, PROP &prop);
void set_photon_init_position(PHOTON &photon, PROP &prop);
void set_photon_init_focusing(PHOTON &photon, PROP &prop);
void set_photon_init_velocity(PHOTON &photon, PROP &prop);
template <typename T> int sign (T val);
PROP read_config(char *config_path);
void read_line(FILE *file, char *line);
void safe_snprintf(char *s, int n, const char *format, ...);
void parse_after_equals(char *buff, const char *line);
void clear_buff(char *buff);

#endif
