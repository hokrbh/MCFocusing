// Include standard C libraries
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <cstdbool>
#include <cstdarg>

// Include other files
#include "Allocate.h"
#include "Hybrid_Taus.h"
#include "MC_Focusing.h"

int main(int argc, char *argv[])
{
	/* Parse command line arguments or use the defaults */
	char run_name[STR_SIZE]; // String for the name of the run
	char config_dir[STR_SIZE]; // String for the config directory
	safe_snprintf(run_name, STR_SIZE, "default");
	safe_snprintf(config_dir, STR_SIZE, "");
	if(argc == 1)
	{
		printf("Input parameters not provided, using defaults\n");
	}
	if(argc == 2)
	{
		printf("Only run_name given, using current directory for config_dir\n");
	}
	if(argc > 1)
	{
		safe_snprintf(run_name, STR_SIZE, "%s", argv[1]);
	}
	if(argc > 2)
	{
		safe_snprintf(config_dir, STR_SIZE, "%s", argv[2]);
	}
	
	/* Parse config file */
	PROP prop; // Struct to hold simulation parameters
	char config_path[STR_SIZE]; // String for the full path of the config file
	if(strlen(config_dir) == 0)
	{
		safe_snprintf(config_path, STR_SIZE, "%s.config", run_name);
	}
	else
	{
		safe_snprintf(config_path, STR_SIZE, "%s/%s.config", config_dir, run_name);
	}
	printf("Config file: %s\n", config_path);
	prop = read_config(config_path);
	
	/* Seed random number generator */
	if(prop.seed_from_clock == true)
	{
		prop.global_seed = time(NULL);
	}
	else
	{
		prop.global_seed = 1;
	}
	srand(prop.global_seed);
	prop.global_taus_seed.z1 = rand();
	prop.global_taus_seed.z2 = rand();
	prop.global_taus_seed.z3 = rand();
	prop.global_taus_seed.z4 = rand();
	
	/* Initialize data arrays taking care that all entries are initialized to zero */
	DATA data;
	unsigned int size_pos = (unsigned int)prop.n_x*prop.n_y*sizeof(DATA_ENTRY_2D);
	unsigned int size_time = (unsigned int)prop.n_t*sizeof(DATA_ENTRY_1D);
	data.r_pos = (DATA_ENTRY_2D*)allocate( size_pos );
	memset( data.r_pos, 0, size_pos );
	data.t_pos = (DATA_ENTRY_2D*)allocate( size_pos );
	memset( data.t_pos, 0, size_pos );
	data.r_time = (DATA_ENTRY_1D*)allocate( size_time );
	memset( data.r_time, 0, size_time );
	data.t_time = (DATA_ENTRY_1D*)allocate( size_time );
	memset( data.t_time, 0, size_time );
	data.R = 0.0;
	data.T = 0.0;
	double dx = (prop.x_f-prop.x_i)/((double)prop.n_x);
	double dy = (prop.y_f-prop.y_i)/((double)prop.n_y);
	for(int i = 0; i < prop.n_x; i++)
	{
		for(int j = 0; j < prop.n_y; j++)
		{
			int count = i + prop.n_x*j;
			// Assign x and y values of the center of the bins
			data.r_pos[count].x = prop.x_i + (i+0.5)*dx;
			data.r_pos[count].y = prop.y_i + (j+0.5)*dy;
			data.t_pos[count].x = prop.x_i + (i+0.5)*dx;
			data.t_pos[count].y = prop.y_i + (j+0.5)*dy;
		}
	}
	double dt = (prop.t_f-prop.t_i)/((double)prop.n_t);
	for(int i = 0; i < prop.n_t; i++)
	{
		// Assign t values to the center of the bins
		data.r_time[i].t = prop.t_i + (i+0.5)*dt;
		data.t_time[i].t = prop.t_i + (i+0.5)*dt;
	}
	
	/* Call MC_Focusing Function */
	MC_Focusing(prop, data);
	
	printf("R=% 8.8g\t T=%8.8g\n", data.R/((double)(prop.N_photons+1)), data.T/((double)(prop.N_photons+1)));
	
	/* Write data to files */
	char data_path[STR_SIZE];
	strncpy(data_path, prop.data_dir, STR_SIZE);
	strncat(data_path,"/",1);
	strncat(data_path, run_name, strlen(run_name));
	/* Reflection spatial profile */
	char r_pos_file[STR_SIZE];
	strncpy(r_pos_file, data_path, STR_SIZE);
	strncat(r_pos_file, "_r_pos.dat", 10);
	FILE *r_pos = fopen(r_pos_file, "w");
	if( r_pos == NULL ) // File did not open correctly
	{
		fprintf(stderr, "Error opening file %s, terminating\n", r_pos_file);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(r_pos, "# (x)\t(y)\t(w)\n");
		for(int i = 0; i < prop.n_x; i++)
		{
			for(int j = 0; j < prop.n_y; j++)
			{
				int count = i + prop.n_x*j;
				fprintf(r_pos, "% 8.8g\t% 8.8g\t% 8.8g\n", data.r_pos[count].x, data.r_pos[count].y, data.r_pos[count].w/((double)prop.N_photons) );
			}
			fprintf(r_pos, "\n");
		}
	}
	/* Transmission spatial profile */
	char t_pos_file[STR_SIZE];
	strncpy(t_pos_file, data_path, STR_SIZE);
	strncat(t_pos_file, "_t_pos.dat", 10);
	FILE *t_pos = fopen(t_pos_file, "w");
	if( t_pos == NULL ) // File did not open correctly
	{
		fprintf(stderr, "Error opening file %s, terminating\n", t_pos_file);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(t_pos, "# (x)\t(y)\t(w)\n");
		for(int i = 0; i < prop.n_x; i++)
		{
			for(int j = 0; j < prop.n_y; j++)
			{
				int count = i + prop.n_x*j;
				fprintf(t_pos, "% 8.8g\t% 8.8g\t% 8.8g\n", data.t_pos[count].x, data.t_pos[count].y, data.t_pos[count].w/((double)prop.N_photons) );
			}
			fprintf(t_pos, "\n");
		}
	}
	/* Reflection temporal profile */
	char r_time_file[STR_SIZE];
	strncpy(r_time_file, data_path, STR_SIZE);
	strncat(r_time_file, "_r_time.dat", 10);
	FILE *r_time = fopen(r_time_file, "w");
	if( r_time == NULL ) // File did not open correctly
	{
		fprintf(stderr, "Error opening file %s, terminating\n", r_time_file);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(r_time, "# (t)\t(w)\n");
		for(int i = 0; i < prop.n_t; i++)
		{
			fprintf(r_time, "% 8.8g\t% 8.8g\n", data.r_time[i].t, data.r_time[i].w/((double)prop.N_photons) );
		}
	}
	/* Transmission temporal profile */
	char t_time_file[STR_SIZE];
	strncpy(t_time_file, data_path, STR_SIZE);
	strncat(t_time_file, "_t_time.dat", 10);
	FILE *t_time = fopen(t_time_file, "w");
	if( t_time == NULL ) // File did not open correctly
	{
		fprintf(stderr, "Error opening file %s, terminating\n", t_time_file);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(t_time, "# (t)\t(w)\n");
		for(int i = 0; i < prop.n_t; i++)
		{
			fprintf(t_time, "% 8.8g\t% 8.8g\n", data.t_time[i].t, data.t_time[i].w/((double)prop.N_photons) );
		}
	}
	/* Reflection and Transmission coefficients */
	char RT_file[STR_SIZE];
	strncpy(RT_file, data_path, STR_SIZE);
	strncat(RT_file, "_RT.dat", 10);
	FILE *RT = fopen(RT_file, "w");
	if( RT == NULL ) // File did not open correctly
	{
		fprintf(stderr, "Error opening file %s, terminating\n", RT_file);
		exit(EXIT_FAILURE);
	}
	else
	{
		fprintf(RT, "# (R)\t(T)\n");
		fprintf(RT, "% 8.8g\t% 8.8g\n", data.R/((double)(prop.N_photons+1)), data.T/((double)(prop.N_photons+1)));
	}
	
	/*FILE *r_time;
	FILE *t_time
	File *RT;*/
	printf("%s\n", r_pos_file);
	
	/* Free data arrays */
	free(data.r_pos);
	free(data.t_pos);
	free(data.r_time);
	free(data.t_time);
}

/* Runs the MC_Focusing method using the properties given by prop. 
   This implementation runs in serial on the CPU					*/
void MC_Focusing(PROP &prop, DATA &data)
{
	/* Start for loop to loop over all of the photons */
	for(int i = 0; i < prop.N_photons; i++)
	{
		/* Initialize photon */
		PHOTON photon;
		/* Set photons individual seed */
		photon.seed.z1 = hybrid_Taus_int(prop.global_taus_seed);
		photon.seed.z2 = hybrid_Taus_int(prop.global_taus_seed);
		photon.seed.z3 = hybrid_Taus_int(prop.global_taus_seed);
		photon.seed.z4 = hybrid_Taus_int(prop.global_taus_seed);
		/* Initialize photon position */
		set_photon_init_position(photon, prop);
		/* Initialize photon focusing parameters */
		set_photon_init_focusing(photon, prop);
		/* Initialize photon velocity */
		set_photon_init_velocity(photon, prop);
		photon.r_id = 0; // Photon starts outside the sample
		photon.scat = false; // Photon starts out unscattered
		photon.det = 0; // Photon has not been detected
		photon.t = 0.0;
		photon.w = 1.0; // Photon starts with unit weight
		
		/* Propagate photon */
		while(photon.t < prop.t_max && photon.det == 0)
		{
			/* Find the distance the photon will travel during this step */
			if( photon.s == 0.0 )
			{
				// Compute dimensionless step size
				photon.s = -log( hybrid_Taus(photon.seed) ); // eq 3.16 Wang1995
			}
			double dr;
			if(prop.layer[photon.r_id].u_s == 0.0) // No scattering
			{
				dr = prop.max_step;
			}
			else
			{
				dr = photon.s/prop.layer[photon.r_id].u_s;
			}
			/* Ensure that photon trajectories do not take too large of a step */
			if( photon.scat == 0 )
			{
				double dr_curved = prop.epsilon*traj_curvature(photon,prop);
				if( dr >= dr_curved )
				{
					dr = dr_curved;
				}
			}
			/* See if the photon will hit the boundary within step */
			/* Note that this method is derived assuming photons travel in straight
			   paths. This is not the case in the focusing case when the photons
			   travel along curved paths, but so long as the path is not severely
			   curved this won't be a problem.										*/
			bool boundary = false; // Flag that tells if photon will hit a boundary
			int new_region; // Place holder for the region id of the neighbor layer
			if( photon.v.z != 0 ) // Photon must be moving toward a boundary to collide
			{
				double dr_bound;
				if( photon.v.z < 0 ) // Photon is moving to the left
				{
					dr_bound = (prop.layer[photon.r_id].z_0-photon.r.z)/photon.v.z;
					new_region = photon.r_id - 1; // The region the photon could refract into
				}
				else // Photon is moving to the right
				{
					dr_bound = (prop.layer[photon.r_id].z_f-photon.r.z)/photon.v.z;
					new_region = photon.r_id + 1; // The region the photon could refract into
					// Check if the photon will exit the transmission side 
					if( new_region > prop.N_layers )
					{
						new_region = 0; // Transmission side region has same properties as the background region
					}
				}
				if( dr >= dr_bound )
				{
					dr = dr_bound;
					boundary = true; // Photon will reach interface
				}
			}
			
			/* Move photon */
			if(photon.scat == false) // Move along trajectories
			{
				// Euler's method (not used but included for illustration)
				/*
				photon.r.x = photon.r.x + photon.v.x*dr;
				photon.r.y = photon.r.y + photon.v.y*dr;
				photon.r.z = photon.r.z + photon.v.z*dr;
				*/
			
				// RK4 method
				REAL3 K1, K2, K3, K4, pos, vel;
				K1.x = dr*photon.v.x;
				K1.y = dr*photon.v.y;
				K1.z = dr*photon.v.z;
				pos.x = photon.r.x + 0.5*K1.x;
				pos.y = photon.r.y + 0.5*K1.y;
				pos.z = photon.r.z + 0.5*K1.z;
				vel = evaluate_velocity(pos, photon, prop);
				K2.x = dr*vel.x;
				K2.y = dr*vel.y;
				K2.z = dr*vel.z;
				pos.x = photon.r.x + 0.5*K2.x;
				pos.y = photon.r.y + 0.5*K2.y;
				pos.z = photon.r.z + 0.5*K2.z;
				vel = evaluate_velocity(pos, photon, prop);
				K3.x = dr*vel.x;
				K3.y = dr*vel.y;
				K3.z = dr*vel.z;
				pos.x = photon.r.x + 0.5*K3.x;
				pos.y = photon.r.y + 0.5*K3.y;
				pos.z = photon.r.z + 0.5*K3.z;
				vel = evaluate_velocity(pos, photon, prop);
				K4.x = dr*vel.x;
				K4.y = dr*vel.y;
				K4.z = dr*vel.z;
				// Set new positions
				photon.r.x = photon.r.x + 0.166666666667*(K1.x + 2.0*K2.x + 2.0*K3.x + K4.x);
				photon.r.y = photon.r.y + 0.166666666667*(K1.y + 2.0*K2.y + 2.0*K3.y + K4.y);
				photon.r.z = photon.r.z + 0.166666666667*(K1.z + 2.0*K2.z + 2.0*K3.z + K4.z);
			}
			else // Move in straight lines like traditional Monte Carlo
			{
				photon.r.x = photon.r.x + photon.v.x*dr;
				photon.r.y = photon.r.y + photon.v.y*dr;
				photon.r.z = photon.r.z + photon.v.z*dr;
			}
			// Decrease weight to account for absorption
			photon.w = photon.w*exp(-dr*prop.layer[photon.r_id].u_a);
			// Update dimensionless step size
			photon.s = photon.s - dr*prop.layer[photon.r_id].u_s;
			// Update photons time
			photon.t = photon.t + dr*prop.layer[photon.r_id].n/C_0;
			// If photon hit the boundary find its new direction
			
			/* Compute what happens at the boundary */
			/* The photon is assigned a probability of reflecting between 0 and 1 by Fresnel reflection and then a random number determines if the entire photon transmits or reflects */
			if( boundary == true )
			{
				double alpha_i = acos(fabs(photon.v.z)); // Incident angle to the boundary
				double R; // Probability that a photon reflects at the boundary
				
				/* If the photon is in a region of higher refractive index, total internal reflection is possible, and if the incident angle of the photon is larger than the critical angle there will be total internal reflection */
				if( prop.layer[photon.r_id].n > prop.layer[new_region].n && alpha_i >= asin(prop.layer[new_region].n/prop.layer[photon.r_id].n) )
				{
					R = 1.0;
				}
				else // Find probability that the photon will transmit
				{
					// Transmission angle
					double alpha_t = asin(prop.layer[photon.r_id].n/prop.layer[new_region].n*sin(alpha_i)); // From Snell's Law
					double cos_i = cos(alpha_i); // Cosine of the incident angle
					double cos_t = cos(alpha_t); // Cosine of the transmission angle
					double temp1 = (prop.layer[photon.r_id].n*cos_t - prop.layer[new_region].n*cos_i);
					double temp2 = (prop.layer[photon.r_id].n*cos_t + prop.layer[new_region].n*cos_i);
					double temp3 = (prop.layer[photon.r_id].n*cos_i - prop.layer[new_region].n*cos_t);
					double temp4 = (prop.layer[photon.r_id].n*cos_i + prop.layer[new_region].n*cos_t);
					R = 0.5*( (temp1*temp1)/(temp2*temp2) + (temp3*temp3)/(temp4*temp4) );
				}
				if( hybrid_Taus(photon.seed) > R ) // Photon will transmit at the boundary
				{
					// Note that the focusing code will currently only handle index matched boundaries
					photon.v.x = photon.v.x*prop.layer[photon.r_id].n/prop.layer[new_region].n;
					photon.v.y = photon.v.y*prop.layer[photon.r_id].n/prop.layer[new_region].n;
					photon.v.z = sign(photon.v.z)*cos(alpha_i);
					photon.r_id = new_region;
					if( photon.r_id == 0 ) // Photon has exited the sample
					{
						if(photon.v.z < 0.0) // Photon has exited in reflection
						{
							photon.det = 1;
						}
						if(photon.v.z > 0.0) // Photon has exited in transmission
						{
							photon.det = 2;
						}
					}
				}
				else // Photon reflects at the boundary
				{
					// Only the z-direction is reversed in a reflection
					photon.v.z = -photon.v.z;
					if(new_region == -1) // Specular reflection at the surface
					{
						// Detect photon
						photon.det = 3; // Set detection flag for specular reflection
					}
				}
			}
			
			/* Find new velocities */
			/* Check if the photon has reached a scattering center */
			if(photon.s == 0.0)
			{
				double phi = 2.0*PI*hybrid_Taus(photon.seed);
				double xi = hybrid_Taus(photon.seed);
				double c_theta = cos_theta( xi, prop.layer[photon.r_id].g);
				double s_theta = sin_theta( xi, prop.layer[photon.r_id].g);
				double c_phi = cos( phi );
				double s_phi = sin( phi );
				// Check if velocity vector lies along the z-axis so we do not divide by zero
				if( fabs( photon.v.z ) > 0.99999999 )
				{
					photon.v.x = s_theta*c_phi;
					photon.v.y = s_theta*s_phi;
					photon.v.z = sign(photon.v.z)*c_theta;
				}
				else
				{
					double temp = sqrt( 1.0 - photon.v.z*photon.v.z );
					REAL3 temp_vel; // Need a place holder so that we do not update a value while we need the old value
					temp_vel.x = s_theta/temp*(photon.v.y*s_phi-photon.v.z*photon.v.x*c_phi) + photon.v.x*c_theta;
					temp_vel.y = s_theta/temp*(-photon.v.x*s_phi-photon.v.z*photon.v.y*c_phi) + photon.v.y*c_theta;
					temp_vel.z = s_theta*temp*c_phi + photon.v.z*c_theta;
					photon.v = temp_vel; // Replace old velocity with new velocity
				}
				photon.scat = true; // Photon has now been scattered
			}
			/* If the photon was not scattered, we need to update it's velocities according to Gaussian beam propagation */
			if( photon.scat == false )
			{
				photon.v = evaluate_velocity(photon.r, photon, prop);
			}
			
			/* Play Russian Roulette to see if the photon survives */
			if( photon.w < prop.w_thresh )
			{
				if( hybrid_Taus( photon.seed ) < 0.1 ) // Photon survives
				{
					photon.w = 10.0*photon.w;
				}
				else // Photon is killed by the roulette process
				{
					photon.w = 0.0;
					photon.det = 4;
				}
			}
		}
		
		/* Process data into data array */
		/* Find spatial distribution of reflected and transmitted photons */
		/* Find x_bin and y_bin if the photon is in the ROI */
		if( photon.r.x >= prop.x_i && photon.r.x < prop.x_f && photon.r.y >= prop.y_i && photon.r.y < prop.y_f )
		{
			int x_bin = (int) floor(prop.n_x*(photon.r.x-prop.x_i)/(prop.x_f-prop.x_i));
			int y_bin = (int) floor(prop.n_y*(photon.r.y-prop.y_i)/(prop.y_f-prop.y_i));
			int count = x_bin + y_bin*prop.n_x; // Values are stored in row-major orientation
			// Add the photons weight to the appropriate array
			if(photon.det == 1 || photon.det == 3) // Detected in reflection
			{
				data.r_pos[count].w += photon.w;
			}
			else if(photon.det == 2) // Detected in transmission
			{
				data.t_pos[count].w += photon.w;
			}
		}
		
		/* Find temporal distribution of reflected and transmitted photons */
		/* Find the t_bin if the photon is in the ROI */
		if( photon.t >= prop.t_i && photon.t < prop.t_f )
		{
			int t_bin = (int) floor(prop.n_t*(photon.t-prop.t_i)/(prop.t_f-prop.t_i));
			// Add the photons weight to the appropriate array
			if(photon.det == 1 || photon.det == 3) // Detected in reflection
			{
				data.r_time[t_bin].w += photon.w;
			}
			else if(photon.det == 2) // Detected in transmission
			{
				data.t_time[t_bin].w += photon.w;
			}
		}
		
		/* Find transmission and reflection coefficients */
		if(photon.det == 1 || photon.det == 3) // Detected in reflection
		{
			data.R += photon.w;
		}
		else if(photon.det == 2) // Detected in transmission
		{
			data.T += photon.w;
		}
		printf("Photon: % 8d\t R=% 8.8g\t T=%8.8g\n", i, data.R/((double)(i+1)), data.T/((double)(i+1)));
	}
}

/* Computes the cosine of the scattering angle for a given g using the Heyney-Greenstien 
distribution for a given uniform random number xi */
double cos_theta( double xi, double g )
{
	double ans;
	if( g == 0.0 )
	{
		ans = 2*xi - 1;
	}
	else
	{
		double temp = (1.0-g*g)/(1.0-g+2.0*g*xi);
		ans = 1.0/(2*g)*( 1.0 + g*g - temp*temp );
	}
	return(ans);
}

/* Computes the sine of the scattering angle for a given g using the Heyney-Greenstien 
distribution for a given uniform random number xi */
double sin_theta( double xi, double g )
{
	double c_theta = cos_theta( xi, g );
	return( sqrt( 1.0 - c_theta*c_theta ) );
}

/* Evaluates the velocity at the point pos using Eq. 5 Hokr et al. (2015) */
REAL3 evaluate_velocity(REAL3 pos, PHOTON &photon, PROP &prop)
{
	// T(z) = 1/R(z)
	double T = (photon.z_f-pos.z)/( (photon.z_f-pos.z)*(photon.z_f-pos.z) + photon.z_R*photon.z_R );
	double temp = 1.0/sqrt( 1.0 + T*T*(pos.x*pos.x+pos.y*pos.y) );
	REAL3 vel;
	vel.x = -temp*T*pos.x;
	vel.y = -temp*T*pos.y;
	vel.z = temp;
	return( vel );
}

/* Returns (c/n)^2*min(1/\ddot{x}, 1/\ddot{y}, 1/\ddot{z}) */
double traj_curvature(PHOTON &photon, PROP &prop)
{
	// T(z) = 1/R(z)
	double T = (photon.z_f-photon.r.z)/( (photon.z_f-photon.r.z)*(photon.z_f-photon.r.z) + photon.z_R*photon.z_R );
	// T^2*(x^2+y^2)
	double temp = T*T*( photon.r.x*photon.r.x + photon.r.y*photon.r.y );
	// F(x,y,z) Eq. 7 Hokr et al. (2015)
	double F = photon.z_R*photon.z_R/( (1.0+temp)*(1.0+temp)*( (photon.r.z-photon.z_f)*(photon.r.z-photon.z_f)+photon.z_R*photon.z_R )*( (photon.r.z-photon.z_f)*(photon.r.z-photon.z_f)+photon.z_R*photon.z_R ) );
	// \ddot{z}/(c/n)^2
	double z_dot_dot = T*F*(photon.r.x*photon.r.x+photon.r.y*photon.r.y);
	// if |x| > |y| choose x for the transverse coordinate
	double max_tran;
	if( fabs(photon.r.x)>fabs(photon.r.y) )
	{
		max_tran = photon.r.x*F; // \ddot{x}/(c/n)^2
	}
	else
	{
		max_tran = photon.r.y*F; // \ddot{y}/(c/n)^2
	}
	// Check if the z-direction is bigger than the maximum transverse curvature
	double curve;
	if( fabs(z_dot_dot) > fabs(max_tran) )
	{
		curve = fabs(1.0/z_dot_dot); // Choose the z-curvature
	}
	else
	{
		curve = fabs(1.0/max_tran); // Choose the transverse curvature
	}
	return( curve );
}

/* Set initial positions based on prop for a photon */
void set_photon_init_position(PHOTON &photon, PROP &prop)
{
	REAL3 temp_pos; // Place holder to transform to focusing geometry
	/* Currently assumes a Gaussian temporal profile for the beam,
	   but will be generalized to include arbitrary distributions
	   in future releases.											*/
	// Convert FWHM into standard deviation
	double sigma_t = 0.42466090014400953*prop.beam_d_tau*C_0;
											// 1/( 2*sqrt( 2*ln(2) ) )
	/* Use Box-Muller method to get 2 normally distributed numbers
	   from 2 uniformly distributed numbers							*/
	double R = sqrt( 2.0*sigma_t*sigma_t*log( 1.0/( 1.0-hybrid_Taus(photon.seed) ) ) );
	double theta = 2.0*PI*hybrid_Taus(photon.seed);
	temp_pos.z = R*sin(theta) - prop.beam_delay*C_0;
	
	double sigma_r = 0.25*prop.beam_d_rho;
	R = sqrt( 2.0*sigma_r*sigma_r*log( 1.0/( 1.0-hybrid_Taus(photon.seed) ) ) );
	theta = 2.0*PI*hybrid_Taus(photon.seed);
	temp_pos.x = R*cos(theta);
	temp_pos.y = R*sin(theta);
	
	/* Convert photons position from collimated space to focussing space */
	double temp = (prop.focal_depth - temp_pos.z)/sqrt(temp_pos.x*temp_pos.x + temp_pos.y*temp_pos.y + prop.focal_length*prop.focal_length); // Eq. 4 from Hokr et al. (2015)
	photon.r.x = temp_pos.x*temp;
	photon.r.y = temp_pos.y*temp;
	photon.r.z = prop.focal_depth - prop.focal_length*temp;
	if( photon.r.z >= 0.0 ) // Ensure that the photon is not initialized in the medium
	{
		fprintf(stderr, "Error: there is a photon initialized in the medium\n");
	}
}

/* Set initial photon focal depth and Rayleigh length for the photon */
void set_photon_init_focusing(PHOTON &photon, PROP &prop)
{
	/* Currently this information is stored for each individual photon 
	   so that abberations due to index mismatch may be included in the
	   future.																				*/
	photon.z_f = prop.focal_depth;
	photon.z_R = 0.001*4.0*prop.beam_wavelength*prop.focal_length*prop.focal_length/(PI*prop.beam_d_rho*prop.beam_d_rho);
}

/* Set initial photon velocities based on prop for a photon */
void set_photon_init_velocity(PHOTON &photon, PROP &prop)
{
	double T = (photon.z_f-photon.r.z)/( (photon.z_f-photon.r.z)*(photon.z_f-photon.r.z) + photon.z_R*photon.z_R ); // T = 1/R
	photon.v.x = -photon.r.x*T;
	photon.v.y = -photon.r.y*T;
	photon.v.z = sqrt(1.0-photon.r.x*photon.r.x*T*T - photon.r.y*photon.r.y*T*T);
}

/* Returns the 0 if val = 0, 1 if val > 0, and -1 if val < 0 */
template <typename T> int sign (T val)
{
	return( (T(0) < val) - (val < T(0)) );
}

/* Parses the config file pointed to by config_path to obtain the
   parameters needed for the simulation 							*/
PROP read_config(char *config_path)
{
	FILE *config;
	char buff[STR_SIZE]; // Buffer string
	char line[STR_SIZE]; // Buffer string for a line read in
	PROP prop; // Struct to hold the simulation parameters
	config = fopen(config_path, "r");
	if(config!=NULL)
	{
		/* Begin parsing simulation parameters */
		read_line(config, line); // Read first line
		read_line(config, line);
		parse_after_equals(buff,line);
		if(strcmp(buff, "here") == 0) // data_dir is assumed to be in this directory if "here" appears in the config file
		{
			safe_snprintf(prop.data_dir, STR_SIZE, "");
		}
		else
		{
			safe_snprintf(prop.data_dir, STR_SIZE, "%s", buff);
		}
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.x_i = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.x_f = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.n_x = atoi(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.y_i = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.y_f = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.n_y = atoi(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.t_i = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.t_f = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.n_t = atoi(buff);
		
		read_line(config, line);
		read_line(config, line);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.t_max = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		if(strncmp(buff,"true",4)==0)
		{
			prop.seed_from_clock = true;
			printf("Seeding random number generation from system clock\n");
		}
		else if(strncmp(buff,"false",5)==0)
		{
			prop.seed_from_clock = false;
			printf("Using predefined seed for random number generation\n");
		}
		else
		{
			fprintf(stderr,"Error: Seed from clock must be true or false\n");
			exit(EXIT_FAILURE);
		}
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.N_photons = atoi(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.N_layers = atoi(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.beam_d_rho = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.beam_d_tau = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.beam_delay = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.focal_depth = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.focal_length = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.beam_wavelength = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.max_step = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.epsilon = atof(buff);
		
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.w_thresh = atof(buff);
		
		read_line(config, line);
		read_line(config, line);
		
		/* Begin parsing the background layer properties */
		read_line(config, line);
		parse_after_equals(buff,line);
		prop.layer[0].n = atof(buff);
		
		read_line(config, line);
		
		/* Begin parsing the layer properties */
		for(int i = 1; i <= prop.N_layers; i++)
		{
			read_line(config, line);
			
			read_line(config, line);
			parse_after_equals(buff,line);
			prop.layer[i].z_0 = atof(buff);
			
			read_line(config, line);
			parse_after_equals(buff,line);
			prop.layer[i].z_f = atof(buff);
			
			read_line(config, line);
			parse_after_equals(buff,line);
			prop.layer[i].n = atof(buff);
			
			read_line(config, line);
			parse_after_equals(buff,line);
			prop.layer[i].g = atof(buff);
			
			read_line(config, line);
			parse_after_equals(buff,line);
			prop.layer[i].u_s = atof(buff);
			
			read_line(config, line);
			parse_after_equals(buff,line);
			prop.layer[i].u_a = atof(buff);
		}
		
		// Set hardcoded values for the background layer
		prop.layer[0].z_f = prop.layer[1].z_0;
	}
	return(prop);
}

/* Reads a single line from the file pointed to by file and returns a string */
void read_line(FILE *file, char *line)
{
	char buff[STR_SIZE]; // Buffer for the string being constructed
	int ch; // Int type to allow for EOF from fgetc
	bool flag = false;
	int count = 0;
	clear_buff(buff); // Clear buff to make sure that the memory is clean
	clear_buff(line); // Make sure that line is an empty string before starting
	while(flag == false)
	{
		ch = fgetc(file); // Get a character
		if( ch != EOF && ch != '\n' )
		{
			buff[count] = (char) ch;
			count ++;
		}
		else
		{
			flag = true; // fgetc has returned the line break character
		}
		if(count >= STR_SIZE - 1)
		{
			fprintf(stderr, "String buffer overflow in read_line()\n");
		}
	}
	strcpy(line,buff);
}

void safe_snprintf(char *s, int n, const char *format, ...)
{
	int err;
	va_list args;
	va_start(args, format);
	err = vsnprintf(s, n, format, args);
	if(err < 0)
	{
		fprintf(stderr, "Error formatting string\n");
	}
	if(err >= n)
	{
		fprintf(stderr, "Buffer overflow error\n");
	}
}

/* buff is a string containging what follows the first '=' of the string line */
void parse_after_equals(char *buff, const char *line)
{
	unsigned int span;
	span = strcspn(line, "=");
	clear_buff(buff);
	strcpy(buff,&line[span+1]);
}

/* clears the buffer and replaces in with an empty string */
void clear_buff(char *buff)
{
	memset(buff, 0, STR_SIZE*sizeof(char)); // Clears the entire memory area of the string
}
