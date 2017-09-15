//-----------------------------------------------------------------------------------
// To set OpenMP environment variable:
// $ export OMP_NUM_THREADS=6
// $ echo $OMP_NUM_THREADS
#include <stdio.h>
#include <stdlib.h> // malloc
#include <math.h> // sqrt
#include <string.h> // memcpy
//#include <time.h> // timing
#include <omp.h>
#include "LLGsolver.h"
//#include "gsl_sf.h"

// Global variables
int LLG_dim;
int skyrmion_transits;

// Global functions
void LLG_ini(int Nx_in, int Ny_in, double bapp_in, double dmi_in, double alpha_in, double beta_in, double K_in, double jx_in, double jy_in, double J0value_in);
void LLG_evolve(double *mag, double *mag_new, double dt);
void normalize(double *mag_new);
void reset_corners(double *mag_new);
void find_max(double *mag, int *skcoord, int sign);
double find_skyrmion_max(double *mag, int sign);
void monitor_transits(double *mag, int sign);
void monitor_new_record(double *mag, int sign);
double magenergy(double *mag);

void open_bc(double *mag);
void periodic_bc(double *mag);
void notch_bc(double *mag, double notch_x0, double notch_x1, double notch_y0);
void notch_zero_bc(double *mag, double notch_x0, double notch_x1, double notch_y0);
void zero_bc(double *mag);

void x_dwall_bc(double *mag);
void y_dwall_bc(double *mag);
void z_dwall_bc(double *mag);

// Internal global variables
static int Nx;
static int Ny;
static double dmi;
static double bapp;
static double alpha, alpha2;
static double beta;
static double K;
static double jx;
static double jy;
static double J0value; // Bessel function value J0(\gamma)
static double dmi_renorm;
//static double Jx[3], Jy[3]; // renormalized exchange couplings
static const double th_prev = 0.5*(THETA-1.);
static const double th = -THETA;
static const double th_next = 0.5*(THETA+1.);
#ifdef TSST_MICROSCOPIC
    static double mu;
#endif
static double grad;

// Internal functions
int ind(int i, int j);
int indx(int i, int j);
int indy(int i, int j);
int indz(int i, int j);
void add_exchange(double *mag, int i, int j, double *beff);
//void add_exchange_renormalized(double *mag, int i, int j, double *beff);
void add_anisotropy_x(double *mag, int i, int j, double *beff);
void add_anisotropy_y(double *mag, int i, int j, double *beff);
void add_anisotropy_z(double *mag, int i, int j, double *beff);
void add_DMI(double *mag, int i, int j, double *beff);
void add_DMI_DY(double *mag, int i, int j, double *beff);
void add_DMI_renorm_J0value(double *mag, int i, int j, double *beff);
//void add_DMI_renorm_test(double *mag, int i, int j, double *beff);
void cross(const double *a, const double *b, double *axb);
void dj(double *mag, int i, int j, double jx, double jy, double *result);
void sst_torque(double *mag, int i, int j, double jx, double jy, double *tsst);
void sst_torque_inplane(double *mag, int i, int j, double jx, double jy, double *tsst);
void sst_torque_debug(double *mag, int i, int j, double jx, double jy, double *tsst);
void normalize(double *mag_new);
double find_smz_max(double *mag, int sign, int isec);

int ind(int i, int j) {
	return 3*((Ny+2)*i+j);
}

int indx(int i, int j) {
	return 3*((Ny+2)*i+j);
}

int indy(int i, int j) {
	return 3*((Ny+2)*i+j)+1;
}

int indz(int i, int j) {
	return 3*((Ny+2)*i+j)+2;
}


void LLG_ini(int Nx_in, int Ny_in, double bapp_in, double dmi_in, double alpha_in, double beta_in, double K_in, double jx_in, double jy_in, double J0value_in) {

    Nx=Nx_in;
    Ny=Ny_in;
    LLG_dim=(Nx+2)*(Ny+2)*3;
    bapp=bapp_in;
    dmi=dmi_in;
    alpha=alpha_in;
	alpha2 = alpha*alpha;    
    beta=beta_in;
	K = K_in;    
    jx=jx_in;
    jy=jy_in;
    skyrmion_transits=0;
	grad=0.;
    J0value = J0value_in;
	dmi_renorm = dmi_in*2.*J0value/(1.+J0value);
		
#ifdef TSST_MICROSCOPIC
    mu=1.;
#endif

}


void LLG_evolve(double *mag, double *mag_new, double dt) {
    int i, j, n;
#pragma omp parallel default(none) private(i,j,n) shared(Nx,Ny,bapp,mag,mag_new,dt,alpha,jx,jy)
	{
#pragma omp for schedule(static)    
	for(i=1;i<=Nx;i++)
		for(j=1;j<=Ny;j++) {

			double beff[3] = {0., 0., bapp};
			double magxbeff[3], magxmagxbeff[3], rhs[3];
			
			add_exchange(mag, i, j, beff);		

//			add_anisotropy_x(mag, i, j, beff);			
//			add_anisotropy_y(mag, i, j, beff);				
//			add_anisotropy_z(mag, i, j, beff);
			
			//add_DMI(mag, i, j, beff);
			//add_DMI_DY(mag, i, j, beff);
			add_DMI_renorm_J0value(mag, i, j, beff);
			
			cross(&mag[ind(i,j)], beff, magxbeff);
			cross(&mag[ind(i,j)], magxbeff, magxmagxbeff);

			double tsst[3], magxtsst[3];
			
			sst_torque(mag, i, j, jx, jy, tsst);
//			sst_torque_inplane(mag, i, j, jx, jy, tsst);	// Inplane torques is on		
//			sst_torque_debug(mag, i, j, jx, jy, tsst); // debug
			
			cross(&mag[ind(i,j)],tsst,magxtsst);

			for(n=0;n<3;n++)
				rhs[n] = - magxbeff[n] - alpha*magxmagxbeff[n] + tsst[n] + alpha*magxtsst[n];

			for(n=0;n<3;n++)
				mag_new[ind(i,j)+n] = mag[ind(i,j)+n] + dt*rhs[n];
		}
    }
}



void LLG_inplane_evolve(double *mag, double *mag_new, double dt) {
    int i, j, n;
#pragma omp parallel default(none) private(i,j,n) shared(Nx,Ny,bapp,mag,mag_new,dt,alpha,alpha2,jx,jy)
	{
#pragma omp for schedule(static)    
	for(i=1;i<=Nx;i++)
		for(j=1;j<=Ny;j++) {

			double beff[3] = {0., 0., bapp};
			double magxbeff[3], magxmagxbeff[3], rhs[3];
			
			add_exchange(mag, i, j, beff);

//			add_anisotropy_x(mag, i, j, beff);			
//			add_anisotropy_y(mag, i, j, beff);			
			add_anisotropy_z(mag, i, j, beff);			
			
			//add_DMI(mag, i, j, beff);
			//add_DMI_DY(mag, i, j, beff);
			add_DMI_renorm_J0value(mag, i, j, beff);
			
			cross(&mag[ind(i,j)], beff, magxbeff);
			cross(&mag[ind(i,j)], magxbeff, magxmagxbeff);

			double tsst[3], magxtsst[3];
			sst_torque_inplane(mag, i, j, jx, jy, tsst);			
			cross(&mag[ind(i,j)],tsst,magxtsst);

			for(n=0;n<3;n++)
				rhs[n] = - magxbeff[n] - alpha*magxmagxbeff[n] + tsst[n] + alpha*magxtsst[n];

			double mx = mag[ind(i,j)];
			double my = mag[ind(i,j)+1];									
			double mz = mag[ind(i,j)+2];
			double Z = rhs[2]/(1.+alpha2*mz*mz);
			rhs[0] += -alpha*my*Z - alpha2*mx*mz*Z;
			rhs[1] += alpha*mx*Z - alpha2*my*mz*Z;
			rhs[2] = (1.+alpha2)*Z;
			
			for(n=0;n<3;n++)
				mag_new[ind(i,j)+n] = mag[ind(i,j)+n] + dt*rhs[n];
		}
    }
}


void normalize(double *mag_new) { 
    int i, j, n;
   	for(i=1;i<=Nx;i++)
   		for(j=1;j<=Ny;j++) {		
            double sum=0.;
       		for(n=0;n<3;n++)
                sum += mag_new[ind(i,j)+n]*mag_new[ind(i,j)+n];
            double normalization = 1./sqrt(sum);
          	for(n=0;n<3;n++)
                mag_new[ind(i,j)+n] *= normalization;
        }
}


// Open boundary conditions
void open_bc(double *mag) {
	int i, j, n;
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = mag[ind(i,1)+n];

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = mag[ind(i,Ny)+n];

	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(0,j)+n] = mag[ind(1,j)+n];

	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(Nx+1,j)+n] = mag[ind(Nx,j)+n];
}


// Periodic boundary conditions
void periodic_bc(double *mag) {
	int i, j, n;
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = mag[ind(i,Ny)+n];

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = mag[ind(i,1)+n];

	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(0,j)+n] = mag[ind(Nx,j)+n];

	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(Nx+1,j)+n] = mag[ind(1,j)+n];
}


// Notch boundary conditions, periodic along x, open along y
void notch_bc(double *mag, double notch_x0, double notch_x1, double notch_y0) {
	int i, j, n;

// Open bc
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = mag[ind(i,1)+n];

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = mag[ind(i,Ny)+n];

// Periodic bc
	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(0,j)+n] = mag[ind(Nx,j)+n];

	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(Nx+1,j)+n] = mag[ind(1,j)+n];

// Notch
	for(i=notch_x0+1;i<=notch_x1-1;i++)
		for(n=0;n<3;n++)
			mag[ind(i,notch_y0)+n] = mag[ind(i,notch_y0-1)+n];

	for(j=notch_y0+1;j<=Ny;j++)
		for(n=0;n<3;n++) {
			mag[ind(notch_x0,j)+n] = mag[ind(notch_x0-1,j)+n];
			mag[ind(notch_x1,j)+n] = mag[ind(notch_x1+1,j)+n];
		}

	for(n=0;n<3;n++) {
		mag[ind(notch_x0,notch_y0)+n] = 0.5*(mag[ind(notch_x0-1,notch_y0)+n]+mag[ind(notch_x0,notch_y0-1)+n]);
		mag[ind(notch_x1,notch_y0)+n] = 0.5*(mag[ind(notch_x1+1,notch_y0)+n]+mag[ind(notch_x1,notch_y0-1)+n]);
	}
}


// Notch boundary conditions, periodic along x, zero otherwise
void notch_zero_bc(double *mag, double notch_x0, double notch_x1, double notch_y0) {
	int i, j, n;

// Zero bc
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = 0.;

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = 0.;

// Periodic bc
	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(0,j)+n] = mag[ind(Nx,j)+n];

	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(Nx+1,j)+n] = mag[ind(1,j)+n];

// Notch
	for(i=notch_x0+1;i<=notch_x1-1;i++)
		for(n=0;n<3;n++)
			mag[ind(i,notch_y0)+n] = 0.;

	for(j=notch_y0;j<=Ny;j++)
		for(n=0;n<3;n++) {
			mag[ind(notch_x0,j)+n] = 0.;
			mag[ind(notch_x1,j)+n] = 0.;
		}
}


// Zero boundary conditions, periodic along x
void zero_bc(double *mag) {
	int i, j, n;

// Zero bc
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = 0.;

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = 0.;

// Periodic bc
	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(0,j)+n] = mag[ind(Nx,j)+n];

	for(j=1;j<=Ny;j++)
		for(n=0;n<3;n++)
			mag[ind(Nx+1,j)+n] = mag[ind(1,j)+n];
}


// z domain wall boundary conditions
void z_dwall_bc(double *mag) {
	int i, j, n;

// Zero bc: top and bottom
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = 0.;

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = 0.;

// dwall
	for(j=1;j<=Ny;j++) {
		mag[ind(0,j)] = 0.;
		mag[ind(0,j)+1] = 0.;
		mag[ind(0,j)+2] = 1.;				
	}

	for(j=1;j<=Ny;j++) {
		mag[ind(Nx+1,j)] = 0.;
		mag[ind(Nx+1,j)+1] = 0.;
		mag[ind(Nx+1,j)+2] = -1.;				
	}
}


// in-plane domain wall boundary conditions
void y_dwall_bc(double *mag) {
	int i, j, n;

// Zero bc: top and bottom
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = 0.;

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = 0.;

// in-plane dwall
	for(j=1;j<=Ny;j++) {
		mag[ind(0,j)] = 0.;
		mag[ind(0,j)+1] = 1.;
		mag[ind(0,j)+2] = 0.;				
	}

	for(j=1;j<=Ny;j++) {
		mag[ind(Nx+1,j)] = 0.;
		mag[ind(Nx+1,j)+1] = -1.;
		mag[ind(Nx+1,j)+2] = 0.;				
	}
}


void x_dwall_bc(double *mag) {
	int i, j, n;

// Zero bc: top and bottom
	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,0)+n] = 0.;

	for(i=1;i<=Nx;i++)
		for(n=0;n<3;n++)
			mag[ind(i,Ny+1)+n] = 0.;

// in-plane x dwall
	for(j=1;j<=Ny;j++) {
		mag[ind(0,j)] = 1.;
		mag[ind(0,j)+1] = 0.;
		mag[ind(0,j)+2] = 0.;				
	}

	for(j=1;j<=Ny;j++) {
		mag[ind(Nx+1,j)] = -1.;
		mag[ind(Nx+1,j)+1] = 0.;
		mag[ind(Nx+1,j)+2] = 0.;				
	}
}

// Add exchange interaction
void add_exchange(double *mag, int i, int j, double *beff) {
	int n;
	for(n=0;n<3;n++)
		beff[n] += EXCHANGE_SIGN*(mag[n+ind(i-1,j)] + mag[n+ind(i+1,j)] + mag[n+ind(i,j-1)] + mag[n+ind(i,j+1)]);
}


// Add anisotropy along x
void add_anisotropy_x(double *mag, int i, int j, double *beff) {
	beff[0] += 2.*K*mag[ind(i,j)]; // along x
}

// Add anisotropy along y
void add_anisotropy_y(double *mag, int i, int j, double *beff) {
	beff[1] += 2.*K*mag[ind(i,j)+1]; // along y
}

// Add anisotropy along z
void add_anisotropy_z(double *mag, int i, int j, double *beff) {
	beff[2] += 2.*K*mag[ind(i,j)+2]; // along z
}


// Add exchange interaction
/*void add_exchange_renormalized(double *mag, int i, int j, double *beff) {
	int n;
	for(n=0;n<3;n++)
		beff[n] += Jx[n]*(mag[n+ind(i-1,j)] + mag[n+ind(i+1,j)]) + Jy[n]*(mag[n+ind(i,j-1)] + mag[n+ind(i,j+1)]);
}*/


// Add Dzyaloshinskii-Moriya interaction
void add_DMI(double *mag, int i, int j, double *beff) {
	beff[0] += dmi*(mag[indz(i,j-1)]-mag[indz(i,j+1)]);
	beff[1] += dmi*(mag[indz(i+1,j)]-mag[indz(i-1,j)]);
	beff[2] += dmi*(mag[indx(i,j+1)]-mag[indy(i+1,j)]-mag[indx(i,j-1)]+mag[indy(i-1,j)]);
}

// Add Dzyaloshinskii-Moriya interaction, Dima Yudin form
void add_DMI_DY(double *mag, int i, int j, double *beff) {
	beff[0] += dmi*(mag[indz(i-1,j)]-mag[indz(i+1,j)]);
	beff[1] += dmi*(mag[indz(i,j-1)]-mag[indz(i,j+1)]);
	beff[2] += dmi*(mag[indx(i+1,j)]-mag[indy(i,j-1)]-mag[indx(i-1,j)]+mag[indy(i,j+1)]);
}

// Add Dzyaloshinskii-Moriya interaction, renormalized
/*void add_DMI_renorm_absJ0(double *mag, int i, int j, double *beff) {
    double coeff = 2.*absJ0/(1.+absJ0);
	beff[0] += coeff*dmi*(mag[indz(i-1,j)]-mag[indz(i+1,j)]);
	beff[1] += coeff*dmi*absJ0*(mag[indz(i,j-1)]-mag[indz(i,j+1)]);
	beff[2] += coeff*dmi*(mag[indx(i+1,j)]-absJ0*mag[indy(i,j-1)]-mag[indx(i-1,j)]+absJ0*mag[indy(i,j+1)]);
}*/

// Add Dzyaloshinskii-Moriya interaction, renormalized
void add_DMI_renorm_J0value(double *mag, int i, int j, double *beff) {
//    double coeff = -2.*J0value/(1.+J0value); // new skyrmion at J0value=-0.35
	beff[0] += dmi_renorm*(mag[indz(i-1,j)]-mag[indz(i+1,j)]);
	beff[1] += dmi_renorm*J0value*(mag[indz(i,j-1)]-mag[indz(i,j+1)]);
	beff[2] += dmi_renorm*(mag[indx(i+1,j)]-J0value*mag[indy(i,j-1)]-mag[indx(i-1,j)]+J0value*mag[indy(i,j+1)]);
}

// Add Dzyaloshinskii-Moriya interaction, renormalized
/*void add_DMI_renorm_test(double *mag, int i, int j, double *beff) {
    double coeff = 2.*J0value/(1.+J0value);
	beff[0] += coeff*dmi*(mag[indz(i-1,j)]-mag[indz(i+1,j)]);
	beff[1] += coeff*dmi*(-J0value)*(mag[indz(i,j-1)]-mag[indz(i,j+1)]);
	beff[2] += coeff*dmi*(mag[indx(i+1,j)]-(-J0value)*mag[indy(i,j-1)]-mag[indx(i-1,j)]+(-J0value)*mag[indy(i,j+1)]);
}*/


// Cross product
void cross(const double *a, const double *b, double *axb) {
	axb[0] = a[1]*b[2] - a[2]*b[1];
	axb[1] = a[2]*b[0] - a[0]*b[2];
	axb[2] = a[0]*b[1] - a[1]*b[0];
}

// Derivative along the direction jdir
void dj(double *mag, int i, int j, double jx, double jy, double *result) {
	int n;
	for(n=0;n<3;n++)
		result[n] = jx*(th_prev*mag[n+ind(i-1,j)] + th_next*mag[n+ind(i+1,j)]) +
				   jy*(th_prev*mag[n+ind(i,j-1)] + th_next*mag[n+ind(i,j+1)]) + 
					(jx+jy)*th*mag[n+ind(i,j)];

}

// SST torque
void sst_torque(double *mag, int i, int j, double jx, double jy, double *tsst) {
	int n;
	double djmag[3], magxdjmag[3], magxdjmag_z[3];
	double mx = mag[ind(i,j)];
	double my = mag[ind(i,j)+1];

	dj(mag, i,j,jx,jy,djmag);
	cross(&mag[ind(i,j)], djmag, magxdjmag);

	magxdjmag_z[0] = djmag[2]*my;
	magxdjmag_z[1] = -djmag[2]*mx;
	magxdjmag_z[2] = 0.;

#ifdef TSST_MICROSCOPIC	
	double mz = mag[ind(i,j)+2];
	double xifun = (1.-mz*mz)*beta/(1.+4.*mu*mu);
	for(n=0;n<3;n++)
		tsst[n] = mu*(1.+2.*xifun) * djmag[n] - (xifun+beta) * magxdjmag[n] + beta * magxdjmag_z[n];
#else
	for(n=0;n<3;n++)
		tsst[n] = djmag[n] - beta*magxdjmag[n];
#endif
}


// SST torque inplane
void sst_torque_inplane(double *mag, int i, int j, double jx, double jy, double *tsst) {
	int n;
	double djmag[3], magxdjmag[3], magxdjmag_z[3];
	double mx = mag[ind(i,j)];
	double my = mag[ind(i,j)+1];

	dj(mag, i,j,jx,jy,djmag);
	cross(&mag[ind(i,j)], djmag, magxdjmag);

	magxdjmag_z[0] = djmag[2]*my;
	magxdjmag_z[1] = -djmag[2]*mx;
	magxdjmag_z[2] = 0.;

	for(n=0;n<3;n++)
		tsst[n] = djmag[n] - beta*magxdjmag[n] + beta * magxdjmag_z[n];
}



// SST torque debug
void sst_torque_debug(double *mag, int i, int j, double jx, double jy, double *tsst) {
	int n;
	double djmag[3], magxdjmag[3], magxdjmag_z[3];
	double mx = mag[ind(i,j)];
	double my = mag[ind(i,j)+1];

	dj(mag, i,j,jx,jy,djmag);
	cross(&mag[ind(i,j)], djmag, magxdjmag);

	magxdjmag_z[0] = djmag[2]*my;
	magxdjmag_z[1] = -djmag[2]*mx;
	magxdjmag_z[2] = 0.;

	for(n=0;n<3;n++)
//		tsst[n] = djmag[n] - beta*magxdjmag[n] + beta * magxdjmag_z[n];
		tsst[n] = djmag[n] - beta * magxdjmag_z[n];
}


// Set unused virtual nodes to zero
void reset_corners(double *mag_new) {
    int n;
    for(n=0;n<3;n++) {
        mag_new[ind(0,0)+n]=0.;
        mag_new[ind(0,Ny+1)+n]=0.;
        mag_new[ind(Nx+1,0)+n]=0.;
        mag_new[ind(Nx+1,Ny+1)+n]=0.;
    }
}


void find_max(double *mag, int *ids, int sign) {
    int i, j;
    double smz_max=0.;
    ids[0]=0;
    ids[1]=0;
    for(i=1;i<=Nx;i++)
        for(j=1;j<=Ny;j++) {
            double smz = sign*mag[ind(i,j)+2];
            if(smz > 0.5 && smz > smz_max) {
                ids[0]=i;
                ids[1]=j;
                smz_max = smz;
            }
        }
}

// sign: sign(mz) in skyrmion center
double find_skyrmion_max(double *mag, int sign) {
    int i, j;
    double smz_max = -1.;
    for(i=1;i<=Nx;i++)
        for(j=1;j<=Ny;j++) {
            double smz = sign*mag[ind(i,j)+2];
            if(smz > smz_max) {
                smz_max = smz;
            }
        }
return smz_max;
}

void skyrmion_position(double *mag, double *skpos, int sign) {
	int ids[2];
	find_max(mag, ids, sign);
	int i=ids[0];
	int j=ids[1];
	double mz0 = mag[ind(i,j)+2];
	double mz1 = mag[ind(i+1,j)+2];
	double mz2 = mag[ind(i-1,j)+2];
	double mz3 = mag[ind(i,j+1)+2];
	double mz4 = mag[ind(i,j-1)+2];
	double gx = mz0 - 0.5*(mz1+mz2);
	double gy = mz0 - 0.5*(mz3+mz4);
	skpos[0] = i + (mz1-mz2)/(4.*gx);
	skpos[1] = j + (mz3-mz4)/(4.*gy);
}



// Find smz_max at fixed i=isec
double find_smz_max(double *mag, int sign, int isec) {
	int j;
	double smz_max=-1;
    for(j=1;j<=Ny;j++) {
        double smz = sign*mag[ind(isec,j)+2];
        if(smz > smz_max)
            smz_max = smz;
    }
return smz_max;
}


void monitor_new_record(double *mag, int sign) {
	double smz_max_1 = find_smz_max(mag, sign, 1);
	double smz_max_Nx = find_smz_max(mag, sign, Nx);
	grad = smz_max_1-smz_max_Nx;
}

void monitor_transits(double *mag, int sign) {
	double smz_max_1 = find_smz_max(mag, sign, 1);
	double smz_max_Nx = find_smz_max(mag, sign, Nx);
	double grad_new = smz_max_1-smz_max_Nx;

    if(smz_max_Nx > 0.5 && grad_new > 0 && grad < 0) {
        skyrmion_transits++;
//        printf("# Skyrmion passed\n");
    }

	grad = grad_new;
}


// Warning: virtual nodes are in use
// make sure virtual nodes have the correct values (run the corresponding b/c routine)
double magenergy(double *mag) {
	int i, j;
	double Henergy, sum;

	// applied field contribution	
	sum = 0.;
    for(i=1;i<=Nx;i++)
        for(j=1;j<=Ny;j++)
            sum += mag[ind(i,j)+2];
    Henergy = - bapp*sum;
//    printf("# bapp*sum: %f\n",bapp*sum);
    
	// dmi contribution    
    sum = 0.;
    for(i=1;i<=Nx;i++)
        for(j=1;j<=Ny;j++) {
	        double mx = mag[ind(i,j)];
			double my = mag[ind(i,j)+1];	        
        	double mz = mag[ind(i,j)+2];    
        	sum += mz*mag[ind(i+1,j)] - mx*mag[ind(i+1,j)+2] + J0value*(mz*mag[ind(i,j+1)+1] - my*mag[ind(i,j+1)+2]);
        	}
    Henergy -= dmi_renorm*sum;
//    printf("# dmi_renorm*sum: %f\n",dmi_renorm*sum);    
    
   	// exchange contribution
    sum = 0.;
    for(i=1;i<=Nx;i++)
        for(j=1;j<=Ny;j++) {
	        double mx = mag[ind(i,j)];
			double my = mag[ind(i,j)+1];	        
        	double mz = mag[ind(i,j)+2];    
			sum += mx*(mag[ind(i+1,j)]+mag[ind(i,j+1)]);
			sum += my*(mag[ind(i+1,j)+1]+mag[ind(i,j+1)+1]);
			sum += mz*(mag[ind(i+1,j)+2]+mag[ind(i,j+1)+2]);
        }
    Henergy -= EXCHANGE_SIGN*sum;
//    printf("# sum: %f\n",sum);

   	// anisotropy contribution (along z)
    sum = 0.;
    for(i=1;i<=Nx;i++)
        for(j=1;j<=Ny;j++)
   	        sum += mag[ind(i,j)+2]*mag[ind(i,j)+2];
	Henergy -= K*sum;
	
return Henergy;           
}
