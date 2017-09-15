#ifndef LLGSOLVER_HEADER
#define LLGSOLVER_HEADER

// THETA = 1: stability for jx>0
// THETA = -1: stability for jx<0
#define THETA -1. // takes values in [-1,1]
//#define TSST_MICROSCOPIC
# define EXCHANGE_SIGN 1. // FM
//# define EXCHANGE_SIGN -1. // AFM

extern int LLG_dim;
extern int skyrmion_transits;

extern void LLG_ini(int Nx_in, int Ny_in, double bapp_in, double dmi_in, double alpha_in, double beta_in, double K_in, double jx_in, double jy_in, double J0value_in);
extern void open_bc(double *mag);
extern void periodic_bc(double *mag);
extern void notch_bc(double *mag, double notch_x0, double notch_x1, double notch_y0);
extern void notch_zero_bc(double *mag, double notch_x0, double notch_x1, double notch_y0);
extern void zero_bc(double *mag);
extern void z_dwall_bc(double *mag);
extern void inplane_dwall_bc(double *mag);
extern void LLG_evolve(double *mag, double *mag_new, double dt);
extern void LLG_inplane_evolve(double *mag, double *mag_new, double dt);
extern void normalize(double *mag_new);
extern void reset_corners(double *mag_new);
extern void find_max(double *mag, int *ids, int sign);
extern double find_skyrmion_max(double *mag, int sign);
extern void skyrmion_position(double *mag, double *skpos, int sign);
extern void monitor_transits(double *mag, int sign);
extern void monitor_new_record(double *mag, int sign);
extern double magenergy(double *mag);

#endif
