#ifndef MCFL_3DF_LIB_H_
#define MCFL_3DF_LIB_H_

#endif /*MCFLLIB_H_*/


#define Boolean char


/* DECLARE FUNCTIONS */
void mcxyz(double Nphotons, float dx, float dy, float dz, int Nx, int Ny, int Nz, float xs, float ys, float zs, float ux0, float uy0, float uz0, float mod_freq,
			int mcflag, int launchflag, int boundaryflag, int pfuncflag, int outflag, int printprogress, float *n_in, float *n_out, float xfocus, float yfocus, float zfocus,	float radius, float	waist, 
		   	float *muav, float *musv, float *gv, float *gammav, float *mua_multv, float *mus_multv, float *muaxfv, float *muaf_multv, float *qeffv, float *tauv, double *FX_real, double *FX_imag, 
			double *RX_real, double *RX_imag, double *FM_real, double *FM_imag, double *RM_real, double *RM_imag, char *vX, char *vF, char *xface);
double RandomGen(char Type, long Seed, long *Status);  
/* Random number generator */
Boolean SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2, double dx,double dy,double dz);
/* Asks,"In the same voxel?" */
double max2(double a, double b);
double min2(double a, double b);
double min3(double a, double b, double c);
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz);
double FindVoxelFace2(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz);
/* How much step size will the photon take to get the first voxel crossing in one single long step? */
double RFresnel(float n1, float n2, double ca1, double *ca2_Ptr);
void SortScatEvents(int nofevents,int *outarray, double xcoord, double ycoord);
