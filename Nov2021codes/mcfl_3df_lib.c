#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


#define ls          1.0E-7      /* Moving photon a little bit off the voxel face */
#define NANOSEC		1.0E-9
#define MHZ			1.0E6
#define	PI          3.1415926
#define	LIGHTSPEED	2.997925E10 /* in vacuo speed of light [cm/s] */
#define ALIVE       1   		/* if photon not yet terminated */
#define DEAD        0    		/* if photon is to be terminated */
#define THRESHOLDX  0.01		/* used in roulette */
#define THRESHOLDM	1e-9
#define CHANCE      0.1  		/* used in roulette */
#define Boolean     char
#define SQR(x)		(x*x) 
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   (double) RandomGen(1, 0, NULL) /* Calls for a random number. */
#define COS90D      1.0E-6          /* If cos(theta) <= COS90D, theta >= PI/2 - 1e-6 rad. */
#define ONE_MINUS_COSZERO 1.0E-12   /* If 1-cos(theta) <= ONE_MINUS_COSZERO, fabs(theta) <= 1e-6 rad. */

/* DECLARE FUNCTIONS */
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

void mcxyz(double Nphotons, float dx, float dy, float dz, int Nx, int Ny, int Nz, float xs, float ys, float zs, float ux0, float uy0, float uz0, float mod_freq,
			int mcflag, int launchflag, int boundaryflag, int pfuncflag, int outflag, int printprogress, float *n_in, float *n_out, float	xfocus, float yfocus, float zfocus,	float radius, float	waist, 
		   	float *muav, float *musv, float *gv, float *gammav, float *mua_multv, float *mus_multv, float *muaxfv, float *muaf_multv, float *qeffv, float *tauv, double *FX_real, double *FX_imag, 
			double *RX_real, double *RX_imag, double *FM_real, double *FM_imag, double *RM_real, double *RM_imag, char *vX, char *vF, char *xface){
		/* Propagation parameters */
	double	x, y, z;        /* photon position */
	double	ux, uy, uz;     /* photon trajectory as cosines */
	double  uxx, uyy, uzz;	/* temporary values used during SPIN */
	double	s;              /* step sizes. s = -log(RND)/mus [cm] */
	double  sleft;          /* dimensionless */
	double	costheta;       /* cos(theta) */
	double  sintheta;       /* sin(theta) */
	double	cospsi;         /* cos(psi) */
	double  sinpsi;         /* sin(psi) */
	double	psi;            /* azimuthal angle */
	long	i_photon;       /* current photon */
	double	W, Wphase, Wtemp;              /* photon weight */
	double	absorb;         /* weighted deposited in a step due to absorption */
	short   photon_status;  /* flag = ALIVE=1 or DEAD=0 */
	Boolean sv;             /* Are they in the same voxel? */
	double 	Ref;			/*Reflection coefficient */
	double  uz1;			/*transmit angle */
	
	
	/* other variables */
	double	mua, muaxf;            /* absorption coefficient [cm^-1] intrinsic, fluorophore*/
	double	mus;            /* scattering coefficient [cm^-1] */
	double	g;              /* anisotropy [-] */
	double 	gamma;
	double 	ang_freq;
	double 	omega, omega_by_c0;
	double	mua_mult, mus_mult, muaf_mult; 		/* multiplication factor for emission wavelength */
	double	qeff;			/* fluorophore quantum efficiency */
	double 	tau; 			/* fluorescence lifetime */
	
	
	/* dummy variables */
	double  rnd;            /* assigned random value 0-1 */
	double	r, phi;			/* dummy values */
	long	i,j,NN,Nxy;         /* dummy indices */
	long 	ri;
	double	tempx, tempy, tempz; /* temporary variables, used during photon step. */
	int 	ix, iy, iz;     /* Added. Used to track photons */
	double 	temp;           /* dummy variable */
    int     bflag;          /* boundary flag:  0 = photon inside volume. 1 = outside volume */
    int 	update_flag; 
	int		CNT;
	int 	type, ftype; 	/* tissue type and fluorophore type*/
	int 	xtype;			/* boundary type */
	double	alpha, num_alpha, denom_alpha;
	double 	P; 
	int		emflag; 
    //int 	fo_index;
    int  	prev_det;
    double 	temp_real, temp_imag; /* temporary variable for storing fluence */
    double 	un_in, un_out;  /* normal component at boundary */
    int 	rx_or_tx; /*reflect or transmit */
    
    /* time */
	float	time_min;               // Requested time duration of computation.
	time_t	now;
	double	start_time, finish_time, temp_time; /* for clock() */
	
	start_time = clock();
	now = time(NULL);
//	printf("%s\n", ctime(&now));	
	
	/**** INITIALIZATIONS 
	 *****/

	RandomGen(0, -(int)time(NULL)%(1<<15), NULL); /* initiate with seed = 1, or any long integer. */
	NN =Nx*Ny*Nz;
	Nxy = Ny*Nx;
	omega = 2*PI*mod_freq*MHZ;
	omega_by_c0 = omega/LIGHTSPEED;

	for(j=0; j<NN;j++) 	{
					FX_real[j] = 0; FX_imag[j] = 0; 
					FM_real[j] = 0; FM_imag[j] = 0; 
					//Forigin[j]=0;
				} // ensure F[] starts empty.	
	for(j=0; j<Nxy;j++){ RX_real[j] = 0;	RX_imag[j]=0;
					RM_real[j] = 0;	RM_imag[j]=0;} // ensure F[] starts empty.	
	
	/**** RUN
	 Launch N photons, initializing each one before progation.
	 *****/
	time_min =10;
	i_photon = 0;
	CNT = 0;
	do {
		/**** LAUNCH 
		 Initialize photon position and trajectory.
		 *****/
		//if (fmod(i_photon,10)==0) printf("photon %ld took %d steps\n",i_photon,CNT);

		i_photon += 1;				/* increment photon count */
		W = 1.0;                    /* set photon weight to one */
		Wphase =0.0;				/* set photon phase */ 
		photon_status = ALIVE;      /* Launch an ALIVE photon */
		CNT = 0;
		emflag = 0; // set to 0 if current photon is an excitation photon, 1 = emission photon
		
		// Print out message about progress.
		
		if (printprogress==1){
		if ((i_photon>1000) & (fmod(i_photon, (int)(Nphotons/100))  == 0)) {
            temp = i_photon/Nphotons*100;
            //printf("%0.1f%% \t\tfmod = %0.3f\n", temp,fmod(temp, 10.0));
            if ((temp<10) | (temp>90)){
                printf("%0.0f%% done\n", i_photon/Nphotons*100);
            }
            else if(fmod(temp, 10.0)>9)
                printf("%0.0f%% done\n", i_photon/Nphotons*100);
        }
    }
        
		// At 1000th photon, update Nphotons to achieve desired runtime (time_min)
		/*if (i_photon==1)
			temp_time = clock();
		if (i_photon==1000) {    
			finish_time = clock();
			Nphotons = (long)( time_min*60*999*CLOCKS_PER_SEC/(finish_time-temp_time) );
			printf("Nphotons = %0.0f for simulation time = %0.2f min\n",Nphotons,time_min);
		} */
	
		/**** SET SOURCE 
		 * Launch collimated beam at x,y center. 
		 ****/
        
		/****************************/
		/* Initial position. */		
		
		/* trajectory */
		// launchfglag 1-3 -- manually set launch vector through matlab script
		// launchflag = 0 - set launch vector through code.
		if (launchflag==1) {
			// launch from z = constant
			x	= xs + radius*(2*RandomGen(1,0,NULL)-1); 
			y	= ys + radius*(2*RandomGen(1,0,NULL)-1); 
			z	= zs;
			ux	= ux0;
			uy	= uy0;
			uz	= uz0;
		}
		else if (launchflag==2){
			// launch from x = constant
			x	= xs;
			y	= ys + radius*(2*RandomGen(1,0,NULL)-1); 
			z	= zs + radius*(2*RandomGen(1,0,NULL)-1); 
			ux	= ux0;
			uy	= uy0;
			uz	= uz0;
		}
		else if (launchflag==3){
			// launch from y = constant
			x	= xs + radius*(2*RandomGen(1,0,NULL)-1); 
			y	= ys;
			z	= zs + radius*(2*RandomGen(1,0,NULL)-1); 
			ux	= ux0;
			uy	= uy0;
			uz	= uz0;
		}
		else { 
		// For theses cases, currently only supports source on z = constant plane. 
		// use mcflag
			if (mcflag==0) { // uniform beam
				// set launch point and width of beam
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				r		= radius*sqrt(rnd); // radius of beam at launch point
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				phi		= rnd*2.0*PI;
				x		= xs + r*cos(phi);
				y		= ys + r*sin(phi);
				z		= zs;
				// set trajectory toward focus
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				r		= waist*sqrt(rnd); // radius of beam at focus
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				phi		= rnd*2.0*PI;
				xfocus	= r*cos(phi);
				yfocus	= r*sin(phi);
				temp	= sqrt((x - xfocus)*(x - xfocus) + (y - yfocus)*(y - yfocus) + zfocus*zfocus);
				ux		= -(x - xfocus)/temp;
				uy		= -(y - yfocus)/temp;
				uz		= sqrt(1 - ux*ux - uy*uy);
			}
			else if (mcflag==2) { // isotropic pt source
				costheta = 1.0 - 2.0*RandomGen(1,0,NULL);
				sintheta = sqrt(1.0 - costheta*costheta);
				psi = 2.0*PI*RandomGen(1,0,NULL);
				cospsi = cos(psi);
				if (psi < PI)
					sinpsi = sqrt(1.0 - cospsi*cospsi); 
				else
					sinpsi = -sqrt(1.0 - cospsi*cospsi);
				x = xs;
				y = ys;
				z = zs;
				ux = sintheta*cospsi;
				uy = sintheta*sinpsi;
				uz = costheta;
			}
			else if (mcflag==3) { // rectangular source collimated
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				x = radius*(rnd*2-1); // use radius to specify x-halfwidth of rectangle
				while ((rnd = RandomGen(1,0,NULL)) <= 0.0); // avoids rnd = 0
				y = radius*(rnd*2-1); // use radius to specify y-halfwidth of rectangle
				z = zs;
				ux = 0.0;
				uy = 0.0;
				uz = 1.0; // collimated beam
			}
		} // end  use mcflag
		/****************************/
		
		/* Get tissue voxel properties of launchpoint.
		 * If photon beyond outer edge of defined voxels, 
		 * the tissue equals properties of outermost voxels.
		 * Therefore, set outermost voxels to infinite background value.
		 */
		ix = (int)(Nx/2 + x/dx);
		iy = (int)(Ny/2 + y/dy);
		iz = (int)(z/dz);        
		if (ix>=Nx) ix=Nx-1;
		if (iy>=Ny) iy=Ny-1;
		if (iz>=Nz) iz=Nz-1;
		if (ix<0)   ix=0;
		if (iy<0)   iy=0;
		if (iz<0)   iz=0;		
		/* Get the tissue type of located voxel */
		i		= (long)(iz*Ny*Nx + ix*Ny + iy);
		type	= vX[i];
		mua 	= muav[type];
		mus 	= musv[type];
		g 		= gv[type];
		mua_mult = mua_multv[type];
		mus_mult = mus_multv[type];
		
		ftype = vF[i];
		qeff  = qeffv[ftype];
		tau   = tauv[ftype]*NANOSEC;
		muaxf = muaxfv[ftype];
		muaf_mult = muaf_multv[ftype];
		
		
		if (pfuncflag==1){
			// Determine alpha for the mhg phase function from the value of gamma
			num_alpha = gammav[type] - 3.0/5.0;
			denom_alpha = g*gammav[type]-g*g + 2.0/5.0;
			alpha = num_alpha/denom_alpha;
		}
		
        bflag = 1; // initialize as 1 = inside volume, but later check as photon propagates.
        
        
		/* HOP_DROP_SPIN_CHECK
		 Propagate one photon until it dies as determined by ROULETTE.
		 *******/
		do {
			
			/**** HOP
			 Take step to new position
			 s = dimensionless stepsize
			 x, uy, uz are cosines of current photon trajectory
			 *****/
			while ((rnd = RandomNum) <= 0.0);   /* yields 0 < rnd <= 1 */
			sleft	= -log(rnd);				/* dimensionless step */
			CNT += 1;
			
			do{  // while sleft>0
				update_flag =1;	   
				if (emflag==1) mus*=mus_mult;
				s     = sleft/mus;				/* Step size [cm].*/
				tempx = x + s*ux;				/* Update positions. [cm] */
				tempy = y + s*uy;	
				tempz = z + s*uz;
				
				sv = SameVoxel(x,y,z, tempx, tempy, tempz, dx,dy,dz);
				if (sv) /* photon in same voxel */
				{  
					x=tempx;					/* Update positions. */
					y=tempy;
					z=tempz;
					
					/** EXCITATION or EMISSION **/
					if ((qeff>0)&&(emflag==0)){ /* If current voxel is fluorescing */
						P = 1- exp(-muaxf*s);
						rnd = RandomNum;
						if (rnd<P){
								if (bflag) {
										// Save excitation fluence upto current point
								//NP: Since the probability of fluorescence already takes this into account, we need not include it here again. mua = mua+muaxf; // updated on 03052020
            					absorb = exp(-mua*s); // added on 01052020
	                    		temp_real = W*(cos(Wphase)-cos(-Wphase-omega_by_c0*s*n_in[xface[i]])*absorb);
	            				temp_imag = W*(-sin(Wphase)+sin(Wphase+omega_by_c0*s*n_in[xface[i]])*absorb);
								FX_real[i] += temp_real;	// only save data if blag==1, i.e., photon inside simulation cube
								FX_imag[i] += temp_imag;
							
								// added on 01052020 ----
								W *= absorb;	/* decrement WEIGHT by amount absorbed */		
								Wphase+=omega_by_c0*s*n_in[xface[i]];	
								//---------------- //	
							}
							
							emflag=1;
					//		fo_index=i;
								// choose random point in the current pixel
							x =  (ix-Nx/2)*dx + dx*(2*RandomGen(1,0,NULL)-1)/2;
							y =  (iy-Ny/2)*dy + dy*(2*RandomGen(1,0,NULL)-1)/2;
							z =  (iz)*dz + dz*(2*RandomGen(1,0,NULL)-1)/2;
							// choose random direction as isotropic point source
							costheta = 1.0 - 2.0*RandomGen(1,0,NULL);
							sintheta = sqrt(1.0 - costheta*costheta);
							psi = 2*PI*RandomNum; // choose a random direction over pi
							cospsi = cos(psi);
							if (psi < PI)
								sinpsi = sqrt(1.0 - cospsi*cospsi); 
							else
								sinpsi = -sqrt(1.0 - cospsi*cospsi);
							ux = sintheta*cospsi;
							uy = sintheta*sinpsi;
							uz = costheta;
							// Update weight
							rnd = -log(RandomNum);
							temp = omega*tau;
							W *= qeff/sqrt(1+temp*temp);
							Wphase -= atan2(temp,1);
							//Wphase = atan2(cos(Wphase)*temp + sin(Wphase),cos(Wphase)-temp*sin(Wphase));
							
							 //s=0; 
						 	sleft = -log(RandomNum);//sleft=0;
							 sleft -= s*(mus*mus_mult);
							 update_flag = 0;
						}
					}
					
								
					/**** DROP
					 Drop photon weight (W) into local bin.
					 *****/
					 
					 if (update_flag){
						 
						 if (emflag==1)
						 {
						 	mua =mua_mult*mua + muaxf*muaf_mult;
						 	absorb = exp(-mua*s);
						 
						 	if (bflag) {
						 			temp_real = W*(cos(Wphase)-cos(-Wphase-omega_by_c0*s*n_in[xface[i]])*absorb);
	            					temp_imag = W*(-sin(Wphase)+sin(Wphase+omega_by_c0*s*n_in[xface[i]])*absorb);
	            					FM_real[i] +=  temp_real;
	            					FM_imag[i] +=  temp_imag;
								 
							}
								W*=absorb;
								Wphase-=omega_by_c0*s*n_in[xface[i]];
						 	
						 }
						 else {
						 	//mua = mua+muaxf;
						 	absorb =  exp(-mua*s);	/* photon weight absorbed at this step */
	
	            
						// If photon within volume of heterogeneity, deposit energy in F[]. 
						// Normalize F[] later, when save output. 
	                    	if (bflag) {
	                    				temp_real = W*(cos(Wphase)-cos(-Wphase-omega_by_c0*s*n_in[xface[i]])*absorb);
	            				temp_imag = W*(-sin(Wphase)+sin(Wphase+omega_by_c0*s*n_in[xface[i]])*absorb);
						
									FX_real[i] += temp_real;
	            					FX_imag[i] += temp_imag;
						
							}
	                		W *= absorb;	/* decrement WEIGHT by amount absorbed */		
							Wphase+=omega_by_c0*s*n_in[xface[i]];
	               		}
	               		sleft =0;
					}
					
					/* Update sleft */
				//	sleft-=s*mus;
				//	if (sleft<=ls) sleft = 0;		/* dimensionless step remaining */
				}
				else /* photon has crossed voxel boundary */
				{
					/* step to voxel face + "littlest step" so just inside new voxel. */
					s = ls + FindVoxelFace2(x,y,z, tempx,tempy,tempz, dx,dy,dz, ux,uy,uz);
					
					if ((qeff>0)&&(emflag==0)){ /* If current voxel is fluorescing */
						P = 1- exp(-muaxf*s);
						rnd = RandomNum;
						if (rnd<P){
								if (bflag) {
										// Save excitation fluence upto current point
								//NP: Since the probability of fluorescence already takes this into account, we need not include it here again. mua = mua+muaxf; // updated on 03052020
            					absorb = exp(-mua*s); // added on 01052020
	                    		temp_real = W*(cos(Wphase)-cos(-Wphase-omega_by_c0*s*n_in[xface[i]])*absorb);
	            				temp_imag = W*(-sin(Wphase)+sin(Wphase+omega_by_c0*s*n_in[xface[i]])*absorb);
								FX_real[i] += temp_real;	// only save data if blag==1, i.e., photon inside simulation cube
								FX_imag[i] += temp_imag;
							
								// added on 01052020 ----
								W *= absorb;	/* decrement WEIGHT by amount absorbed */		
								Wphase+=omega_by_c0*s*n_in[xface[i]];	
								//---------------- //
								
							}
							
							emflag=1;
					//		fo_index=i;
								// choose random point in the current pixel
							x =  (ix-Nx/2)*dx + dx*(2*RandomGen(1,0,NULL)-1)/2;
							y =  (iy-Ny/2)*dy + dy*(2*RandomGen(1,0,NULL)-1)/2;
							z =  (iz)*dz + dz*(2*RandomGen(1,0,NULL)-1)/2;
							// choose random direction as isotropic point source
							costheta = 1.0 - 2.0*RandomGen(1,0,NULL);
							sintheta = sqrt(1.0 - costheta*costheta);
							psi = 2*PI*RandomNum; // choose a random direction over pi
							cospsi = cos(psi);
							if (psi < PI)
								sinpsi = sqrt(1.0 - cospsi*cospsi); 
							else
								sinpsi = -sqrt(1.0 - cospsi*cospsi);
							ux = sintheta*cospsi;
							uy = sintheta*sinpsi;
							uz = costheta;
							// Update weight
							rnd = -log(RandomNum);
							temp = omega*tau;
							W *= qeff/sqrt(1+temp*temp);
							Wphase   -= atan2(temp,1);
							//Wphase = atan2(cos(Wphase)*temp + sin(Wphase),cos(Wphase)-temp*sin(Wphase));	
						
							 //s=0; 
							 //sleft=0;
							 
							 sleft = -log(RandomNum);//sleft=0;
							 sleft -= s*(mus*mus_mult);
							 update_flag = 0;
						}
				}
					
					/**** DROP
					 Drop photon weight (W) into local bin.
					 *****/
					 if (update_flag==1){
					 
					  if (emflag==1)
					 {
					 	mua = mua_mult*mua + muaxf*muaf_mult;
					 	absorb = exp(-mua*s);
					 	if (bflag) {
					 			temp_real = W*(cos(Wphase)-cos(-Wphase-omega_by_c0*s*n_in[xface[i]])*absorb);
            					temp_imag = W*(-sin(Wphase)+sin(Wphase+omega_by_c0*s*n_in[xface[i]])*absorb);
            				
            					FM_real[i] += temp_real;
            					FM_imag[i] += temp_imag;
						
						}
						W*=absorb;
						Wphase-=omega_by_c0*s*n_in[xface[i]];
					 }
					 else { 
					 
					 	//	mua = mua+muaxf;
            			absorb = exp(-mua*s);              
					// If photon within volume of heterogeneity, deposit energy in F[]. 
					// Normalize F[] later, when save output. 
				
                    	if (bflag) {
                    				temp_real = W*(cos(Wphase)-cos(-Wphase-omega_by_c0*s*n_in[xface[i]])*absorb);
            				temp_imag = W*(-sin(Wphase)+sin(Wphase+omega_by_c0*s*n_in[xface[i]])*absorb);
					
								FX_real[i] += temp_real;
            					FX_imag[i] += temp_imag;
					
						}
							W *= absorb;	/* decrement WEIGHT by amount absorbed */		
							Wphase+=omega_by_c0*s*n_in[xface[i]];		
					}
                
					
					/* Update sleft */
					sleft -= s*mus;  /* dimensionless step remaining */
					if (sleft<=ls) sleft = 0;
					
					/* Update positions. */
					x += s*ux;
					y += s*uy;
					z += s*uz;
					}
					// pointers to voxel containing optical properties
                    ix = (int)(Nx/2 + x/dx);
                    iy = (int)(Ny/2 + y/dy);
                    iz = (int)(z/dz);
                    
                    bflag = 1;  // Boundary flag. Initialize as 1 = inside volume, then check.
                    if (boundaryflag==0) { // Infinite medium.
								// Check if photon has wandered outside volume.
                                // If so, set tissue type to boundary value, but let photon wander.
                                // Set blag to zero, so DROP does not deposit energy.
						if (iz>=Nz) {iz=Nz-1; bflag = 0; photon_status = DEAD;}
						if (ix>=Nx) {ix=Nx-1; bflag = 0;photon_status = DEAD;}
						if (iy>=Ny) {iy=Ny-1; bflag = 0;photon_status = DEAD;}
						if (iz<0)   {iz=0;    bflag = 0;photon_status = DEAD;}
						if (ix<0)   {ix=0;    bflag = 0;photon_status = DEAD;}
						if (iy<0)   {iy=0;    bflag = 0;photon_status = DEAD;}
					}
					else if (boundaryflag==1) { // Escape at boundaries
						rx_or_tx=0;
						
						if (iz<0)   {iz=0;    ri = (long) (ix*Ny + iy); un_in = -uz;	prev_det = 0;
							// These will be used if reflection occurs
						uz = -uz; z=-z;  rx_or_tx=1; }//iz = (int)(z/dz);
						
						if (iy<0)   {iy=0;    ri = (long) (iz*Nx + ix); un_in = -uy; prev_det =Nx*Ny;
							// These will be used if reflection occurs
						uy = -uy; y=-Ny*dy-y; rx_or_tx=1;}//iy = (int)(Ny/2 + y/dy)
						
						if (ix<0)   {ix=0;    ri = (long) (iz*Ny + iy); un_in = -ux; prev_det =Nx*Ny + Nx*Nz;
							// These will be used if reflection occurs
						ux = -ux; x=-Nx*dx-x; rx_or_tx=1;}//ix = (int)(Nx/2 + x/dx)
										
						if (iy>=Ny) {iy=Ny-1; ri = (long) (iz*Nx + ix);	un_in = uy;	prev_det = Nx*Ny + Ny*Nz + Nx*Nz; 	
						// These will be used if reflection occurs
						uy = -uy; y=Ny*dy-y;  rx_or_tx=1; } //iy = (int)(Ny/2 + y/dy);
						
						if (ix>=Nx) {ix=Nx-1; ri = (long) (iz*Ny + iy); un_in = ux;	prev_det = Nx*Ny + Ny*Nz + 2*Nx*Nz; 
						// These will be used if reflection occurs
						ux = -ux; x=Nx*dx-x;  rx_or_tx=1; } //ix = (int)(Nx/2 + x/dx);
						
						
						if (iz>=Nz) {iz=Nz-1; ri = (long) (ix*Ny + iy);	un_in = uz;	prev_det = Nx*Ny + 2*Ny*Nz + 2*Nx*Nz; 		
						// These will be used if reflection occurs
						uz = -uz; z=2*Nz*dz-z; rx_or_tx=1; }// iz = (int)(z/dz);
						
						
						if (rx_or_tx==1){
						
							rnd = RandomGen(1,0,NULL);														
							xtype = xface[i];
				
							Ref = RFresnel(n_in[xtype],n_out[xtype],un_in,&un_out);
							

							if (rnd>Ref)
							{
								// transmit current photon
								if (outflag==1){
									
										if (emflag==1){
												
										 RM_real[prev_det + ri]+= W*cos(Wphase)*un_out;
										 RM_imag[prev_det + ri]-= W*sin(Wphase)*un_out;
										// Forigin[fo_index]+=1;
										}
										else{
											RX_real[prev_det + ri]+= W*cos(Wphase)*un_out;
										 	RX_imag[prev_det + ri]-= W*sin(Wphase)*un_out;
										}		 
								}
								else {
										if (emflag==1){
												
										 RM_real[prev_det + ri]+= W*cos(Wphase);
										 RM_imag[prev_det + ri]-= W*sin(Wphase);
									//	 Forigin[fo_index]+=1;
										}
										else{
											RX_real[prev_det + ri]+= W*cos(Wphase);
										 RX_imag[prev_det + ri]-= W*sin(Wphase);
							 
										}
											
								}
								photon_status = DEAD; sleft=0;
								  bflag =0;
								//  SortScatEvents(num_scat, scatevents, x,y);
								  
							}
							else
							{
								sleft-=s*mus;
								bflag=1;
							}
					}
						
					}
					else if (boundaryflag==2) { // Escape at top surface, no x,y bottom z boundaries
						if (iz>=Nz) {iz=Nz-1; bflag = 0;photon_status = DEAD; sleft=0;} // setting photon to dead if it walks out of computational domain may be better idea in terms of comp. time
						if (ix>=Nx) {ix=Nx-1; bflag = 0;photon_status = DEAD; sleft=0;} 
						if (iy>=Ny) {iy=Ny-1; bflag = 0;photon_status = DEAD; sleft=0;} 
						if (ix<0)   {ix=0;    bflag = 0;photon_status = DEAD; sleft=0;}
						if (iy<0)   {iy=0;    bflag = 0;photon_status = DEAD; sleft=0;}
						// Only exiting on top surface
						if (iz<0)   {iz=0; 	
						ri = (long) (ix*Ny +iy);
						xtype = xface[ri];
						rnd = RandomGen(1,0,NULL);
						Ref = RFresnel(n_in[xtype],n_out[xtype],-uz,&uz1);
					
					
						
						//printf("look \t");
						if(rnd>=Ref)
						{
							if (outflag==1){
									if (emflag==1){
											
									 RM_real[ri]+= W*cos(Wphase)*uz1;
									 RM_imag[ri]-= W*sin(Wphase)*uz1;
								//	 Forigin[fo_index]+=1;
									}
									else{
										RX_real[ri]+= W*cos(Wphase)*uz1;
									 RX_imag[ri]-= W*sin(Wphase)*uz1;
									}
							}
								
							else {
							if (emflag==1){
											
									 RM_real[ri]+= W*cos(Wphase);
									 RM_imag[ri]-= W*sin(Wphase);
								//	 Forigin[fo_index]+=1;
									}
									else{
										RX_real[ri]+= W*cos(Wphase);
									 RX_imag[ri]-= W*sin(Wphase);
						 
									}
							}
							photon_status = DEAD; sleft=0;
							  bflag =0;
							//  SortScatEvents(num_scat, scatevents, x,y);
							  emflag =0; 			// reset emission flag after photon is lost to reflection
						//	printf("quack \n");
						}
						else
						{
							uz = -uz; z=-z;
							sleft-=s*mus;
							iz = (int)(z/dz);
							bflag=1;
						//	printf("duck \n");
						}
						}
					}
					
                    // update pointer to tissue type
					i    = (long)(iz*Ny*Nx + ix*Ny + iy);
					type = vX[i];
                    mua  = muav[type];
                    mus  = musv[type];
                    g    = gv[type];
                    gamma = gammav[type];
                    mua_mult = mua_multv[type];
                    mus_mult = mus_multv[type];
					
		
					ftype = vF[i];
					qeff  = qeffv[ftype];
					tau   = tauv[ftype]*NANOSEC;
					muaxf = muaxfv[ftype];
					muaf_mult = muaf_multv[ftype];
					
                    
                    if (pfuncflag==1){
						// Determine alpha for the mhg phase function from the value of gamma
						num_alpha = gammav[type] - 3.0/5.0;
						denom_alpha = g*gammav[type]-g*g + 2.0/5.0;
						alpha = num_alpha/denom_alpha;
					}
                    
				} //(sv) /* same voxel */
                
			} while(sleft>0); //do...while
			
			/**** SPIN 
			 Scatter photon into new trajectory defined by theta and psi.
			 Theta is specified by cos(theta), which is determined 
			 based on the Henyey-Greenstein scattering function.
			 Convert theta and psi into cosines ux, uy, uz. 
			 *****/
			/* Sample for costheta */
		//	num_scat+=1;
			rnd = RandomNum;
			if (pfuncflag==0){
					if (g == 0.0)
					costheta = 2.0*rnd - 1.0;
				else {
					double temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
					costheta = (1.0 + g*g - temp*temp)/(2.0*g);
				}
				sintheta = sqrt(1.0 - costheta*costheta); /* sqrt() is faster than sin(). */
			}
			else if (pfuncflag==1){
				if (rnd<=alpha)
				{		
					rnd = RandomNum;
					if (g == 0.0)
						costheta = 2.0*rnd - 1.0;
					else {
						double temp = (1.0 - g*g)/(1.0 - g + 2*g*rnd);
						costheta = (1.0 + g*g - temp*temp)/(2.0*g);
					}
					sintheta = sqrt(1.0 - costheta*costheta); /* sqrt() is faster than sin(). */
				}
				else
				{
					
					int loop=1;
					while(loop==1){
						rnd = RandomNum;
						costheta =2*RandomNum-1;
						if (rnd>0.5*(costheta*costheta)) loop=0;	 //(0.5 is normalisation constant. cross check use of this
						// Ref. importance sampling by Frisvald
					}
					sintheta = sqrt(1.0 - costheta*costheta); /* sqrt() is faster than sin(). */
				}
			}
			else{
				printf("Mie scattering not yet ready... Press Ctrl + C to abort.... \n");
				getchar();
			}	
			
			
			/* Sample psi. */
			psi = 2.0*PI*RandomNum;
			cospsi = cos(psi);
			if (psi < PI)
				sinpsi = sqrt(1.0 - cospsi*cospsi);     /* sqrt() is faster than sin(). */
			else
				sinpsi = -sqrt(1.0 - cospsi*cospsi);
			
			/* New trajectory. */
			if (1 - fabs(uz) <= ONE_MINUS_COSZERO) {      /* close to perpendicular. */
				uxx = sintheta * cospsi;
				uyy = sintheta * sinpsi;
				uzz = costheta * SIGN(uz);   /* SIGN() is faster than division. */
			} 
			else {					/* usually use this option */
				temp = sqrt(1.0 - uz * uz);
				uxx = sintheta * (ux * uz * cospsi - uy * sinpsi) / temp + ux * costheta;
				uyy = sintheta * (uy * uz * cospsi + ux * sinpsi) / temp + uy * costheta;
				uzz = -sintheta * cospsi * temp + uz * costheta;
			}
			
			/* Update trajectory */
			ux = uxx;
			uy = uyy;
			uz = uzz;
			
			/**** CHECK ROULETTE 
			 If photon weight below THRESHOLD, then terminate photon using Roulette technique.
			 Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
			 and 1-CHANCE probability of terminating.
			 *****/
			 if (emflag==1)
			 	{
				 if ((float)W < THRESHOLDM) {
					if (RandomNum <= CHANCE)
						W /= CHANCE;
					else photon_status = DEAD;
			 	}
			 }
			else
			{
				if ((float)W < THRESHOLDX) {
					if (RandomNum <= CHANCE)
						W /= CHANCE;
					else photon_status = DEAD;
				}				
			}
			
            
		} while (photon_status == ALIVE);  /* end STEP_CHECK_HOP_SPIN */
        /* if ALIVE, continue propagating */
		/* If photon DEAD, then launch new photon. */	
		if (i_photon==1000){
			temp_time = clock();
			time_min = (double)(temp_time-start_time)/CLOCKS_PER_SEC/60;
			printf("Time per photon is %5.6f min  \n", time_min/1000);
			printf("Estimated Time for %0.3e photons = %5.3f min\n",Nphotons,time_min*(Nphotons-1000)/1000);
			
		}
        
	} while (i_photon < Nphotons);  /* end RUN */
	
    if (printprogress==1)
    {
	printf("------------------------------------------------------\n");
	finish_time = clock();
	time_min = (double)(finish_time-start_time)/CLOCKS_PER_SEC/60;
	printf("Elapsed Time for %0.3e photons = %5.3f min\n",Nphotons,time_min);
	printf("%0.2e photons per minute\n", Nphotons/time_min);}
	
}



/* SUBROUTINES */

/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double RandomGen(char Type, long Seed, long *Status){
    static long i1, i2, ma[56];   /* ma[0] is not used. */
    long        mj, mk;
    short       i, ii;
    
    if (Type == 0) {              /* set seed. */
        mj = MSEED - (Seed < 0 ? -Seed : Seed);
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        for (i = 1; i <= 54; i++) {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }
        for (ii = 1; ii <= 4; ii++)
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }
        i1 = 0;
        i2 = 31;
    } else if (Type == 1) {       /* get a number. */
        if (++i1 == 56)
            i1 = 1;
        if (++i2 == 56)
            i2 = 1;
        mj = ma[i1] - ma[i2];
        if (mj < MZ)
            mj += MBIG;
        ma[i1] = mj;
        return (mj * FAC);
    } else if (Type == 2) {       /* get status. */
        for (i = 0; i < 55; i++)
            Status[i] = ma[i + 1];
        Status[55] = i1;
        Status[56] = i2;
    } else if (Type == 3) {       /* restore status. */
        for (i = 0; i < 55; i++)
            ma[i + 1] = Status[i];
        i1 = Status[55];
        i2 = Status[56];
    } else
        puts("Wrong parameter to RandomGen().");
    return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/***********************************************************
 *  Determine if the two position are located in the same voxel
 *	Returns 1 if same voxel, 0 if not same voxel.
 ****/				
Boolean SameVoxel(double x1,double y1,double z1, double x2, double y2, double z2, double dx,double dy,double dz)
{
    double xmin=min2((floor)(x1/dx),(floor)(x2/dx))*dx;
    double ymin=min2((floor)(y1/dy),(floor)(y2/dy))*dy;
    double zmin=min2((floor)(z1/dz),(floor)(z2/dz))*dz;
    double xmax = xmin+dx;
    double ymax = ymin+dy;
    double zmax = zmin+dz;
    Boolean sv=0;
    
    sv=(x1<=xmax && x2<=xmax && y1<=ymax && y2<=ymax && z1<zmax && z2<=zmax);
    return (sv);
}

/***********************************************************
 * max2
 ****/
double max2(double a, double b) {
    double m;
    if (a > b)
        m = a;
    else
        m = b;
    return m;
}

/***********************************************************
 * min2
 ****/
double min2(double a, double b) {
    double m;
    if (a >= b)
        m = b;
    else
        m = a;
    return m;
}
/***********************************************************
 * min3
 ****/
double min3(double a, double b, double c) {
    double m;
    if (a <=  min2(b, c))
        m = a;
    else if (b <= min2(a, c))
        m = b;
    else
        m = c;
    return m;
}

/********************
 * my version of FindVoxelFace for no scattering.
 * s = ls + FindVoxelFace2(x,y,z, tempx, tempy, tempz, dx, dy, dz, ux, uy, uz);
 ****/
double FindVoxelFace2(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz)
{	
    int ix1 = floor(x1/dx);
    int iy1 = floor(y1/dy);
    int iz1 = floor(z1/dz);
    
    int ix2,iy2,iz2;
    if (ux>=0)
        ix2=ix1+1;
    else
        ix2 = ix1;
    
    if (uy>=0)
        iy2=iy1+1;
    else
        iy2 = iy1;
    
    if (uz>=0)
        iz2=iz1+1;
    else
        iz2 = iz1;
    
    double xs = fabs( (ix2*dx - x1)/ux);
    double ys = fabs( (iy2*dy - y1)/uy);
    double zs = fabs( (iz2*dz - z1)/uz);
    
    double s = min3(xs,ys,zs);
    
    return (s);
}


/***********************************************************
 *	FRESNEL REFLECTANCE
 * Computes reflectance as photon passes from medium 1 to 
 * medium 2 with refractive indices n1,n2. Incident
 * angle a1 is specified by cosine value ca1 = cos(a1).
 * Program returns value of transmitted angle a1 as
 * value in *ca2_Ptr = cos(a2).
 ****/
double RFresnel(float n1,		/* incident refractive index.*/
                float n2,		/* transmit refractive index.*/
                double ca1,		/* cosine of the incident */
                /* angle a1, 0<a1<90 degrees. */
                double *ca2_Ptr) 	/* pointer to the cosine */
/* of the transmission */
/* angle a2, a2>0. */
{
    double r;
    
    if(n1==n2) { /** matched boundary. **/
        *ca2_Ptr = ca1;
        r = 0.0;
	}
    else if(ca1>(1.0 - 1.0e-12)) { /** normal incidence. **/
        *ca2_Ptr = ca1;
        r = (n2-n1)/(n2+n1);
        r *= r;
	}
    else if(ca1< 1.0e-6)  {	/** very slanted. **/
        *ca2_Ptr = 0.0;
        r = 1.0;
	}
    else  {			  		/** general. **/
        double sa1, sa2; /* sine of incident and transmission angles. */
        double ca2;      /* cosine of transmission angle. */
        sa1 = sqrt(1-ca1*ca1);
        sa2 = n1*sa1/n2;
        if(sa2>=1.0) {	
            /* double check for total internal reflection. */
            *ca2_Ptr = 0.0;
            r = 1.0;
		}
        else {
            double cap, cam;	/* cosines of sum ap or diff am of the two */
            /* angles: ap = a1 + a2, am = a1 - a2. */
            double sap, sam;	/* sines. */
            *ca2_Ptr = ca2 = sqrt(1-sa2*sa2);
            cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
            cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
            sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
            sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
            r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
            /* rearranged for speed. */
		}
	}
    return(r);
} /******** END SUBROUTINE **********/



/***********************************************************
 * the boundary is the face of some voxel
 * find the the photon's hitting position on the nearest face of the voxel and update the step size.
 ****/
double FindVoxelFace(double x1,double y1,double z1, double x2, double y2, double z2,double dx,double dy,double dz, double ux, double uy, double uz)
{
    double x_1 = x1/dx;
    double y_1 = y1/dy;
    double z_1 = z1/dz;
    double x_2 = x2/dx;
    double y_2 = y2/dy;
    double z_2 = z2/dz;
    double fx_1 = floor(x_1) ;
    double fy_1 = floor(y_1) ;
    double fz_1 = floor(z_1) ;
    double fx_2 = floor(x_2) ;
    double fy_2 = floor(y_2) ;
    double fz_2 = floor(z_2) ;
    double x=0, y=0, z=0, x0=0, y0=0, z0=0, s=0;
    
    if ((fx_1 != fx_2) && (fy_1 != fy_2) && (fz_1 != fz_2) ) { //#10
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        
        x = (max2(fx_1,fx_2)-x_1)/ux;
        y = (max2(fy_1,fy_2)-y_1)/uy;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (x == min3(x,y,z)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else if (y == min3(x,y,z)) {
            y0 = max2(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = max2(fz_1,fz_2);
            y0 = (z0-z_1)/uz*uy+y_1;
            x0 = (z0-z_1)/uz*ux+x_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fy_1 != fy_2) ) { //#2
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        x = (max2(fx_1,fx_2)-x_1)/ux;
        y = (max2(fy_1,fy_2)-y_1)/uy;
        if (x == min2(x,y)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            y0 = max2(fy_1, fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
    }
    else if ( (fy_1 != fy_2) &&(fz_1 != fz_2) ) { //#3
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        y = (max2(fy_1,fy_2)-y_1)/uy;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (y == min2(y,z)) {
            y0 = max2(fy_1,fy_2);
            x0 = (y0-y_1)/uy*ux+x_1;
            z0 = (y0-y_1)/uy*uz+z_1;
        }
        else {
            z0 = max2(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if ( (fx_1 != fx_2) && (fz_1 != fz_2) ) { //#4
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        x = (max2(fx_1,fx_2)-x_1)/ux;
        z = (max2(fz_1,fz_2)-z_1)/uz;
        if (x == min2(x,z)) {
            x0 = max2(fx_1,fx_2);
            y0 = (x0-x_1)/ux*uy+y_1;
            z0 = (x0-x_1)/ux*uz+z_1;
        }
        else {
            z0 = max2(fz_1, fz_2);
            x0 = (z0-z_1)/uz*ux+x_1;
            y0 = (z0-z_1)/uz*uy+y_1;
        }
    }
    else if (fx_1 != fx_2) { //#5
        fx_2=fx_1+SIGN(fx_2-fx_1);//added
        x0 = max2(fx_1,fx_2);
        y0 = (x0-x_1)/ux*uy+y_1;
        z0 = (x0-x_1)/ux*uz+z_1;
    }
    else if (fy_1 != fy_2) { //#6
        fy_2=fy_1+SIGN(fy_2-fy_1);//added
        y0 = max2(fy_1, fy_2);
        x0 = (y0-y_1)/uy*ux+x_1;
        z0 = (y0-y_1)/uy*uz+z_1;
    }
    else { //#7 
        z0 = max2(fz_1, fz_2);
        fz_2=fz_1+SIGN(fz_2-fz_1);//added
        x0 = (z0-z_1)/uz*ux+x_1;
        y0 = (z0-z_1)/uz*uy+y_1;
    }
    //s = ( (x0-fx_1)*dx + (y0-fy_1)*dy + (z0-fz_1)*dz )/3;
    //s = sqrt( SQR((x0-x_1)*dx) + SQR((y0-y_1)*dy) + SQR((z0-z_1)*dz) );
    //s = sqrt(SQR(x0-x_1)*SQR(dx) + SQR(y0-y_1)*SQR(dy) + SQR(z0-z_1)*SQR(dz));
    s = sqrt( SQR((x0-x_1)*dx) + SQR((y0-y_1)*dy) + SQR((z0-z_1)*dz));
    return (s);
}
void SortScatEvents(int nofevents,int *outarray, double xcoord, double ycoord)
{
	float dist_frm_src;
	float dist_bins[]= {0,0.01,0.02,0.03,0.04,0.05,0.06,0.08,0.1};
	int Nbins =8;
	int i,j;
	
	dist_frm_src = sqrt(xcoord*xcoord + ycoord*ycoord);
	
	for (i=0;i<Nbins+1;i++)
	{
		if ((dist_frm_src<=dist_bins[i+1])&&(dist_frm_src>dist_bins[i]))
			{
				if (nofevents<=2) outarray[i*5+0]+=1;
				else if ((nofevents<=4)&&(nofevents>2)) outarray[i*5+1]+=1;
				else if ((nofevents<=6)&&(nofevents>4)) outarray[i*5+2]+=1;
				else if ((nofevents<=8)&&(nofevents>6)) outarray[i*5+3]+=1;
				else outarray[i*5+4]+=1;
			}	
	}
	
}

