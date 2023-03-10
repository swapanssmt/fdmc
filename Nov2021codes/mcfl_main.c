#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "mcfl_3df_lib.h"

#define Ntiss		50          /* Number of tissue types. */
#define NXFACE		10			/* Max Number of interface types in library */
#define STRLEN 		32          /* String length. */
#define Boolean     char
#define	LIGHTSPEED	2.997925E10 /* in vacuo speed of light [cm/s] */
#define	PI          3.1415926
#define	MHZ			1.0E6

int main(int argc, const char * argv[]) {
    
    if (argc==0) {
        printf("assuming you've compiled mcxyz.c as gomcxyz ...\n");
		printf("USAGE: gomcxyz name\n");
		printf("which will load the files name_H.mci and name_T.bin\n");
		printf("and run the Monte Carlo program.\n");
		printf("Yields  name_F.bin, which holds the fluence rate distribution.\n");
        return 0;
    }
    
    	
	double	Nphotons, Fphotons;       /* number of photons in simulation */
	int 	emflag;			/*set to 1 to do a fluorescence simulation */
	
	
	/* launch parameters */
	int		mcflag, launchflag, boundaryflag, outflag, pfuncflag;
	int 	printprogress=1;
	float	xfocus, yfocus, zfocus;
	float	ux0, uy0, uz0;
	float	radius;
	float	waist;
	float 	mod_freq; // modulation frequency
	
	/* dummy variables */

	long	i,j,NN, Nxy;         /* dummy indices */
	int 	ix, iy, iz;     /* Added. Used to track photons */
	double 	temp;           /* dummy variable */
    double 	mua;
	long 	nphotons_fluor;
	
	/* mcxyz bin variables */
	float	dx, dy, dz;     /* bin size [cm] */
	int		Nx, Ny, Nz, Nt, Nf; /* # of bins */
	int 	nxface;
	float	xs, ys, zs;		/* launch position */
	
    
    /* time */
	float	time_min;               
	time_t	now;
	double	start_time, finish_time, temp_time, time_1; /* for clock() */
	
	/* tissue parameters */
	char	tissuename[50][32];
	float 	muav[Ntiss];            // muav[0:Ntiss-1], absorption coefficient of ith tissue type
	float 	musv[Ntiss];            // scattering coeff. 
	float 	gv[Ntiss];              // anisotropy of scattering
	float	gammav[Ntiss];			// gamma parameter for modified Henyay Greenstein
	float	muaxfv[Ntiss];			// absorption coefficient of fluorophore
	float 	Qeffv[Ntiss];			// Quantum efficiency of fluorophore
	float	tauv[Ntiss];			// fluorescence lifetime of fluorophore
	float	n_inv[NXFACE];			// Refractive index inside
	float	n_outv[NXFACE];			// refractive index outside
	float	mua_multv[Ntiss];		// multiplication factor for absptn coeff at emission wavelength
	float	mus_multv[Ntiss];		// multiplication factor for mus at emission wavelength
	float	muaf_multv[Ntiss];		// multiplying factor to obtain muaf at emission wavelength
	
    
	/* Input/Output */
	char   	myname[STRLEN];		// Holds the user's choice of myname, used in input and output files. 
	char	filename[STRLEN];     // temporary filename for writing output.
    FILE*	fid=NULL;               // file ID pointer 
    char    buf[32];                // buffer for reading header.dat
    
    
	strcpy(myname, argv[1]);    // acquire name from argument of function call by user.
    //strcpy(myname, "sample");
    
	printf("name = %s\n",myname);
    
	/**** INPUT FILES *****/
    /* IMPORT myname_H.mci */
    strcpy(filename,myname);
    strcat(filename, "_H.mci");
	fid = fopen(filename,"r");
	fgets(buf, 32, fid);
		// run parameters
		sscanf(buf, "%lf", &Nphotons); // desired time duration of run [min]
		fgets(buf, 32, fid);
		sscanf(buf, "%d", &Nx);  // # of bins  
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &Ny);  // # of bins
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &Nz);  // # of bins   
	
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &dx);	 // size of bins [cm]
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &dy);	 // size of bins [cm] 
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &dz);	 // size of bins [cm] 
	
		// launch parameters
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &mcflag);  // mcflag, 0 = uniform, 1 = Gaussian, 2 = iso-pt
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &launchflag);  // launchflag, 0 = ignore, 1 = manually set
        fgets(buf, 32,fid);
        sscanf(buf, "%d", &boundaryflag);  // 0 = no boundaries, 1 = escape at all boundaries, 2 = escape at surface only
        fgets(buf,32,fid);
        sscanf(buf,"%d", &emflag); // 1 = with fluorescence, 0 = without fluorescence
        
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &xs);  // initial launch point
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &ys);  // initial launch point 
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &zs);  // initial launch point
	
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &xfocus);  // xfocus
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &yfocus);  // yfocus
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &zfocus);  // zfocus

		fgets(buf, 32,fid);
		sscanf(buf, "%f", &ux0);  // ux trajectory
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &uy0);  // uy trajectory
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &uz0);  // uz trajectory

		fgets(buf, 32,fid);
		sscanf(buf, "%f", &radius);  // radius
		fgets(buf, 32,fid);
		sscanf(buf, "%f", &waist);  // waist
		
		fgets(buf,32,fid);
    	sscanf(buf,"%f", &mod_freq); // modulation frequency
	
		// tissue optical properties
		fgets(buf,32,fid);
        sscanf(buf,"%d", &pfuncflag); 
		fgets(buf, 32,fid);
		sscanf(buf, "%d", &Nt);				// # of tissue types in tissue list
		for (i=1; i<=Nt; i++) {
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &muav[i]);	// absorption coeff [cm^-1]
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &musv[i]);	// scattering coeff [cm^-1]
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &gv[i]);		// anisotropy of scatter [dimensionless]
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &gammav[i]);		// gamma parameter for modified Henyey Greenstein phase function
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &mua_multv[i]);		// multiplication factor for abs coeff at emission wavelength
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &mus_multv[i]);		// multiplication factor for scat coeff at emission wavelength
		}    
		fgets(buf, 32, fid);
		sscanf(buf, "%d", &nxface);
		for (i=1; i<=nxface; i++){
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &n_inv[i]);
			fgets(buf, 32, fid);
			sscanf(buf, "%f", &n_outv[i]);
		}
	
		
		
		fgets(buf,32,fid);
		sscanf(buf,"%d",&Nf);
			for (i=1; i<=Nf; i++) {
				fgets(buf, 32, fid);
				sscanf(buf, "%f", &muaxfv[i]);	// absorption coeff [cm^-1]
				fgets(buf, 32, fid);
				sscanf(buf, "%f", &muaf_multv[i]);	// mult factor for absorption coeff [cm^-1]
				fgets(buf, 32, fid);
				sscanf(buf, "%f", &Qeffv[i]);	// absorption coeff [cm^-1]
				fgets(buf, 32, fid);
				sscanf(buf, "%f", &tauv[i]);	// absorption coeff [cm^-1]
			} 
		fgets(buf,32,fid);
		sscanf(buf,"%d",&outflag);
		
    fclose(fid);
    
    printf("Number of photons = %lf \n",Nphotons);
    printf("Nx = %d, dx = %0.4f [cm]\n",Nx,dx);
    printf("Ny = %d, dy = %0.4f [cm]\n",Ny,dy);
    printf("Nz = %d, dz = %0.4f [cm]\n",Nz,dz);

    printf("xs = %0.4f [cm]\n",xs);
    printf("ys = %0.4f [cm]\n",ys);
    printf("zs = %0.4f [cm]\n",zs);
    printf("mcflag = %d\n",mcflag);
    if (mcflag==0) printf("launching uniform flat-field beam\n");
    if (mcflag==1) printf("launching Gaissian beam\n");
    if (mcflag==2) printf("launching isotropic point source\n");
    if (mcflag==3) printf("launching square source\n");
    printf("xfocus = %0.4f [cm]\n",xfocus);
    printf("yfocus = %0.4f [cm]\n",yfocus);
    printf("zfocus = %0.2e [cm]\n",zfocus);
	if (launchflag==1) {
		printf("Launchflag ON, so launch the following:\n");
		printf("ux0 = %0.4f [cm]\n",ux0);
		printf("uy0 = %0.4f [cm]\n",uy0);
		printf("uz0 = %0.4f [cm]\n",uz0);
	}
	else {
		printf("Launchflag OFF, so program calculates launch angles.\n");
		printf("radius = %0.4f [cm]\n",radius);
		printf("waist  = %0.4f [cm]\n",waist);
	}
    if (boundaryflag==0)
		printf("boundaryflag = 0, so no boundaries.\n");
    else if (boundaryflag==1)
		printf("boundaryflag = 1, so escape at all boundaries.\n");    
	else if (boundaryflag==2)
		printf("boundaryflag = 2, so escape at surface only.\n");    
	else{
        printf("improper boundaryflag. quit.\n");
        return 0;
    }
    printf("Choice of phase function : \t");
    if (pfuncflag==0) printf("Henyey Greenstein phase function \n");
    else if (pfuncflag==1) printf("Modified Henyey Greenstein phase function \n");
    else if (pfuncflag ==2) printf("Mie scattering phase function \n");
    else {
    	printf("Invalid choice. Quit \n");
    	return 0;
	}
   /* if (outflag==0)
    	printf("output is diffuse reflectance. \n");
    else
    	printf("Output is partial current \n"); */
    	
    printf("# of tissues available, Nt = %d\n",Nt);
    for (i=1; i<=Nt; i++) {
        printf("muav[%ld] = %0.4f [cm^-1]\n",i,muav[i]);
        printf("musv[%ld] = %0.4f [cm^-1]\n",i,musv[i]);
        printf("  gv[%ld] = %0.4f [--]\n\n",i,gv[i]);
        printf("  gammav[%ld] = %0.4f [--]\n\n",i,gammav[i]);
        printf("mua_multv[%ld] = %0.4f [cm^-1]\n",i,mua_multv[i]);
        printf("mus_multv[%ld] = %0.4f [cm^-1]\n",i,mus_multv[i]);
    }
     for (i=1; i<=nxface; i++){
    printf("ref index [%d] in = %f \t",i, n_inv[i]);
    printf("ref index [%d] out = %f \n",i, n_outv[i]);
	}
	
    //printf("ref index in = %lf \n", n_in);
    //printf("ref index in = %lf \n", n_out);
    
   	for (i=1; i<=Nf; i++) {
        printf("muafv[%ld] = %0.4f [cm^-1]\n",i,muaxfv[i]);
        printf("Qeffv[%ld] = %0.4f \n",i,Qeffv[i]);
         printf("Tauv[%ld] = %0.4f [ns]\n",i,tauv[i]);
    }
	
	if (outflag==1) printf("Output is partial current \n");
	else printf("Output is reflectance. \n");


    // SAVE optical properties, for later use by MATLAB.
/*	strcpy(filename,myname);
	strcat(filename,"_props.m");
	fid = fopen(filename,"w");
	for (i=1; i<=Nt; i++) {
		fprintf(fid,"muav(%ld) = %0.4f;\n",i,muav[i]);
		fprintf(fid,"musv(%ld) = %0.4f;\n",i,musv[i]);
		fprintf(fid,"gv(%ld) = %0.4f;\n\n",i,gv[i]);
		fprintf(fid,"gammav(%ld) = %0.4f;\n\n",i,gammav[i]);
	}
	fclose(fid); */
    
    /* IMPORT BINARY TISSUE FILE */
	char 	*vX=NULL;
	char 	*xface = NULL;
	double 	*FX_real=NULL;
	double 	*FX_imag=NULL;
	double	*RX_real=NULL;
	double	*RX_imag=NULL;
	

	NN = Nx*Ny*Nz;
	if (boundaryflag==1) Nxy = 2*(Nx*Ny + Ny*Nz + Nx*Nz);
	else Nxy = Nx*Ny;
	
	vX  = ( char *)malloc(NN*sizeof(char));  /* tissue structure */
	xface  = ( char *)malloc(NN*sizeof(char));  /* tissue structure */
	
	FX_real  = (double *)malloc(NN*sizeof(double));	/* relative fluence rate [W/cm^2/W.delivered] */
	FX_imag  = (double *)malloc(NN*sizeof(double));	/* relative fluence rate [W/cm^2/W.delivered] */
	RX_real  = (double *)malloc(Nxy*sizeof(double));
	RX_imag  = (double *)malloc(Nxy*sizeof(double));
	
	// read binary file
    strcpy(filename,myname);
    strcat(filename, "_TX.bin");
    fid = fopen(filename, "rb");
    fread(vX, sizeof(char), NN, fid);
    fclose(fid);
    
    // read binary file
    strcpy(filename,myname);
    strcat(filename, "_XFACE.bin");
    fid = fopen(filename, "rb");
    fread(xface, sizeof(char), NN, fid);
    fclose(fid);
    
    
    // Show tissue on screen, along central z-axis, by listing tissue type #'s.
    iy = Ny/2;
    ix = Nx/2;
    printf("central axial profile of tissue types:\n");
    for (iz=0; iz<Nz; iz++) {
        i = (long)(iz*Ny*Nx + ix*Ny + iy);
        printf("%d",vX[i]);
    }
    printf("\n\n");
	
	// Fluorescence simulation 
			
	char	*vF=NULL;
	double 	*FM_real=NULL;
	double 	*FM_imag=NULL;
	double	*RM_real=NULL;
	double	*RM_imag=NULL;
	//double 	*Forigin=NULL;
	
	vF	= (char *)malloc(NN*sizeof(char)); /*location of fluorophore */
	FM_real  = (double *)malloc(NN*sizeof(double));	/* relative fluence rate [W/cm^2/W.delivered] */
	FM_imag  = (double *)malloc(NN*sizeof(double));	/* relative fluence rate [W/cm^2/W.delivered] */
	RM_real  = (double *)malloc(Nxy*sizeof(double));
	RM_imag  = (double *)malloc(Nxy*sizeof(double));	

	//Forigin = (double *)malloc(NN*sizeof(double));

	// read binary file
	strcpy(filename,myname);
    strcat(filename, "_TF.bin");
    fid = fopen(filename, "rb");
    fread(vF, sizeof(char), NN, fid);
    fclose(fid);	
    
    printf("central axial profile of fluorophore types:\n");
    for (iz=0; iz<Nz; iz++) {
        i = (long)(iz*Ny*Nx + ix*Ny + iy);
        printf("%d",vF[i]);
    }
    printf("\n\n");
    
    printf("------------- Begin FreDMC (w/ Fluorescence) -------------\n");
    printf("%s\n",myname);
    //Nphotons = 1000;
    	
    mcxyz(Nphotons, dx, dy, dz, Nx, Ny, Nz, xs, ys, zs, ux0, uy0, uz0, mod_freq, mcflag, launchflag, boundaryflag, pfuncflag, outflag, printprogress, n_inv, n_outv, xfocus, yfocus, zfocus,
	    	 radius, waist, muav, musv, gv, gammav,mua_multv, mus_multv, muaxfv, muaf_multv, Qeffv, tauv, FX_real, FX_imag, RX_real, RX_imag, FM_real, FM_imag, RM_real, RM_imag, vX, vF, xface);
	    	
		 
	/**** SAVE
     Convert data to relative fluence rate [cm^-2] and save.
     *****/
    printf("Back to main loop \n");
 //   getchar();
    // Normalize deposition (A) to yield fluence rate (F).
    temp = dx*dy*dz*Nphotons;
    double atemp,btemp,omega_by_c;
    double omega_by_c0 = 2*PI*mod_freq*MHZ/LIGHTSPEED;
    for (i=0; i<NN;i++){
        	atemp = FX_real[i]; btemp = FX_imag[i];
    	omega_by_c = omega_by_c0*n_inv[xface[i]];
    
		// Impt use only mua here, not mua + muaxf
      FX_real[i] = (double) (btemp*omega_by_c + atemp*muav[vX[i]])/(omega_by_c*omega_by_c + muav[vX[i]]*muav[vX[i]])/temp;
      FX_imag[i] = (double) -(atemp*omega_by_c - btemp*muav[vX[i]])/(omega_by_c*omega_by_c + muav[vX[i]]*muav[vX[i]])/temp;
      
    }
    
	temp = dx*dy*Nphotons;
    for (i=0; i<Nxy; i++)
    {
    	RX_real[i]/=temp; RX_imag[i]/=temp; 
	}
    
    // Save the binary file
    strcpy(filename,myname);
    strcat(filename,"_FXreal.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(FX_real, sizeof(double), NN, fid);
    fclose(fid);
    
    strcpy(filename,myname);
    strcat(filename,"_FXimag.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(FX_imag, sizeof(double), NN, fid);
    fclose(fid);
    
      /* save reflectance */

	strcpy(filename,myname);
    strcat(filename,"_RXreal.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 2D voxel output */
    fwrite(RX_real, sizeof(double), Nxy, fid);
    fclose(fid);
    
    strcpy(filename,myname);
    strcat(filename,"_RXimag.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 2D voxel output */
    fwrite(RX_imag, sizeof(double), Nxy, fid);
    fclose(fid);
   	
    
    if (emflag==1){
	
    temp = dx*dy*dz*Nphotons;
    for (i=0; i<NN; i++){
    	  atemp = FM_real[i]; btemp = FM_imag[i];
    	  	omega_by_c = omega_by_c0*n_inv[xface[i]];
    	  	
    	  	// Impt must use mua*mua_mult + muaxf*muaf_mult here. 
      mua = muav[vX[i]]*mua_multv[vX[i]] + muaxfv[vF[i]]*muaf_multv[vF[i]];
      FM_real[i] = (btemp*omega_by_c + atemp*mua)/(omega_by_c*omega_by_c + mua*mua)/temp;
      FM_imag[i] = -(atemp*omega_by_c - btemp*mua)/(omega_by_c*omega_by_c + mua*mua)/temp;
       
	}
	
	temp = dx*dy*Nphotons;
    for (i=0; i<Nxy; i++)
    {
    	RM_real[i]/=temp; RM_imag[i]/=temp;
	}
  
   // Emission data
    strcpy(filename,myname);
    strcat(filename,"_FMreal.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(FM_real, sizeof(double), NN, fid);
    fclose(fid);
    
    strcpy(filename,myname);
    strcat(filename,"_FMimag.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 3D voxel output */
    fwrite(FM_imag, sizeof(double), NN, fid);
    fclose(fid);
    
   
    /* save reflectance */

	strcpy(filename,myname);
    strcat(filename,"_RMreal.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 2D voxel output */
//	int Nyx = Ny*Nx;
    fwrite(RM_real, sizeof(double), Nxy, fid);
    fclose(fid);
    
    strcpy(filename,myname);
    strcat(filename,"_RMimag.bin");
    printf("saving %s\n",filename);
    fid = fopen(filename, "wb");   /* 2D voxel output */
    fwrite(RM_imag, sizeof(double), Nxy, fid);
    fclose(fid);
	
	}

    free(vF);
    free(FM_real);
    free(FM_imag);
    free(RM_real);
    free(RM_imag);
    //free(Forigin);
     
    
    free(xface);
    free(vX);
    free(FX_real);
    free(FX_imag);
    free(RX_real);
    free(RX_imag);
  
  	printf("%s is done.\n",myname);	
    printf("------------------------------------------------------\n");
    now = time(NULL);
    printf("%s\n", ctime(&now));
    return 0;   
}
