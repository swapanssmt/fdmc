# Frequency Domain Monte Carlo code for light propagation with fluorescence (FreDMC)

The **FreDMC** program simulates light propagation through a voxelised medium. 

Present version features

1. Simulation of heterogenous media.
2. Variable refractive index on different sections of the boundary. 
3. Light sources: Boundary sources; Collimated, Gaussian, isotropic; broad beam and point sources; modulation
4. Measurements: Partial Current and Exiting fluence; reflectance/ transmission mode; 
5. Supports fluorescence. 

A beta version **FreDMC_ang** is also included. In addition to the above features, it supports specification detector fiber properties such as Numerical aperture, fiber orientation, refractive index and fiber diameter, in reflectance mode only. 

## Package inclusions

The package contains the following files 

**I. MATLAB files**

1. `maketissue.m` 
2. `maketissue_ang.m`
3. `maketissueList.m`
4. `makeFluorList.m`
5. `makeXfaceList.m`
6. `lookmcxyz_m_freq.m`
7. `lookmcxyz_ang_freq.m`

These are used for generating the tissue geometry and for post-processing and viewing the **FreDMC** & **FreDMC_ang** output.

**II. C- files** 

1. `mcfl_main.c`
2. `mcfl_3df_lib.c`
3. `mcfl_3df_lib.h`
4. `main_ang.c` (beta)
5. `mcfl_3df_ang_lib.c` (beta)
6. `mcfl_3df_ang_lib.h`

These are the main c files that of the **FreDMC**  and **FredMC_ang** program. *Always make copies before modifying.* 

**III. Sample files **

1. `FredMC.exe`
2. `FredMC_ang.exe`
3. `sample_tissue.mci`
4. `sample_tissue_TX.bin`
5. `sample_tissue_TF.bin`
6. `sample_tissue_XFACE.bin`
7. `sample_ang.mci`
8. `sample_ang_TX.bin`
9. `sample_ang_TF.bin`
10. `sample_ang_XFACE.bin`

We include a pre-compiled `exe` file for the program. The remaining files provide a sample input for a cuboidal domain with a single fluorescing heterogeneity. 

## Installation

### 1. Compiling the code via command prompt

- Open command prompt and navigate to the folder which contains the `c` files. 

- To compile the file with the gcc-compiler type in the following commands and press enter.

  > gcc -o mymcfile mcfl_main.c mcfl_3df_lib.c -lm

- `mymcfile.exe` will be the name of the compiled file. 

  

### 2. Compiling using an IDE

- Open the program `mcfl_main.c` in your IDE. 

- Navigate to the **Compiler options** from the **Taskbar** or **Tools menu**. 

- Choose the gcc compiler. In the General setting, check the tickbox for "Add the following commands when calling the compiler"

- In the space below the checkbox, type in 

  > -g3 mcfl_3df_lib.c -lm

- Save the options and the compile the file `mcfl_main.c`. 

- This will generate `mcfl_main.exe`. 

*Note: these instructions are based on the DevC++ IDE. There might be slight variation on how to set additional calls for the compiler for different IDEs.

### 3. Setting up the simulation

The MATLAB interface is used for setting up the simulation through the script `maketissue.m`. 

The script takes in the user specified inputs such as the bin spacing, # of bins, # of photons, modulation frequency, source types etc. and rewrites them into a monte carlo input (MCI) file, which is generally named as `myname_H.mci`. 
The code also allows the user to specify the properties of each pixel through a tissue-type , fluorophore-type and interface-type pointer, which are then written to the files `myname_TX.bin`, `myname_TF.bin` and `myname_XFACE.bin` respectively.  These files are read by the FreDMC program.

The tissue-types, fluorophore-types and interface types, are stored in the scripts `maketissueList.m, makeFluorList.m` and `makeXfaceList.m`. Additional tissue types can be created by adding to the list in these files.   

By assigning different tissue-types to different voxels, heterogeneous tissue can be simulated. 

### 4. Running the program

To run the program, navigate to the folder that contains the compiled `exe` file, as well as the input files created by the MATLAB script. At the command prompt type

> mymcfile myname

where `mymcfile` is the name of the exe file, and `myname` is the base-name for the input files.

### 5. Post-processing and visualisation

A successful run of the FreDMC program, generates the following output files 

`myname_FXreal.bin`, `myname_FMreal.bin`, `myname_FXimag.bin`, `myname_FMimag.bin`, `myname_RXreal.bin`,`myname_RXimag.bin`, `myname_RMreal.bin`, `myname_RMimag.bin`. 

These files are read by the MATLAB script `loomcxyz_m_freq.m` , and plots of the total fluence and total reflectance at both excitation and emission are generated. 

## Usage

The simulation is set up through the file `maketissue.m`. It allows the user to specify the following parameters

1. Nphotons - \# of photons

2. binsize - Bin/ voxel size.  The default shape of the voxel is a cube. Scaled values of binsize can be used to specify cuboidal elements. 

3. Nbins - \# of bins. The computational domain will have a physical dimension of (Nbins X binsize) units^3. 

   **Source parameters **

4. mod_freq - Modulation frequency of the source in MHz. 

5. mcflag - specifies type of source. Choose between, (0)- Uniform beam , (1) - Gaussian,  (2) - Isotropic point.  

6. emflag - choose to save (1) or reject (0) emission data. 

7. (xs,ys,zs) - source position. Currently support only a single source in a given simulation. 

8. (xfocus,yfocus,zfocus) - focusing of the beam , set to inf in direction of propagation for simulating a uniform/ collimated beam.

9. radius - radius of the beam. 

10. waist - beam-waist. 

11. launchflag - choose launch trajectory. (0) - computed in the C-program. Always set to 0 , for mcflag = 2. 

    (1) - Launch from z = zmin/zmax, (2) - launch from x = xmin/xmax, (3) - launch from y = ymin/ymax. 

    When the launchflag is non-zero, the matlab script computes the direction cosines (ux0,uy0,uz0) to be used at launch. For the same, the user can specify:

    **theta_inc** -polar angle at incidence (outside the domain, wrt surface normal)

    **phi_inc** - azimuthal angle at incidence

    refr_source - refractive index at the source position. (this is used to compute the launch angle inside the domain). 

    **Boundaries**

12. boundaryflag - specify refractive index mismatch at boundaries. (0) - no mismatch/ no- boundaries,  (1) - mismatch at all boundaries, (2) - mismatch only at the z= 0 face. 

    In all cases the photon is marked as dead, once it exits the computational domain. In order to account for re-entry of the photon from outside the computational domain when simulating open boundaries (cases 0,2), recommend setting the computational domain to slightly larger than that desired, i.e. instead of simulating a 1cm^3 domain, set a domain of 1.2 cm^3. In most cases this will account for the small fraction of re-entry of photons from outside the desired computational domain. 

13.  Interface list - Interfaces with different refractive index combinations (refr_in & refr_out) can be specifies through the interface library `makeXfaceList`.  Corresponding to each voxel, an interface pointer is created and saved as `_XFACE.bin`. The default value for each voxel is 3, which assumes no mismatch at the current voxel. 

    **Optical Properties**

14. pfuncflag - choose the phase function. Currently supports  (0) - Henyey Greenstein phase function , and (1) - Modified Henyey Greenstein phase function. 

15.  Optical property list - The optical p

    

