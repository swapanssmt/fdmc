% Modified on 21st August 2020
% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use,
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%
%Edited to include fluorescence properties.

clear
format compact
clc
home

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci
% 0 = don't save. Just check the program.

myname      = 'sample_ang';% name for files: myname_T.bin, myname_H.mci
time_min    = 10;      	% time duration of the simulation [min] <----- run time -----
Nphotons    = 1e6;
nm          = 532;   	% desired wavelength of simulation
Nbins       = 50;    	% # of bins in each dimension of cube
binsize     = 0.0125; 	% size of each bin, eg. [cm] or [mm]

mod_freq = 100; % modulation frequency in [MHZ]
pfuncflag =0; % Choice of phase function  0 - HG phase function. 1- Modified HG phase function

% Set Monte Carlo launch flags
mcflag      = 0;     	% launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt.
% 3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
launchflag  = 1;        % 0 = let mcxyz.c calculate launch trajectory
% 1 = manually set launch vector. % 2 = set through code. works for expt
% only.
boundaryflag = 2;       % 0 = no boundaries, 1 = escape at boundaries
% 2 = escape at surface only. No x, y, bottom z
% boundaries
emflag = 1; % 1 = w/ fluorescence, 0 = w/o fluorescence

% Sets position of source
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 0;  	% z of source

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus      = 0;        % set x,position of focus
yfocus      = 0;        % set y,position of focus
zfocus      = inf;    	% set z,position of focus (=inf for collimated beam)

% only used if mcflag == 0 or 1 or 3 (not 2=isotropic pt.)
radius      = 0.000000005; %0.00400;   % 1/e radius of beam at tissue surface
waist       = 0.000005; %0.00400;  	% 1/e radius of beam at focus

% only used if launchflag == 1 (manually set launch trajectory):
ux0 = 0; uy0 = 0; uz0 = 0;
if launchflag~=0
    % For oblique incidence
    theta_inc = 45; % angle of incidence in degrees
    phi_inc = 0;
%     uz=cos(theta_inc*pi/180);
%     ux = sqrt(1-uz^2)*cos(phi_inc*pi/180);
%     uy = sqrt(1-uz^2)*sin(phi_inc*pi/180);
    refr_source =1; refr_medium = 1.37;
    sin_theta_in =  (refr_source/refr_medium)*sin(theta_inc*pi/180); % angle inside the medium
%     sin_theta_in = sin(theta_inc*pi/180);
    
%     % ref. monte carlo report by wang for these expressions

switch (launchflag)
    case 1
    ux0 = sin_theta_in;
    uy0 = 0;
    uz0 = sqrt(1-sin_theta_in^2);
    case 2
    ux0 = sqrt(1-sin_theta_in^2);
    uy0 = 0;
    uz0 = sin_theta_in;
   
    case 3
        ux0 = 0;
        uy0 = sin_theta_in;
        uz0 = sqrt(1-sin_theta_in^2);
    case 4
        ux0 = 0;
        uy0 = sqrt(1-sin_theta_in^2);
        uz0 = sin_theta_in;
end
%     ux0         = (sin_theta_in/sqrt(1-uz^2))*(ux*uz*cos(phi_inc*pi/180)-uy*sin(phi_inc*pi/180)) - ux*sqrt(1-sin_theta_in^2);      % trajectory projected onto x axis
%     uy0         = (sin_theta_in/sqrt(1-uz^2))*(uy*uz*cos(phi_inc*pi/180)+ux*sin(phi_inc*pi/180)) - uy*sqrt(1-sin_theta_in^2);      % trajectory projected onto y axis
%     uz0         = -sin_theta_in*cos(phi_inc*pi/180)*sqrt(1-uz^2) - uz*sqrt(1-sin_theta_in^2); % such that ux^2 + uy^2 + uz^2 = 1
% else
%     ux0=0;
%     uy0=0;
%     uz0=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%
% Prepare Monte Carlo
%%%
nm = 440;
% Create tissue properties
tissueX = makeTissueList_pan(); % also --> global tissue(1:Nt).s %Tissue list at excitation
% tissueX = makeTissueList2();
Nt = length(tissueX);
for i=1:Nt
    muavX(i)  = tissueX(i).mua;
    musvX(i)  = tissueX(i).mus;
    gvX(i)    = tissueX(i).g;
    gammaX(i) = tissueX(i).gamma;
    mua_mult(i)= tissueX(i).mua_mult;
    mus_mult(i) = tissueX(i).mus_mult;
end

InterfaceList = makeXfaceList();
Nofinterfaces = length(InterfaceList);
for i=1:Nofinterfaces
    n_in(i) = InterfaceList(i).n_in;
    n_out(i) = InterfaceList(i).n_out;
end


Fluor=makeFluorList();
Nf = length(Fluor);
for i =1:Nf
    muaf(i)=Fluor(i).muaf; 
    muaf_mult(i) = Fluor(i).muaf_mult;
    Qeff(i)=Fluor(i).quanteff;
    tau(i) = Fluor(i).tau;
end
% Specify Monte Carlo parameters
Nx = Nbins;
Ny = Nbins;
Nz = Nbins;
dx = binsize;
dy = binsize; 
dz = binsize;
x  = ([1:Nx]'-Nx/2+1/2)*dx;
y  = ([1:Ny]'-Ny/2+1/2)*dy;
z  = [[1:Nz]-1]'*dz;
% x  = ([1:Nx]'-Nx/2+1/2)*dx;
% y  = ([1:Ny]'-Ny/2+1/2)*dy;
% z  = [[1:Nz]+1]'*dz;
% 


zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);

if isinf(zfocus), zfocus = 1e12; end

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

% At excitation
TX = double(zeros(Ny,Nx,Nz));
% TM = double(zeros(Ny,Nx,Nz));
TF = double(zeros(Ny,Nx,Nz));
% z_int=0.03;
% for iz=1:Nz
%     if (iz<=round(z_int/dz))
%     TX(:,:,iz) =13 ;      % fill background with standard tissue
%     else
%         TX(:,:,iz)=14;
%     end
% end
%
TX=TX+2;
% TM=TM+16;

%

% At emission
%     TM = double(zeros(Ny,Nx,Nz));
%     TF = double(zeros(Ny,Nx,Nz));
%     ztop=0.02;
%     zbot = 0.03;
%     for iz = 1:Nz
%         if (iz<round(ztop/dz))||(iz>round(zbot/dz))
%                TM(:,:,iz) = 10;  % fill background with standard tissue at emission
%         else
%             TM(:,:,iz)=12;
%             TF(:,:,iz)=1;
%         end
%     end
%
TF = TF+3;

xc1 = 0;  
zc1 = 0.15; %0.625; 
yc1 = 0;
% xc2 = 0; zc2 = 0.1; yc2 = 0;
% xc3 = 0.05; zc3 = 0.15; yc3 = 0;


r1 = 0.075; %r2 = 0.01; r3 = 0.01;
% % % % % 
for iz = 1:Nz
 for ix=1:Nx
     for iy = 1:Ny
            xd = x(ix) - xc1;	
            zd = z(iz) - zc1;
            yd = y(iy) - yc1;
            r  = sqrt(xd^2 + zd^2 + yd^2);	% r from inhomogeneity centre
            if (r<=r1)     	% if r is within inhomogeneity
                TF(iy,ix,iz) = 4; % fluorophore
%                 TX(iy,ix,iz) = 22;
            end
            
%             xd = x(ix) - xc2;	
%             zd = z(iz) - zc2;
%             yd = y(iy) - yc2;
%             r  = sqrt(xd^2 + zd^2 + yd^2);	% r from inhomogeneity centre
%             if (r<=r3)     	% if r is within inhomogeneity
%                 TF(iy,ix,iz) = 10; % fluorophore
%             end
            
%             xd = x(ix) - xc3;	
%             zd = z(iz) - zc3;
%             yd = y(iy) - yc3;
%             r  = sqrt(xd^2 + zd^2 + yd^2);	% r from inhomogeneity centre
%             if (r<=r3)     	% if r is within inhomogeneity
%                 TF(iy,ix,iz) = 10; % fluorophore
%             end

% %             r  = sqrt(xd^2 + zd^2);	% r from vessel center
%             if (r<=r2)     	% if r is within vessel
%                 TF(iy,ix,iz) = 11; % blood
%             end
     end % iy
 end %ix
end %iz


%zsurf = 0.0100;  % position of air/skin surface

% for iz=1:Nz % for every depth z(iz)

%     % air
%     if iz<=round(zsurf/dz)
%         T(:,:,iz) = 2;
%     end

% epidermis (60 um thick)
%     if iz>round(zsurf/dz) & iz<=round((zsurf+0.0060)/dz)
%         T(:,:,iz) = 5;
%     end
%
%     % blood vessel @ xc, zc, radius, oriented along y axis
%     xc      = 0;            % [cm], center of blood vessel
%     zc      = Nz/2*dz;     	% [cm], center of blood vessel
%     vesselradius  = 0.0100;      	% blood vessel radius [cm]
%     for ix=1:Nx
%             xd = x(ix) - xc;	% vessel, x distance from vessel center
%             zd = z(iz) - zc;   	% vessel, z distance from vessel center
%             r  = sqrt(xd^2 + zd^2);	% r from vessel center
%             if (r<=vesselradius)     	% if r is within vessel
%                 T(:,ix,iz) = 3; % blood
%             end
%
%     end %ix

% end % iz
%%
DetectorLoc=[-0.15,-0.15,0
    -0.1,-0.15,0
    0,-0.15,0
    0.1,-0.15,0
    0.15,-0.15,0
    -0.15,-0.1,0
    -0.05,-0.1,0
    0.05,-0.1,0
    0.15,-0.1,0
    -0.1,-0.05,0
    0,-0.05,0
    0.1,-0.05,0
    -0.15,0,0
    -0.05,0,0
    0.05,0,0
    0.15,0,0
    -0.15,0.15,0
    -0.1,0.15,0
    0,0.15,0
    0.1,0.15,0
    0.15,0.15,0
    -0.15,0.1,0
    -0.05,0.1,0
    0.05,0.1,0
    0.15,0.1,0
    -0.1,0.05,0
    0,0.05,0
    0.1,0.05,0];

NofDetectors = size(DetectorLoc, 1);
Detfiber_dia = 0.0125*ones(NofDetectors,1);
det_ang_min =  40;
det_ang_max =  61;



%% Create interfaces
Xface = double(3*ones(Ny,Nx,Nz));
% Check interfaces and apply appropriate values 
    Xface(:,:,1) = 1; % Only top surface has different R.I. in this case
%     Xface(:,:,end) =2;
%     Xface(:,1,:)=2;
%     Xface(:,end,:)=2;
%     Xface(1,:,:) = 2;
%     Xface(end,:,:) = 2;
for ii=1:NofDetectors
    for ix = 1:Nx
        for iy = 1:Ny
            r = ((x(ix) - DetectorLoc(ii,1))^2 + (y(iy) - DetectorLoc(ii,2))^2) - Detfiber_dia(ii)^2;
            if (r<=0)
                Xface(ix,iy,1) = 6;
            end
        end
    end
end


%     
%%
% OUTPUT FLAGS
eval_diffuse_fluence =1; % Evaluate diffuse fluence, set to zero if only total fluence is requires
Nscatbins = 1; % #  of scattering bins to evaluate fluence. Will store fluence > n1, >n2,  >n3... scattering events respectively
Nscatevents = 1; %[n1,n2,n3...]
eval_partcurr = 1; % 1 = evaluate partial current. 0 = total diffuse reflectance.

%%
if SAVEON
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    vx = uint8(reshape(TX,Ny*Nx*Nz,1));
%     vm = uint8(reshape(TM,Ny*Nx*Nz,1));
    vf = uint8(reshape(TF,Ny*Nx*Nz,1));
    vxface = uint8(reshape(Xface,Ny*Nx*Nz,1));
    
    %% WRITE FILES
    % Write myname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
    % run parameters
    fprintf(fid,'%ld\n',Nphotons);
    fprintf(fid,'%d\n'   ,Nx);
    fprintf(fid,'%d\n'   ,Ny);
    fprintf(fid,'%d\n'   ,Nz);
    fprintf(fid,'%0.4f\n',dx);
    fprintf(fid,'%0.4f\n',dy);
    fprintf(fid,'%0.4f\n',dz);
    % launch parameters
    fprintf(fid,'%d\n'   ,mcflag);
    fprintf(fid,'%d\n'   ,launchflag);
    fprintf(fid,'%d\n'   ,boundaryflag);
    fprintf(fid,'%d\n'   ,emflag);
    fprintf(fid,'%0.4f\n',xs);
    fprintf(fid,'%0.4f\n',ys);
    fprintf(fid,'%0.4f\n',zs);
    fprintf(fid,'%0.4f\n',xfocus);
    fprintf(fid,'%0.4f\n',yfocus);
    fprintf(fid,'%0.4f\n',zfocus);
    fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
    fprintf(fid,'%0.4f\n',uy0);
    fprintf(fid,'%0.4f\n',uz0);
    fprintf(fid,'%0.4f\n',radius);
    fprintf(fid,'%0.4f\n',waist);
    fprintf(fid,'%0.4f\n',mod_freq);
    fprintf(fid,'%d\n',pfuncflag);
    % tissue optical properties
    fprintf(fid,'%d\n',Nt);
    for i=1:Nt
        fprintf(fid,'%0.4f\n',muavX(i));
        fprintf(fid,'%0.4f\n',musvX(i));
        fprintf(fid,'%0.4f\n',gvX(i));
        fprintf(fid,'%1.1f\n',gammaX(i));
        fprintf(fid,'%0.4f\n',mua_mult(i));
        fprintf(fid,'%0.4f\n',mus_mult(i));
    end
    %Refractive Index
    fprintf(fid,'%d\n',Nofinterfaces);
    for i=1:Nofinterfaces
        fprintf(fid,'%0.4f\n',n_in(i));
        fprintf(fid,'%0.4f\n',n_out(i));
    end
    
    %Fluorophore parameters
    fprintf(fid,'%d\n',Nf);
    for i = 1:Nf
        fprintf(fid,'%0.4f\n',muaf(i));
        fprintf(fid,'%0.4f\n',muaf_mult(i));    
        fprintf(fid,'%0.4f\n',Qeff(i));
        fprintf(fid,'%0.4f\n',tau(i));
    end
    
    % Output flags
    fprintf(fid,'%d\n',eval_partcurr);
    fprintf(fid,'%0.4f\n',det_ang_min);
    fprintf(fid,'%0.4f\n',det_ang_max);
    fprintf(fid,'%d\n',NofDetectors); 
    for ii=1:NofDetectors
        fprintf(fid,'%0.4f \n',DetectorLoc(ii,1));
        fprintf(fid,'%0.4f \n',DetectorLoc(ii,2));
        fprintf(fid,'%0.4f \n',Detfiber_dia(ii));
    end
    fclose(fid);
    
    %% write myname_T.bin file
    filename = sprintf('%s_TX.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,vx,'uint8');
    fclose(fid);
    
    filename = sprintf('%s_XFACE.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,vxface,'uint8');
    fclose(fid);
    
    
%     filename = sprintf('%s_TM.bin',myname);
%     disp(['create ' filename])
%     fid = fopen(filename,'wb');
%     fwrite(fid,vm,'uint8');
%     fclose(fid);
%     
    
    filename = sprintf('%s_TF.bin',myname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,vf,'uint8');
    fclose(fid);
    
    
    toc
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(TF,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx

%%
figure(1); clf
sz = 12;  fz = 10;
imagesc(x,z,Tzx,[1 Nf])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
colorbar
cmap = makecmap(Nf);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%

for i=1:Nf
    yy = (Nf-i)/(Nf-1)*Nz*dz;
    text(max(x)*1.2,yy, Fluor(i).name,'fontsize',fz)
end

text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])

%%% draw launch
N = 20; % # of beam rays drawn
% switch mcflag
%     case 0 % uniform
%         for i=0:N
%             plot((xs-radius + 2*radius*i/N)*[1 1],[zs min(z)],'r-')
              plot([xs max(x)],(zs-radius + 2*radius*i/N)*[1 1],'r-')
%         end
%
%     case 1 % Gaussian
%         for i=0:N
%             plot([(-radius + 2*radius*i/N) xfocus],[zs zfocus],'r-')
%         end

%     case 2 % iso-point
%         for i=1:N
%             th = (i-1)/19*2*pi;
%             xx = Nx/2*cos(th) + xs;
%             zz = Nx/2*sin(th) + zs;
%             plot([xs xx],[zs zz],'r-')
%         end
%
%     case 3 % rectangle
%         zz = max(z);
%         for i=1:N
%             xx = -radius + 2*radius*i/20;
%             plot([xx xx],[zs zz],'r-')
%         end
% end

disp('done')

