
% lookmcxyz.m
%   Looks at myname_F.bin, created by mcxyz.c
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%   Makes figures:
%       myname_tissue.jpg   = tissue structure (shows tissue types)
%       myname_Fzx.jpg      = fluence rate vs z,x
%       myname_Fzy.jpg      = fluence rate vs z,y
%   Uses:
%       myname_H.mci    = input file from maketissue.m
%       myname_T.bin    = tissue input file from maketissue.m
%       myname_F.bin    = fluence rate output from Monte Carlo
%       reportH_mci.m   = lists input parameters in myname_H.mci
%       makecmap.m      = makes colormap for tissue types
%       makec2f.m       = makes colormap for fluence rate
%
%   This example sets myname = 'skinvessel'.
%
% 7/feb/2017, add boundaryflag (see A(10)).
% 1/june/2017 , no major changes, just clean up display outputs.
% Steven L Jacques
home; %clear
format compact
commandwindow

SAVEPICSON = 0;
SAVEMAT = 1;
if SAVEPICSON
    sz = 10; fz = 7; fz2 = 5; % to use savepic.m
else
    sz = 12; fz = 9; fz2 = 7; % for screen display
end

%%%% USER CHOICES <---------- you must specify -----
myname = 'test_broad'; nm = 532;

%%%%

disp(sprintf('------ mcxyz %s -------',myname))

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])

fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
Nphotons = A(1);
Nx = A(2);
Ny = A(3);
Nz = A(4);
dx = A(5);
dy = A(6);
dz = A(7);
mcflag = A(8);
launchflag = A(9);
boundaryflag = A(10);
emflag=A(11);
xs = A(12);
ys = A(13);
zs = A(14);
xfocus = A(15);
yfocus = A(16);
zfocus = A(17);
ux0 = A(18);
uy0 = A(19);
uz0 = A(20);
radius = A(21);
waist = A(22);
mod_freq = A(23);
pfuncflag = A(24);
Nt = A(25);
j = 25;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
    j=j+1;
    gammav(i,1)=A(j);
    j = j+1; 
    mua_multv(i,1) = A(j);
    j = j+1;
    mus_multv(i,1) = A(j);
end
j=j+1;
nxface=A(j); j = j+1;
for i=1:nxface
    n_in(i) = A(j);j=j+1;
    n_out(i) = A(j);j=j+1;
end
% j=j+1;
Nf = A(j);
for i=1:Nf
        j=j+1;
        muafv(i,1) = A(j);j=j+1;
        muaf_multv(i,1) = A(j); j=j+1;
        Qeffv(i,1) = A(j); j=j+1;
        tauv(i,1) = A(j);
end
j=j+1;
outflag = A(j); j=j+1;

%reportHmci(myname)


%% Load Excitation data rate F(y,x,z)
% Load fluence
filename = sprintf('%s_FXreal.bin',myname);
disp(['loading ' filename])
tic
fid = fopen(filename, 'rb');
[Data count] = fread(fid, Ny*Nx*Nz, 'double');
fclose(fid);
toc
FX_real = reshape(Data,Ny,Nx,Nz); % F(y,x,z)

filename = sprintf('%s_FXimag.bin',myname);
disp(['loading ' filename])
tic
fid = fopen(filename, 'rb');
[Data count] = fread(fid, Ny*Nx*Nz, 'double');
fclose(fid);
toc
FX_imag = reshape(Data,Ny,Nx,Nz); % F(y,x,z)


% Load Reflectance F(y,x,z)
filename = sprintf('%s_RXreal.bin',myname);
disp(['loading ' filename])
tic
fid = fopen(filename, 'rb');
[Data count] = fread(fid, Ny*Nx, 'double');
fclose(fid);
toc
RX_real = reshape(Data,Nx,Ny); % F(y,x,z)

filename = sprintf('%s_RXimag.bin',myname);
disp(['loading ' filename])
tic
fid = fopen(filename, 'rb');
[Data count] = fread(fid, Ny*Nx, 'double');
fclose(fid);
toc
RX_imag = reshape(Data,Nx,Ny); % F(y,x,z)


%%
% Load tissue structure in voxels, T(y,x,z)
filename = sprintf('%s_TX.bin',myname);
disp(['loading ' filename])
tic
fid = fopen(filename, 'rb');
[Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
fclose(fid);
toc
TX = reshape(Data,Ny,Nx,Nz); % T(y,x,z)

if (emflag)
filename = sprintf('%s_TF.bin',myname);
disp(['loading ' filename])
tic
fid = fopen(filename, 'rb');
[Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
fclose(fid);
toc
TF = reshape(Data,Ny,Nx,Nz);

clear Data
end
%%
x = ([1:Nx]-Nx/2+1/2)*dx;
y = ([1:Ny]-Ny/2+1/2)*dx;
z = ([1:Nz]-1)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

%% Look at structure, Tzx
Tzx = reshape(TX(Ny/2,:,:),Nx,Nz)';
tissue = makeTissueList_pan();
Nt = length(tissue);

figure();
imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',sz)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
title('Tissue','fontweight','normal','fontsize',fz2)
for i=1:Nt
    yy = zmin + (Nt-i)/(Nt-1)*zdiff;
    text(xmax*1.4,yy, sprintf('%d %s',i,tissue(i).name),'fontsize',fz2)
end

if SAVEPICSON
    name = sprintf('%s_tissue.jpg',myname);
    savepic(1,[4 3],name)
end


%% Look at Fluence Fzx @ launch point
FX = FX_real + 1i*FX_imag;

%     mua_cplx = muav(TX) + muafv(TF) + 1i*(2*pi*mod_freq*1.33/(3*10^4));

%     FX = FX.*mua_cplx;
    
Fzx = squeeze(FX(Ny/2+1,:,:))'; % in z,x plane through source

figure(100);
subplot(2,2,1)
imagesc(x(ux),z(uz),(log10(abs(Fzx(uz,ux)))))
colorbar
xlabel('x [cm]')
ylabel('z [cm]')
title('Magnitude -\Phi_x')

subplot(2,2,3)
imagesc(x(ux),z(uz),(angle(Fzx(uz,ux)))*180/pi)
colorbar
xlabel('x [cm]')
ylabel('z [cm]')
title('Phase -\Phi_x')



if SAVEPICSON
    name = sprintf('%s_Fzx.jpg',myname);
    savepic(2,[4 3],name)
end

%% look Fzy
Fzy = squeeze(FX(:,Nx/2,:))';
subplot(2,2,2)
imagesc(y(uy),z(uz),real(log10(Fzy(uz,uy))))
colorbar
xlabel('y [cm]')
ylabel('z [cm]')
title('Magnitude-\Phi_x')

subplot(2,2,4)
imagesc(y(uy),z(uz),imag(log10(Fzy(uz,uy))))
colorbar
xlabel('y [cm]')
ylabel('z [cm]')
title('Phase \Phi_x')


if SAVEPICSON
    name = sprintf('%s_Fzy.jpg',myname);
    savepic(3,[4 3],name)
end
%%
figure(101)
RX = RX_real + 1i*RX_imag;
subplot(1,2,1)
imagesc(x(ux),y(uy),real(log10(RX(ux,uy))))
colorbar
xlabel('x [cm]')
ylabel('y [cm]')
title('Magnitude - RX')

subplot(1,2,2)
imagesc(x(ux),y(uy),imag(log10(RX(ux,uy))))
colorbar
xlabel('x [cm]')
ylabel('y [cm]')
title('Phase - RX')


%% Emission

if (emflag)
    
    
    %% Load Emission Data
    filename = sprintf('%s_FMreal.bin',myname);
    disp(['loading ' filename])
    tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'double');
    fclose(fid);
    toc
    FM_real = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
    filename = sprintf('%s_FMimag.bin',myname);
    disp(['loading ' filename])
    tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'double');
    fclose(fid);
    toc
    FM_imag = reshape(Data,Ny,Nx,Nz); % F(y,x,z)
        
    filename = sprintf('%s_RMreal.bin',myname);
    disp(['loading ' filename])
    tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx, 'double');
    fclose(fid);
    toc
    RM_real = reshape(Data,Ny,Nx); % F(y,x,z)
     filename = sprintf('%s_RMimag.bin',myname);
    disp(['loading ' filename])
    tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx, 'double');
    fclose(fid);
    toc
    RM_imag = reshape(Data,Ny,Nx); % F(y,x,z)
    
    %%
    FM = FM_real + 1i*FM_imag;
    
    
    
%     mua_cplx = muav(TX).*mua_multv(TF) + muafv(TF).*muaf_multv(TF) - 1i*(2*pi*mod_freq*1.33/(3*10^4));
    
%     FM = FM.*mua_cplx;
      
    
    
    figure(200)
    Fzx = squeeze(FM(Ny/2,:,:))'; % in z,x plane through source
    subplot(2,2,1)
    imagesc(x(ux),z(uz),real(log10(Fzx(uz,ux))))
    xlabel('x [cm]')
    ylabel('z [cm]')
    title('Magnitude \Phi_m')
    
    subplot(2,2,3)
    imagesc(x(ux),z(uz),imag(log10(Fzx(uz,ux)))*180/pi)
    xlabel('x [cm]')
    ylabel('z [cm]')
    title('Phase \Phi_m')
    
    
    Fzy = squeeze(FM(:,Nx/2,:))';
    subplot(2,2,2)
    imagesc(y(uy),z(uz),real(log10(Fzy(uz,uy))))
    xlabel('y [cm]')
    ylabel('z [cm]')
    title('Magnitude \Phi_m')
    
    subplot(2,2,4)
    imagesc(y(uy),z(uz),imag(log10(Fzy(uz,uy)))*180/pi)
    xlabel('y [cm]')
    ylabel('z [cm]')
    title('Phase \Phi_m')
    
    %%
    figure (250)
    RM = RM_real + 1i*RM_imag;
    subplot(1,2,1)
    imagesc(x(ux),y(uy),real(log10(RM(ux,uy))));
    xlabel('x [cm]')
    ylabel('y [cm]')
    title('Magnitude Rm')
    
    subplot(1,2,2)
    imagesc(x(ux),y(uy),imag(log10(RM(ux,uy))))
    xlabel('x [cm]')
    ylabel('y [cm]')
    title('Phase Rm')
    
    
    %%
    
    if SAVEPICSON
        name = sprintf('%s_diffzy.jpg',myname);
        savepic(3,[4 3],name)
    end
    
    
end
disp('done')

if (SAVEMAT)
    month_now = month(datetime,'shortname');
    year_now = year(datetime);
fname_save = [myname,'_mc_',month_now{1},num2str(year_now),'.mat'];
save(fname_save,'dx','dy','dz','FM','FX','gammav','gv','boundaryflag','launchflag','mcflag','mod_freq','mua_multv','muaf_multv',...
    'muafv','muav','mus_multv','musv','myname','n_in','n_out','Nf','Nphotons','Nt','Nx','Ny','Nz','outflag','pfuncflag','Qeffv','radius',...
    'RM','RX','tauv','tissue','ux0','uy0','uz0','x','y','z','xs','ys','zs');
end