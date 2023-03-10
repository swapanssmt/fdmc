function tissue = makeTissueList_pan()
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm:
%       tissueProps = [mua; mus; g]';
%       global tissuenames(i).s
%   Uses 
%       SpectralLIB.mat

j=1;
tissue(j).name  = 'standard tissue';
tissue(j).mua   = 1;
tissue(j).mus   = 100;
tissue(j).g     = 0.80;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0;
tissue(j).mus_mult = 0;

j=2;
tissue(j).name  = 'Fedele2';
tissue(j).mua   = 0.031;
tissue(j).mus   = 54.75;
tissue(j).g     = 0.8;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.7987;
tissue(j).mus_mult = 0.732;

j=3;
tissue(j).name  = 'spn_validation';
tissue(j).mua   = 0.05;
tissue(j).mus   = 50;
tissue(j).g     = 0.9;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;


j=4;
tissue(j).name  = 'Kienle_test';
tissue(j).mua   = 2;
tissue(j).mus   = 100;
tissue(j).g     = 0.9;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0;
tissue(j).mus_mult = 0;

j=5;
tissue(j).name  = 'Fedele_ex';
tissue(j).mua   = 0.031;
tissue(j).mus   = 54.75;
tissue(j).g     = 0.8;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0;
tissue(j).mus_mult = 0;

j=6;
tissue(j).name  = 'Kim';
tissue(j).mua   = 0.45;
tissue(j).mus   = 10;
tissue(j).g     = 0.8;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;

j=7;
tissue(j).name  = 'standard tissue';
tissue(j).mua   = 0.001;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0;
tissue(j).mus_mult = 0;

j = 8; 
tissue(j).name  = 'Gardner2';
tissue(j).mua   = 8;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0;
tissue(j).mus_mult = 0;

j = 9;
tissue(j).name  = 'Phantom1';
tissue(j).mua   = 0.1;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;



j = 10;
tissue(j).name  = 'P1';
tissue(j).mua   = 0.1;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;


j = 11;
tissue(j).name  = 'Phantom3';
tissue(j).mua   = 1;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;

j = 12;
tissue(j).name  = 'Phantom4';
tissue(j).mua   = 5;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;

j = 13;
tissue(j).name  = 'Phantom1A';
tissue(j).mua   = 0.1;
tissue(j).mus   = 50;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;

j = 14;
tissue(j).name  = 'Phantom1B';
tissue(j).mua   = 0.1;
tissue(j).mus   = 20;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;

j = 15;
tissue(j).name  = 'Phantom3A';
tissue(j).mua   = 1;
tissue(j).mus   = 50;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;

j = 16;
tissue(j).name  = 'Phantom3B';
tissue(j).mua   = 1;
tissue(j).mus   = 20;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;

j = 17;
tissue(j).name  = 'Phantom4A';
tissue(j).mua   = 5;
tissue(j).mus   = 50;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 0.8;
tissue(j).mus_mult = 0.8;

j = 18;
tissue(j).name  = 'P1A';
tissue(j).mua   = 0.1;
tissue(j).mus   = 50;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;

j = 19;
tissue(j).name  = 'P1B';
tissue(j).mua   = 0.1;
tissue(j).mus   = 20;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;

j = 20;
tissue(j).name  = 'P1_abs';
tissue(j).mua   = 1;
tissue(j).mus   = 50;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;

j = 21; 
tissue(j).name  = 'lu_test';
tissue(j).mua   = 0.05;
tissue(j).mus   = 50;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;

j = 22;
tissue(j).name  = 'P1_inhom';
tissue(j).mua   = 0.3;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;

j = 23;
tissue(j).name  = 'HSHA';
tissue(j).mua   = 2;
tissue(j).mus   = 100;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;

j = 24;
tissue(j).name  = 'LSHA';
tissue(j).mua   = 2;
tissue(j).mus   = 50;
tissue(j).g     = 0.90;
tissue(j).gamma = 1.8;
tissue(j).mua_mult = 1;
tissue(j).mus_mult = 1;


disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(tissue)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name, tissue(i).mua,tissue(i).mus,tissue(i).g,...
        tissue(i).mus*(1-tissue(i).g)))
end
disp(' ')

