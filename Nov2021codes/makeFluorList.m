function fluor = makeFluorList()
%% Create tissueList

j=1;
fluor(j).name  = 'Fluorophore 1';
fluor(j).muaf = 0.006; 
fluor(j).muaf_mult = 0.0846; 
fluor(j).quanteff = 0.016;
fluor(j).tau = 0.56; % in ns

j=2;
fluor(j).name  = 'Fluorophore 2';
fluor(j).muaf = 0.08; 
fluor(j).muaf_mult = 1; 
fluor(j).quanteff = 0.019;
fluor(j).tau = 0.56;

j=3;
fluor(j).name  = 'null_try';
fluor(j).muaf = 0; 
fluor(j).muaf_mult = 0; 
fluor(j).quanteff = 0;
fluor(j).tau = 0; % in ns

j=4;
fluor(j).name  = 'Fluorophore 4';
fluor(j).muaf = 0.3; 
fluor(j).muaf_mult = 0.0846; 
fluor(j).quanteff = 0.016;
fluor(j).tau = 0.56; % in ns

j=5;
fluor(j).name  = 'spn_validation';
fluor(j).muaf = 0.147; 
fluor(j).muaf_mult = 1; 
fluor(j).quanteff = 1;
fluor(j).tau = 1; % in ns

j=6;
fluor(j).name  = 'fibtest';
fluor(j).muaf = 1; 
fluor(j).muaf_mult = 0.0846; 
fluor(j).quanteff = 0.016;
fluor(j).tau = 0.56; % in ns

j=7;
fluor(j).name  = 'Fluorophore 4';
fluor(j).muaf = 0.15; 
fluor(j).muaf_mult = 0.0846; 
fluor(j).quanteff = 0.016;
fluor(j).tau = 0.56; % in ns

j=8;
fluor(j).name  = 'Fluorophore 5';
fluor(j).muaf = 0.2; 
fluor(j).muaf_mult = 0.0846; 
fluor(j).quanteff = 0.016;
fluor(j).tau = 0.56; % in ns

j=9;
fluor(j).name  = 'Phantom1_fl';
fluor(j).muaf = 0.005; 
fluor(j).muaf_mult = 0.08; 
fluor(j).quanteff = 1;
fluor(j).tau = 0.56; % in ns

j=10;
fluor(j).name  = 'P1F1';
fluor(j).muaf = 0.5; 
fluor(j).muaf_mult = 1; 
fluor(j).quanteff = 1;
fluor(j).tau = 0.56; % in ns

j=11;
fluor(j).name  = 'Phantom3_fl';
fluor(j).muaf = 0.05; 
fluor(j).muaf_mult = 0.08; 
fluor(j).quanteff = 1;
fluor(j).tau = 0.56; % in ns

j=12;
fluor(j).name  = 'Phantom4_fl';
fluor(j).muaf = 0.5; 
fluor(j).muaf_mult = 0.08; 
fluor(j).quanteff = 1;
fluor(j).tau = 0.56; % in ns


j=13;
fluor(j).name  = 'P1BF';
fluor(j).muaf = 0.005; 
fluor(j).muaf_mult = 1; 
fluor(j).quanteff = 1;
fluor(j).tau = 0.56; % in ns

j=14;
fluor(j).name  = 'C3F';
fluor(j).muaf = 2; 
fluor(j).muaf_mult = 1; 
fluor(j).quanteff = 1;
fluor(j).tau = 0.56; % in ns


j=15;
fluor(j).name  = 'Fluorophore 4';
fluor(j).muaf = 0.1; 
fluor(j).muaf_mult = 0.0846; 
fluor(j).quanteff = 0.016;
fluor(j).tau = 0.56; % in ns



disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(fluor)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,fluor(i).name, fluor(i).muaf,fluor(i).quanteff))
end
disp(' ')

