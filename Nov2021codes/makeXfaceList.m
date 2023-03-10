function xface = makeXfaceList()
%% Create tissueList

j=1;
xface(j).name  = 'Tissue-Air';
xface(j).n_in = 1.37; 
xface(j).n_out = 1;

j=2;
xface(j).name  = 'air-tissue';
xface(j).n_in = 1.33; 
xface(j).n_out = 1;

j=3;
xface(j).name  = 'Tissue-tissue';
xface(j).n_in = 1.37; 
xface(j).n_out = 1.37;

j=4;
xface(j).name  = 'Gardner';
xface(j).n_in = 1.38; 
xface(j).n_out = 1;

j=5;
xface(j).name  = 'tissue-tissue';
xface(j).n_in = 1.33; 
xface(j).n_out = 1.33;

j = 6; 
xface(j).name  = 'tissue - fiber';
xface(j).n_in = 1.37; 
xface(j).n_out = 1.5;




disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(xface)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,xface(i).name, xface(i).n_in,xface(i).n_out))
end
disp(' ')

