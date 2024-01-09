%% Process GAM data into maps

%%  GAM test load
% Note, GAM data is mm9 not mm10.  we will map it to the mm9 SEs. 
GenomeFolder = 'U:\';
NAS02_Vol4 = '\\169.254.113.81\NAS02_Vol4\';
c= 19;
gam_chr = readtable([GenomeFolder,'GenomeData\ByPublication\GAM_Pombo2017\GSE64881_chr',num2str(c),'_chr',num2str(c),'.1Mb.txt']);
disp(gam_chr(1,:))
saveFolder =  [NAS02_Vol4,'Derek\SE_analyses_supplementary_files\'];

%% get GAM map for our SEs
gamStep = 30e3; % GAM bin size
gam_maps = cell(19,1);
fish_maps = cell(19,1);
for c=1:nChr  % c=19
    isChr = strcmp(mm9coords.Var1,['chr',num2str(c)]);
    chrStarts = mm9coords.Var2(isChr);
    gam_chr_30kb = readtable([GenomeFolder,'GenomeData\ByPublication\GAM_Pombo2017\GSE64881_chr',num2str(c),'_chr',num2str(c),'.30kb.txt']);
    interval = floor(chrStarts/gamStep);
    gam_SEmap = gam_chr_30kb{interval,interval};
    gam_maps{c} = gam_SEmap;
end
