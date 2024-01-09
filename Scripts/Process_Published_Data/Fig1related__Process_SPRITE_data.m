
%% Process published SPRITE data into maps

%% Load all SPRITE data
sprite_maps = cell(19,1);
for c=1:19 % 3
    disp(['loading SPRITE data for chr',num2str(c),'...']);
    isChr = strcmp(mm10coords.Var1,['chr',num2str(c)]);
    chrStarts = mm10coords.Var2(isChr);
    % dlmread is much faster than readtable.  Maybe readmatrix can compete.
    sprite_chr_20kb = dlmread(['U:\GenomeData\ByPublication\SPRITE_Gutmann2018\GSE114242_mouse_intra_20kb\mouse_chr',num2str(c),'_20kb_all_nover2.txt']); % slow
    interval = floor(chrStarts/20e3);  % note 20 kb is already hard-coded in the processed SPRITE data, though we could re-bin the raw reads.  
    sprite_SEmap = sprite_chr_20kb(interval,interval);
    sprite_maps{c} = sprite_SEmap;
end
%%
cfrac = ContactFrac(esMap(1:376,1:376,:),'threshold',.3);
for c=1:19
fish_maps{c} = cfrac(cisCoords{c},cisCoords{c});
end

% save([saveFolder,'cisMaps_all.mat'],"fish_maps","sprite_maps","gam_maps","hic_maps");
