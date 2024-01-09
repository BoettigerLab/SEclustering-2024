%% Process Bonev Hi-C data into maps

%% Load all HiC data
hic_maps = cell(19,1);
mm10coords = readtable([NAS02_Vol4,'Derek\','20210722_L10_SE_MultiplexedCells\Elist_coords_mm10.txt']);
for c=1:19  % up to 8 ran, 9 is missing
    disp(['loading Hi-C data for chr',num2str(c),'...']);
    bonevChrTxt = ['U:\GenomeData\ByPublication\Bonev2017\HiC_Bonev\','GSE161259_ES_chr',num2str(c),'.validPairs.csort.txt'];
    opts = detectImportOptions(bonevChrTxt); % this is super slow;
    %%
    opts.DataLines =[2,Inf];
    opts.SelectedVariableNames = {'Var4','Var8'};
    bonevChrC = readtable(bonevChrTxt,opts)  % 'DataLines',[2,100]
    
    isChr = strcmp(mm10coords.Var1,['chr',num2str(c)]);
    chrSE = mm10coords(isChr,:);
    nB = height(chrSE);
    hicSE = nan(nB,nB);
    for i=1:nB
        for j=1:nB
            read1 = bonevChrC.Var1 > chrSE.Var2(i) & bonevChrC.Var1 < chrSE.Var2(i) + 15e3;
            read2 = bonevChrC.Var2 > chrSE.Var2(j) & bonevChrC.Var2 < chrSE.Var2(j) + 15e3;           
            hicSE(i,j) = sum(read1 & read2);
        end
    end
    % figure(1); clf; imagesc(log10(hicSE)); colorbar; 
    %     colormap(flipud(GetColorMap('redToWhiteK')));
    
    hicSE2 = hicSE + triu(hicSE,1)';
                            hic_maps{c} = hicSE2;
    figure(1); clf; imagesc(log10(hicSE2)); colorbar; 
    colormap(flipud(GetColorMap('redToWhiteK')));
end
