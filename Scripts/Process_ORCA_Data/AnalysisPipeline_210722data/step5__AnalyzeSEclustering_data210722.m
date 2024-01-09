%%
% CELL BARCODES  (grp 1 to 4)
% E14-mESC	381
% JB-mESC		382
% JB-EpiLC	383
% JB-NPC		384


%% Load the data 
% this part is slow
% sort it into single cell distance maps and individual polymers

tic
analysisFolder = 'J:\Derek_TEMP\20210722_L10_SE_MultiplexedCells\Analysis_CT4_v15\';
datTables = FindFiles([analysisFolder,'FOV*_CorrChromDriftGrp.csv']);
nFOV = length(datTables);

% faster to  load 1-fov at a time
%   indexing out of a crazy massive table is unnecessarily slow
%   this is also more memory efficient

allMaps = cell(nFOV,1); 
allPols = cell(nFOV,1); 
cellTypes = cell(nFOV,1);
for f=1:nFOV % f = 1
    spotTable = readtable(datTables{f});
    nCells = max(spotTable.cellID);
    nTraces = 2*nCells;
    nB = max(spotTable.genomeID);
    pols = nan(nB,3,nTraces);
    maps = nan(nB,nB,nTraces);
    cellType = nan(nTraces,1);
    for c=1:nCells % c=10
        cellTable = spotTable(spotTable.cellID == c,:); 
        for a = 0:1
            traceTable = cellTable(cellTable.allele==a,:);
            xyz = [traceTable.x,traceTable.y,traceTable.z];
            b = traceTable.genomeID;
            k = c+a*nCells;
            if size(xyz,1)>1 % need at least 2 points for a trace
                cellType(k,1) = traceTable.groupID(1);
                pols(b,:,k) = xyz;
                maps(b,b,k) = squareform(pdist(xyz)); 
            end
        end
    end
    allMaps{f} = maps;
    allPols{f} = pols;
    cellTypes{f} = cellType;
end

% combining all the data
distMaps = cat(3,allMaps{:});
cellGrps = cat(1,cellTypes{:});
toc
%% Median distance by cell type
medDist = cell(4,1);
figure(3); clf; 
for g=1:4
    subplot(2,2,g);
    medDist{g} = nanmedian(distMaps(:,:,cellGrps==g),3);
    imagesc(nanmedian(distMaps(:,:,cellGrps==g),3)); 
    colorbar;
    title(sum(cellGrps==g));
end
%% compare cell groups
% difference in median distance
figure(4); clf; 
for g=1:4
    subplot(2,2,g);
    imagesc(medDist{g}-medDist{1}); caxis([-.5,.5]); 
    colorbar;
    title(sum(cellGrps==g));
end
GetColorMap('RedWhiteBlue');

figure(4); clf;
for g=1:4
subplot(2,2,g);
    PlotCorr(medDist{g}(:),medDist{2}(:),'MarkerSize',1); 
end

% exclude the on-diagonal values
figure(5); clf;
for g=1:4
subplot(2,2,g);
    cmap1 = triu(medDist{2},100);
    cmapG = triu(medDist{g},100);
    PlotCorr(cmap1(:),cmapG(:),'MarkerSize',1,'log',true); 
end


% that's nice, the two ESC correlate better than the others.
% let's normalize off the distance effects

%%
cFrac = cell(4,1);
distMapCell = cell(4,1);
figure(3); clf; 
for g=1:4
    subplot(2,2,g);
    distMapCell{g} = distMaps(:,:,cellGrps==g);
    cFrac{g} = ContactFrac(distMapCell{g},'threshold',.5);
    imagesc(log2(cFrac{g})); 
    colorbar;
    title(sum(cellGrps==g));
end
%% contact freq
figure(4); clf; 
for g=1:4
    subplot(2,2,g);
    imagesc(log2(cFrac{g}./cFrac{2})); caxis([-2,2]); 
    colorbar;
    title(sum(cellGrps==g));
end
GetColorMap('RedWhiteBlue');


%%
medNorm = cell(4,1);
figure(5); clf;
for g=1:4
    normM = NormMap(medDist{g});
    medNorm{g} = medDist{g}./normM;
    figure(4); subplot(2,2,g);
    imagesc(medNorm{g});
    figure(5); subplot(2,2,g);
    PlotCorr(medNorm{g}(:),medNorm{1}(:),'MarkerSize',1); 
end

medNorm = cell(4,1);
figure(5); clf;
for g=1:4
    normM = nanmedian(medDist{g}(:));
    medNorm{g} = medDist{g}./normM;
    figure(4); subplot(2,2,g);
    imagesc(medNorm{g});
    figure(5); subplot(2,2,g);
    PlotCorr(medNorm{g}(:),medNorm{1}(:),'MarkerSize',1); 
end



figure(6); clf; 
for g=1:4
    subplot(2,2,g);
    imagesc(medNorm{g}-medNorm{1}); caxis([-.1,.1]); 
    colorbar;
    title(sum(cellGrps==g));
end
GetColorMap('RedWhiteBlue');


%% compare to random simulations? 
% I think we need a simple linear separation decay plot



%% normalizing distance 
% power law doesn't really get things correct at chromosomal scales due to
% the flattening effect of chromosomal territories.  Should really fit a
% sigmoid / saturation curve to the distance distribution data to do this
% correctly. 
%
% the relative to the NPC state where they aren't SEs is perhaps the best
% comparison
chrTableHyb = readtable([analysisFolder,'chrTableHyb.csv']);

% I think we need a simple linear separation decay plot
chrNums = unique(chrTableHyb.chrNum,'stable')

normDist = medDist;
for c=1:length(chrNums) % c=1
    isC = chrTableHyb.chrNum==c;
    if sum(isC)>0
        dists = chrTableHyb.start(isC);
        chrNorm = squareform(pdist(dists)).^-.9;
        figure(1); clf; imagesc(log10(chrNorm)); colormap('default'); colorbar;
        for g=1:4
            lowDecile = quantile(chrNorm(~isnan(chrNorm(:))),.1);
            normDist{g}(isC,~isC) = normDist{g}(isC,~isC)*lowDecile;
            normDist{g}(isC,isC) =  medDist{g}(isC,isC).*chrNorm;
        end
    end
end
figure(6); clf; 
for g=1:4
    subplot(2,2,g);
    imagesc(log10(normDist{g}));  colorbar; caxis([-7,-5]);
end

GetColorMap('default');


%% multiway contacts
scMap = cell(4,1);
for g=1:4
    scMap{g} =double(distMapCell{g}(:,:,:) < .3);
    blnk = isnan(distMapCell{g});
    scMap{g}(blnk) = nan;
end

clc;
figure(5); clf;
aveContact = cell(4,1);
for g=1:4
    contactCnt = squeeze(nansum(scMap{g},1)); % genes x cells
    nanmedian(contactCnt);
    figure(5); subplot(2,2,g); semilogy(1,1,'w.'); hold on; boxplot(contactCnt'); ylim([0 100]);
    aveContact{g} = nanmean(contactCnt,2); % ave over cells
    isHub = find(quantile(contactCnt,.99,2)>5);
    % isHub = find(aveContact{g}>3);
    chrTableHyb(isHub,:)
end

figure(6); clf;
for g=1:4
    subplot(2,2,g);
    rat = aveContact{g}./aveContact{2};
    rat(isinf(rat))  = nan;
    semilogy(rat); ylim([.5,1.5])
end
%%
% figure(4); clf; bar(aveContact);  ylabel('number of contacts'); xlabel('gene ID');
% figure(5); clf; imagesc(contactCnt(:,1:600));

sum(aveContact>1)
sum(aveContact>2)
sum(aveContact>3)
clc;

%%

isHub = contactCnt > 5;
contactCntCond = contactCnt;
contactCntCond(~isHub) = nan;
aveCondContact = nanmean(contactCntCond,2);
figure(5); clf; bar(aveCondContact);  ylabel('number of contacts'); xlabel('gene ID')
medContact = nanmedian(contactCnt,3);
% geneNames(medContact > 1)
figure(3); clf; bar(medContact);  ylabel('number of contacts'); xlabel('gene ID')


figure(3); clf; boxplot(contactCnt','PlotStyle','compact','OutlierSize',1);  
ylim([0,50]);
ylabel('number of contacts');
% set(gca,'XTick',1:nG,'XTickLabel',geneNames,'FontSize',6);
xtickangle(gca,80);  xlabel('')
%     
%     figure(3); clf; violin(squeeze(contactCnt),'bandwidth',1);
%     set(gca, 'YScale', 'log')

contactCnt % numContacts genes x cells

g=4;
hubSize = 10;
isGeneGhub = contactCnt(1,g,:) > hubSize;
figure(4); clf; bar(nansum(scMap(g,:,isGeneGhub),3));
title(['partner freq of SE ',num2str(g),' in hub > ',num2str(hubSize),' of N hubs = ',num2str(sum(isGeneGhub))]);

figure(2); clf; imagesc(ContactFrac(distMaps(:,:,isGeneGhub),'threshold',.5)); title(sum(isGeneGhub));



%% questions
% are there seedhbs? -- any associations that predict the presence 


