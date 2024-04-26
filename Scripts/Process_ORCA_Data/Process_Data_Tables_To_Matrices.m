
%% covert data tables per FOV to a stack of single cell matrices 

% CELL BARCODES  (grp 1 to 4)
% E14-mESC	381 =       2
% JB-mESC   382 =       4
% JB-EpiLC	383 =       1
% JB-NPC    384 =       3
% order imaged 383, 381, 384, 382

names = {'EpiLC-JB','mESC-E14','NPC-JB','mESC-JB'};
% EpiLC=2223, mESC-E14=10733, NPC-JB=1888, mESC-JB=4097

saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\ForZenodo\Precomputed_Analysis_Files\'
analysisFolder = 'U:\Manuscripts\SE Clustering Paper\Data\ForZenodo\Corrected_Data_Tables_by_FOV\'
%% Load the data 
% this part is slow
% sort it into single cell distance maps and individual polymers
tic

datTables = FindFiles([analysisFolder,'FOV*_CorrChromDriftGrp.csv']);
nFOV = length(datTables);
chrTableHyb = readtable([saveFolder,'SuperEnhancerLoci.xlsx']) ; % mm10   

% faster to  load 1-fov at a time
%   indexing out of a crazy massive table is unnecessarily slow
%   this is also more memory efficient
nB = 376*2;  % 752x752x20K breaks the computer
allMaps = cell(nFOV,1); 
allPols = cell(nFOV,1); 
cellTypes = cell(nFOV,1);
for f=1:nFOV % f = 1
    disp(['processing FOV ',num2str(f),' of ',num2str(nFOV)]);
    spotTable = readtable(datTables{f});
    nCells = max(spotTable.cellID);
    nTraces = nCells; 
    pols = nan(nB,3,nTraces);
    maps = nan(nB,nB,nTraces,'single');
    cellType = nan(nTraces,1);
    for c=1:nCells % c=10
        cellTable = spotTable(spotTable.cellID == c,:); 
            traceTable0 = cellTable(cellTable.allele==0,:);
            xyz0 = [traceTable0.x,traceTable0.y,traceTable0.z];
            b0 = traceTable0.genomeID;
            traceTable1 = cellTable(cellTable.allele==1,:);
            xyz1 = [traceTable1.x,traceTable1.y,traceTable1.z];
            b1 = nB/2 + traceTable1.genomeID;
            xyz = cat(1,xyz0,xyz1);
            b = cat(1,b0,b1);
            traceTable = cat(1,traceTable0,traceTable1);
            if size(xyz,1)>1 % need at least 2 points for a trace
                cellType(c,1) = traceTable.groupID(1);
                pols(b,:,c) = xyz;
                maps(b,b,c) = squareform(pdist(xyz)); % the SLOW step
            end
    end
    allMaps{f} = maps;
    allPols{f} = pols;
    cellTypes{f} = cellType;
end
% combining all the data
distMaps = cat(3,allMaps{:}); % this is the unfeasible 
clear allMaps; % can't afford 2 copies of this.
polTraces = cat(3,allPols{:});
cellGrps = cat(1,cellTypes{:});
toc
%% gene names
nearestGene = chrTableHyb.Yo_gene;
nearestAllele = [nearestGene; nearestGene];
idx_pp = contains(nearestGene,{'Sox2','Mycn','Tbx3','Pou5f1','Pvt1','Nanog','Prdm14'});

%% divide polTraces into groups
polTraceCell = cell(4,1);
for g=1:4
    polTraceCell{g} = polTraces(:,:,cellGrps==g);
end


%% divide data into groups
% run once! (variables cleared at end
% this is very RAM heavy, it is oscillitory though

if exist('distMaps','var')
    figure(11); clf; figure(10); clf; 
    cFrac2 = cell(4,1);
    cFrac3 = cell(4,1);
    cFrac6 = cell(4,1);
    nObs = cell(4,1);
    medDist = cell(4,1);
    distMapCell = cell(4,1);
    for g=1:4
        distMapCell{g} = distMaps(:,:,cellGrps==g);
        medDist{g} = nanmedian(distMaps(:,:,cellGrps==g),3);
         [cFrac2{g},nObs{g}] = ContactFrac(distMaps(:,:,cellGrps==g),'threshold',.2);
         [cFrac3{g},nObs{g}] = ContactFrac(distMaps(:,:,cellGrps==g),'threshold',.3);
         [cFrac6{g},nObs{g}] = ContactFrac(distMaps(:,:,cellGrps==g),'threshold',.6);
        figure(11); subplot(2,2,g); imagesc(medDist{g}); colorbar; title([names{g},'  n=',num2str(sum(cellGrps==g))]);
        figure(10); subplot(2,2,g); imagesc(nObs{g}); colorbar; title(names{g});
    end
     clear distMaps; % run once
end

%%
cFrac = cFrac2; % 200 nm cutoff 
meanObs = nanmean(nObs{2});
meanObs = meanObs(1:376);  % validated already the two halves are identical in badhybes  
bH = [209:212 ,find(meanObs < .5*max(meanObs))];  % 209:212 had an error in the barcode updates  
bH = [bH, 48, 337]; % 
bH = [bH,bH+376];
% these bad hybes have low detection efficiency and often as not record
% noise instead of true signal. To avoid the resulting confusion, we simply
% exclude these SE from furtehr analysis. 
cxr = [91,213,231,257,291,296,333]; 
cxr = [cxr,cxr+376]; % both alleles;
%  remove bad hybes from single cell data

g = 2; % take just the E14 cells
esMap = distMapCell{g};
esMap(bH,:,:) = nan;
esMap(:,bH,:) = nan;
% remove cross reactive hybes, but not on main diagonal 
idx = sub2ind([2*376,2*376],cxr,cxr);
nC = size(esMap,3);
for c=1:nC  
    temp = esMap(:,:,c);
    temp(cxr,cxr) = nan;
    temp(idx) = 0;
    esMap(:,:,c) = temp; 
end
