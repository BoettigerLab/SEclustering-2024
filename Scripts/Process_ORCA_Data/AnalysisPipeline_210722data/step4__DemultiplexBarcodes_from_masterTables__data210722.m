%% Demultiplex data.  
% 
% The output data is now finally drift corrected, chromatically corrected,
% demultiplexed, and grouped. We name these tables _CorrChromDriftGrp.csv
% 
% This script is written as applied to the 210722 data. 
% 
%% Demulitplex Derek Data

% Find the relevant data:
%  1) analysis folder
%  2) fov folders with cell segmentation maps from cellpose
%  3) drift correction data
%  4) chromatic correction data
%  5) the barcode folders
%  6) the master tables 
nB = 4; % 4 total barcode hybes in this experiment, 2 per folder
analysisFolder = 'J:\Derek_TEMP\20210722_L10_SE_MultiplexedCells\Analysis_CT4_v15\';
cellIDfolders = FindFiles([analysisFolder,'fov*'],'onlyFolders',true);
driftCorrectFiles = FindFiles([analysisFolder,'alignTable_fov','*.csv']);
load([analysisFolder,'tform3D_750to647.mat'],'tform3D'); % chromatic corrections
barFolders = {...
    'J:\Derek_TEMP\20210722_L10_SE_MultiplexedCells\DNA_Expt\barcode_381-383\'
    'J:\Derek_TEMP\20210722_L10_SE_MultiplexedCells\DNA_Expt\barcode_382-384\'};
spotTables = FindFiles([analysisFolder,'Fits/fov*','_masterTable.csv']);

SetFigureSavePath(analysisFolder); 

%% loop over fovs
nFOV = length(cellIDfolders); 
cellBarVs = cell(nFOV,1);
newTables = cell(nFOV,1);
for f=1:nFOV
    spotTable = readtable(spotTables{f});
    % we use the existing cell boundary maps, we just need to load them
    cellIDs = imread([cellIDfolders{f},'fov',num2str(f,'%03d'),'_cellIDs.png']);
    figure(1); clf; imagesc(cellIDs); GetColorMap('distColorsW');

    % Load Barcodes 1,2
    bH =1; % no drift correct for first Hyb (it's not even in the table)
    imBar= LoadDax([barFolders{bH},'ConvZscan_',num2str(f-1,'%02d'),'.dax'],...
                        'maxProject',true,'verbose',false);
    imBars = imBar(:,:,1:2);
    %Load Barcodes 3,4 and apply drift correction
    bH = 2; % second and on
    driftTable = readtable(driftCorrectFiles{f});
    driftCorrData = table2struct(driftTable(driftTable.hyb == bH,:));
    [imBar,imProps] = LoadDax([barFolders{bH},'ConvZscan_',num2str(f-1,'%02d'),'.dax'],...
                        'maxProject',true,'verbose',false);
    imBars(:,:,3)= ApplyReg(imBar(:,:,1),driftCorrData);
    imBars(:,:,4)= ApplyReg(imBar(:,:,2),driftCorrData);

% imaged in order aquired, which is not at all numerical order:
% bar-hyb1  750-B383 + 647-B381, 
% bar-hyb2  750-B384,  647-B382
% order imaged 383, 381, 384, 382

    %% loop over cells per FOV, add cell barcode label
    % barNorm = quantile( cat(3,cellBarVs{:}),.95);
    barNorm =  [0.8428    1.7490    0.4110    3.2671];
    
    nCells = max(cellIDs(:));
    nB = size(imBars,3); 
    spotTable.groupID = nan(height(spotTable),1);
    spotTable.barValues = nan(height(spotTable),nB);
    cellBarVs{f} = zeros(nCells,nB);
    for c=1:nCells
        barValues = zeros(1,nB);
        for b=1:nB
            imBarB = imBars(:,:,b);
            barVs = imBarB(cellIDs==c);
            barValues(b) = median(barVs(:));
        end
        [~,mi] = max(barValues./barNorm);
        cellBarVs{f}(c,:) = barValues;
        isC = spotTable.cellID==c;
        spotTable.barValues(isC,:) = repmat(barValues,sum(isC),1);   
        spotTable.groupID(isC) = mi;   
    end

    [~,cellIdx] =unique(spotTable.cellID,'stable');
    spotTable.groupID(cellIdx);

    % show labels
    xy = [spotTable.x_fid(cellIdx),spotTable.y_fid(cellIdx)]./imProps.xy2um;
    imBars2 = imBars;
    for b=1:nB
        imBars2(:,:,b) = imBars(:,:,b)./barNorm(b);
    end
    f1 = figure(1); clf; Ncolor(5*imBars2);
    % Ncolor(IncreaseContrast(imBars));
    hold on; plot(xy(:,1),xy(:,2),'w.'); 
    grp = cellstr(num2str(spotTable.groupID(cellIdx)));
    text(xy(:,1),xy(:,2),grp,'color','w'); 
    SaveFigure(f1,'name',['demultiplex_fov',num2str(f,'%03d')],'formats',{'png'},'overwrite',true);
    
    % save a new table with a subset of columns
    %    this has the chormatic corrected, allele matched, and now
    %    preliminary grouped values (could normalize and correct better).
    keepCols = {'x','y','z','spt_brightness','spt_bkd','fov','hyb','stageX','stageY','barID','chrID','genomeID','cellID','allele','groupID','barValues'};
    keep = StringFind(spotTable.Properties.VariableNames,keepCols,'exactly',true);
    spotTable= spotTable(:,keep);
    noData = isnan(spotTable.x); % remove NaNs
    spotTable(noData,:) = [];
    newTables{f} = spotTable;    
    saveName = [analysisFolder,'FOV',num2str(f,'%03d'),'_CorrChromDriftGrp.csv'];
    writetable(newTables{f},saveName);
    disp(['wrote ',saveName]);
end


