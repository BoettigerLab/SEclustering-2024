%% analysis
% this script ran the analysis on the mini-DNA data set from 2022-04-26
%  which combines DNA FISH for and RNA FISH for the pluripotency gene
%  network. 
% Here we do the drift correction and spot fitting for the DNA data
% still need to chromatic correction. 

% clear all;
% startup;
% clc;

% rerun to correct error in HybName

%% on the fly version
% load experiment layout
t_total = tic; 

expFolder = 'J:\20220426_E14_TwistRNA_seDNA\';
dataFolder = [expFolder,'DNA',filesep];
analysisFolder = [expFolder,'Analysis_CT4_v2',filesep];
hybFolderRoot = ''; % '';


% old experiment table format is now optional
%   Much of the information about the organization of the dax files is now
%   saved in the associated dax.xml, and can be read automatically. 

eTable = readtable([dataFolder,'2022-04-26_E14-mES_TwistRNA_seDNA_eTable_v2.xlsx']);
hybFolders = strcat(dataFolder,eTable.FolderName,filesep);
hybFolderNames = eTable.FolderName;

isH = strcmp(eTable.DataType,'H');
% hybFolders = cat(1,hybFolders{isH},hybFolders(~isH));
% hybFolderNames = cat(1,hybFolderNames{isH},hybFolderNames(~isH));

% HERE we use the toe folders first, and we use the data channel in the toe
% for calling cell boundaries, because the fiducial channel is blurry
hybFolders = cat(1,hybFolders(~isH),hybFolders{isH});
 hybFolderNames = cat(1,hybFolderNames{~isH},hybFolderNames(isH));

% some options
useParallel = true;
overwriteAlignment = false;
maskOutOtherCells = true;
pars.saveSpotsFig = 0; % figure to save SpotFig images, set to zero for no saving
pars.corrAlignFig = 0; % figure to save CorrAlign data, set to zero for no saving 

% more parameters
% data loading parameters
pars.daxRoot = 'ConvZscan';
pars.maxHyb = inf;  % 0 or inf for all. 
pars.maxFOV = inf;  %  speed up on the fly analysis by restricint FOV 
pars.overwrite = false;
pars.overwriteMaxProject = false;
pars.overwriteAlignment = false;
pars.overwriteCellSegment = false;
pars.verbose = true;
pars.fidChn = inf;  % (inf = auto) should be last, but due to 10 steps and 3 channels 
pars.datChn = inf; % (inf = auto) 

% crop parameters
pars.boxSize = 8; % box-radius of crop area, in pixels (TO DO: should be nanometers for the user)
pars.boxSizeZ = pars.boxSize; % height of fitted area, in z-pixels (To DO: also should be nanometers for the user) 

% % fit parameters  passed to CT4_FitSpots
% % actuall fitting
% pars.nFits = 1; % number of spots to fit per image
% % fine-scale alignment parameters
% pars.maxXYdrift = 8; % in pixels
% pars.maxZdrift = 16; % yes it's measured in z-pixels, but z-drift magnitude is wholly independent of z-drift size
% pars.upsample = 4; % this sets the accuracy of the fine drift correction in xy. It should be converted into nanometers for the user      
% pars.fidRegTheta = 0.6;
% % fine-scale alignment plotting
% pars.saveAlign = false;
% pars.figShowAlign = 10;% ;10; % 0
% pars.figShowFits = 11;% 11;

% initiatializaiton values
currHyb = 1;
prevSpotsCounted = 0;
totSpotsCounted = 0;
refDaxFiles = [];
% % Get hyb folders 
% if ~exist(eTableXLS,'file') % if we have it, may as well use it
%     hybFolders = FindFiles([dataFolder,hybFolderRoot,'*'],'onlyFolders',true);
% else
%     hybFolders = strcat(dataFolder,eTable.FolderName,filesep);
% end

%% background subtractions
toeDax = FindFiles(['J:\20220426_E14_TwistRNA_seDNA\DNA\Toe_002_020\','ConvZscan*dax']);
bkd = cell(length(toeDax),1);
for d=1:length(toeDax)
    bkd{d} = LoadDax(toeDax{d},'maxProject',true,'verbose',false);
end
bkdStk = cat(4,bkd{:});
[fidIms,fidBkd]= FlattenBackground(squeeze(bkdStk(:,:,3,:)),'showPlots',true,'backgroundCorrect','removeData');
[fidIms,fidBkd]= FlattenBackground(squeeze(bkdStk(:,:,3,:)),'showPlots',true,'backgroundCorrect','median');

[dat2Ims,dat2Bkd]= FlattenBackground(squeeze(bkdStk(:,:,2,:)),'showPlots',true,'backgroundCorrect','median');
[dat1Ims,dat1Bkd]= FlattenBackground(squeeze(bkdStk(:,:,1,:)),'showPlots',true,'backgroundCorrect','median');
%%

if pars.maxHyb == 0 || isinf(pars.maxHyb)
    pars.maxHyb = length(hybFolders); 
end

fitFolder = [analysisFolder,'Fits',filesep];
if exist(fitFolder,'dir')==0
    mkdir(fitFolder);
end

% look for existing align table files table files
alignTableFiles = FindFiles([analysisFolder,'alignTable_fov','*.csv']);
if pars.overwriteAlignment
   for f=1:length(alignTableFiles)
       delete(alignTableFiles{f});
   end
end
%=============================== Watcher ================================ %
%
% wait for next hyb
while currHyb+1 < pars.maxHyb
    disp(['waiting for data from Hyb ',num2str(currHyb+1)]); 
    nextHyb = FindFiles([hybFolders{currHyb+1},pars.daxRoot,'*.dax']);
    while isempty(nextHyb) 
        pause(5);
        nextHyb = FindFiles([hybFolders{currHyb+1},pars.daxRoot,'*.dax']);
    end
    currHyb = currHyb + 1;
    %-----------------------------------------------------------------------

    % ========================================================================%
    %                            hyb-to-hyb analysis
    %=========================================================================%
    if pars.verbose
        disp(['now analyzing data from hyb ',num2str(currHyb)]);
    end

    % --------- load refHyb
    if isempty(refDaxFiles)
        refDaxFiles = FindFiles([hybFolders{1},pars.daxRoot,'*.dax']); 
        nFOV = min(pars.maxFOV,length(refDaxFiles));
        refHybFid2D = cell(nFOV,1);
        refHybDat2D = cell(nFOV,1);
        refCrop3D = cell(nFOV,1);
        for f=1:nFOV
            refHybFid2D{f} = LoadDax(refDaxFiles{f},'maxProject',true,'channel',pars.fidChn,'overwrite',pars.overwriteMaxProject);  % this is fast
            refHybDat2D{f} = LoadDax(refDaxFiles{f},'maxProject',true,'channel',2,'overwrite',pars.overwriteMaxProject);
        end
    end


    curDaxFiles = FindFiles([hybFolders{currHyb},filesep,pars.daxRoot,'*.dax']);
    %%  Load the data
    % fid channel is assumed to be the last channel in the stack (typically the
    % most blue shifted, since data is best collected red->blue for
    % photobleaching reasons).  
    % if an alternative channel is used for fiducial, that can be specified as
    % an optional parameter. 

   
    
    for f=1:nFOV % loop over FOVs
        fidAlign =[];
        nSpots = 0;
        
        fovTable = [fitFolder,'spotTable_fov',num2str(f-1,'%03d'),'.csv'];
        if exist(fovTable,'file')
            fovTableData = readtable(fovTable);
            lastHyb = max(fovTableData.hyb);
            if lastHyb >= currHyb && pars.overwrite == false
                spotsInFile = fovTableData.spotID(end)-fovTableData.spotID(1);
                % totSpotsCounted = totSpotsCounted + spotsInFile;
                disp(['fov ',num2str(f),' hyb ',num2str(currHyb),' using existing data with ',num2str(spotsInFile), ' spots.']);
                continue
            end
        end     
        
        
        % read the xml struct
        refHybFid = refHybFid2D{f};
        hybFolderName = hybFolderNames{currHyb}; 
        [curHybIm, imProps] = LoadDax(curDaxFiles{f},'maxProject',false,'verbose',false,'overwrite',pars.overwriteMaxProject);
        if isinf(pars.fidChn)
            fidChn = size(curHybIm,4);
            datChn = 1:fidChn-1;
        elseif isinf(pars.datChn)
            fidChn = pars.fidChn;
            allChns = 1:size(curHybIm,4);
            datChn = allChns;
            datChn(fidChn) = [];
        else
            fidChn = pars.fidChn;
            datChn = pars.datChn;
        end
        curHybFid = max(curHybIm(:,:,:,fidChn),[],3);
        % we only need to load the currHybFid maxProjection for drift
        % correction. However, below we are going to need the full 3D data,
        % and we load it again from scratch, which is slow.  Since it is
        % only one movie, we should load both at full resolution here. 
        
        %%   align drift to hyb 1. 
        % for 'on the fly' hyb 1 is all that exists.  for post-processing, we can
        % have a falg to allow users to choose a different channel
        % 
        
        runDriftCorrect = true;
        if ~isempty(alignTableFiles)
            if length(alignTableFiles)>=f
                if exist(alignTableFiles{f},'file')
                   alignTable = readtable(alignTableFiles{f});
                   fidAlign = table2struct(alignTable(alignTable.hyb==currHyb,:));
                   if ~isempty(fidAlign)
                       runDriftCorrect = false;
                   end
                end
            end
        end
               
       if runDriftCorrect
           if pars.corrAlignFig
            f1 = figure(pars.corrAlignFig); clf;
           end
            fidAlign = CorrAlignFast(refHybFid,curHybFid,'fineUpsample',2,'showplot',pars.corrAlignFig);  

            % save data  
            % save FOV name to table with shift properties 
            % should write an updated version of this
            fidAlign.fov = f;
            fidAlign.hyb = currHyb;
            alignTable = struct2table(fidAlign);
            alignTableFile = [analysisFolder,'alignTable_fov',num2str(f,'%03d'),'.csv'];
            writetable(alignTable,alignTableFile,'WriteMode','Append');
            % 
            % save image too. 
            if pars.corrAlignFig
                SetFigureSavePath([analysisFolder,'CorrAlign',filesep],'makeDir',true);
                imName = ['fov',num2str(f,'%03d'),'_H',num2str(currHyb,'%03d')];
                SaveFigure(f1,'name',imName,'formats',{'png'},'overwrite',true);
            end
        end
        
        % 
        %% segment cells and find ROIs (fiducial spots) 
        % we use refHyb by default. 
        %   this loads saved data by default, so it's okay that we recall it after
        %   every hyb.
        %  used in the SE analysis, we are just getting the cell border
        %  fits here. 
        fovFolder = [analysisFolder,'fov',num2str(f,'%02d'),filesep];
        refHybFlat =  cast(quantile(dat2Bkd(:),.5)*double(refHybDat2D{f})./dat2Bkd,class(refHybDat2D{f}));

        [~,cellID_full,figH] = SegmentSpotsPerNucleus(refHybFlat,'f',f,'saveFolder',fovFolder,'figFormats',{'png'},...
                                        'overwrite',pars.overwrite,'rerunCellpose',pars.overwriteCellSegment,...
                                        'figFormats',{'png'},...
                                        'figName',['hyb',num2str(1,'%03d'),'_fid'],...
                                        'figShowResult',pars.saveSpotsFig);
       cellIDmap = cellID_full; % to crop  % figure(21); clf; imagesc(cellID_full); GetColorMap('distColorsW');
        % troubleshoot
        if currHyb==2
            mask = boundarymask(cellID_full); 
            im = labeloverlay(IncreaseContrast(refHybFlat,'high',.9999),mask,'Transparency',0);
            f4= figure(4); clf; imagesc(im); colormap(gray); colorbar;  hold on; pause(.1);


                SetFigureSavePath(analysisFolder);
                imName = ['cellSegment_fov',num2str(f,'%03d'),'_H',num2str(currHyb,'%03d')];
                SaveFigure(f4,'name',imName,'formats',{'png'},'overwrite',true);
        end



       % cellIDmap = []; % to not crop other cells;
        if currHyb == 2
        % fit spots for hyb 1.  
        %    since the fitter doesn't start until there are 2 hybs
        %    worth of data, we need to catch up and fit hyb 1 now
            [hyb1Im,hy1Props] = LoadDax(refDaxFiles{f},'maxProject',false,'verbose',false,'overwrite',pars.overwriteMaxProject);
            nDataChns = length(datChn); 
            spotTableChn = cell(nDataChns,1);
            for d=datChn % d=2             
               curHybDat = max(hyb1Im(:,:,:,d),[],3);
               % a quick find brightest 2 x,y pixels in cell
               sptPixTable = SegmentSpotsPerNucleus(curHybDat,'f',f,'saveFolder',fovFolder,...
                                                'overwrite',true,'rerunCellpose',false,...
                                                'figFormats',{'png'},...
                                                'figName',['hyb',num2str(1,'%03d'),'_dat',num2str(d)],...
                                                'figShowResult',pars.saveSpotsFig);
               % currHybIm is already loaded and ready to go
               imCrop = CT4_CropSpots(hyb1Im,hy1Props,sptPixTable,...
                    'cellIDmap',cellIDmap,... % cellID_full if you want to crop
                    'boxSize',pars.boxSize);
                % loop over the cropped images 
                nSpots = size(imCrop,1);
                sptTable = cell(nSpots,1);
                parfor t=1:nSpots % t = 5  % can parfor this. 
                    cropIm = imCrop{t,d};
                    sptTable{t} = CT4_SimpleFit(cropIm,hy1Props,sptPixTable(t,:),d,nDataChns,1,'hybName',hybFolderName);  
                end
                spotTableChn{d} = cat(1,sptTable{:});      
            end 
            spotTableCur = cat(1,spotTableChn{:});
            if nSpots > 0
                % totSpotsCounted =  height(spotTableCur);
                % spotTableCur.spotID = (1:totSpotsCounted)';
                fovTable = [fitFolder,'spotTable_fov',num2str(f-1,'%03d'),'.csv'];
                if exist(fovTable,'file')
                    answer = questdlg([fovTable ' exists, overwrite?'], ... % question
                                'Prompt ', ...  % pop-up label
                                'Yes','No','No'); % op1 op2 default
                    if strcmp(answer,'Yes')
                        warning(['deleting existing file ', fovTable]);
                        writetable(spotTableCur,fovTable); 
                    else
                        error('cancelled by user');
                    end
                else
                    writetable(spotTableCur,fovTable);  % start new table
                end  
            end
        end
        
        
        nDataChns = length(datChn); 
        spotTableChn = cell(nDataChns,1);
        for d=datChn % d=1
           curHybDat = max(curHybIm(:,:,:,d),[],3);
           curHybDat = ApplyReg(curHybDat,fidAlign);  % figure(1); clf; imagesc(curHybDat);
    % a quick find brightest 2 x,y pixels in cell
           sptPixTable = SegmentSpotsPerNucleus(curHybDat,'f',f,'saveFolder',fovFolder,...
                                            'overwrite',true,'rerunCellpose',false,...
                                            'figFormats',{'png'},...
                                            'figName',['hyb',num2str(currHyb,'%03d'),'_dat',num2str(d)],...
                                            'figShowResult',pars.saveSpotsFig);
           % currHybIm is already loaded and ready to go
           imCrop = CT4_CropSpots(curHybIm,imProps,sptPixTable,...
                'align',fidAlign,...
                'cellIDmap',cellIDmap,... % cellID_full if you want to crop
                'boxSize',pars.boxSize);
            % loop over the cropped images 
            nSpots = size(imCrop,1);
            sptTable = cell(nSpots,1);
            parfor t=1:nSpots % t = 5  % can parfor this. 
                cropIm = imCrop{t,d};
                sptTable{t} = CT4_SimpleFit(cropIm,imProps,sptPixTable(t,:),d,nDataChns,currHyb,'hybName',hybFolderName);  
            end
            spotTableChn{d} = cat(1,sptTable{:});      
        end
        % figure(2); clf; ProjectIm3D(curHybIm(:,:,:,1)); colormap(gray);
        % figure(3); clf; imagesc(cellIDmap);
        if nSpots>0
            spotTableCur = cat(1,spotTableChn{:});
            % prevSpotsCounted = totSpotsCounted;
            % totSpotsCounted = totSpotsCounted + height(spotTableCur);
            % spotTableCur.spotID = ((prevSpotsCounted+1):totSpotsCounted)'; % unique across all fov, since it just counts up
            fovTable = [fitFolder,'spotTable_fov',num2str(f-1,'%03d'),'.csv'];
            writetable(spotTableCur,fovTable,'WriteMode','Append');  % append to the fov table (this is convienent and compact) 
        end   
    end 
end

elapsed = toc(t_total);
disp(['completed in ', char(duration([0, 0,elapsed]))]);
%%
