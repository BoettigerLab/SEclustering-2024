
%% Overlay RNA and DNA data
derekData = '\\BLabServer1\DerekData01\';
% %% Run the drift correction on the RNA data
% 
rnaFolder = [derekData,'20220426_E14_TwistRNA_seDNA\RNA\'];
rnaTable = readtable([rnaFolder,'RNA_Hyb_Codebook_v2.xlsx'])
rnaHybFolders = FindFiles([rnaFolder,'Hyb*'],'onlyFolders',true);

nFOV =  45;
nH = height(rnaTable);
rna_ims = cell(nH,nFOV,3);
rna_names = cell(nH,2);
for h=1:nH
   % [im,improps] = LoadDax(rnaHybFolders{h},'maxProject',true);
   hybFolder = [rnaFolder,'Hyb_',num2str(rnaTable{h,1},'%03d'),'_',num2str(rnaTable{h,2},'%03d'),'\'];
   daxFiles = FindFiles([hybFolder,'Conv*.dax']);
   for f=1:nFOV
       [im,improps] = LoadDax(daxFiles{f},'maxProject',true,'verbose',false);
       rna_names{h,1} = rnaTable{h,3}{1}; % 750 channel name
       rna_names{h,2} = rnaTable{h,4}{1}; % 647 channel name
       rna_ims{h,f,1} = im(:,:,1); % 750 data
       rna_ims{h,f,2} = im(:,:,2); % 647 data
       rna_ims{h,f,3} = im(:,:,3); % 561 data
   end
end
f=45;
stack_intron = cat(3,rna_ims{:,f,2});
figure(1); clf; Ncolor( .2*IncreaseContrast( stack_intron,'high',.9999,'low',.2));
figure(2); clf; Ncolor( .7*IncreaseContrast( LaplaceFilterImage(  stack_intron,'hsize',5,'sigma',1),'high',.999995,'low',.5));

stack_intron = cat(3,rna_ims{[5,1,3],f,2});
figure(1); clf; Ncolor( .2*IncreaseContrast( stack_intron,'high',.9999,'low',.2));
figure(2); clf; Ncolor( .7*IncreaseContrast( LaplaceFilterImage(  stack_intron,'hsize',5,'sigma',1),'high',.999995,'low',.5));

% 
% 
% stack_exon = cat(3,rna_ims{:,f,1});
% figure(1); clf; Ncolor( .2*IncreaseContrast( stack_exon,'high',.9999,'low',.2));
% figure(3); clf; Ncolor( .4*IncreaseContrast( LaplaceFilterImage(  stack_exon,'hsize',5,'sigma',1),'high',.99995,'low',.5));
% 
% 
% stack_exon = cat(3,rna_ims{[5,4,3,7],f,1});
% figure(1); clf; Ncolor( IncreaseContrast( stack_exon,'high',.99999,'low',.2));
% figure(3); clf; Ncolor( IncreaseContrast( LaplaceFilterImage(  stack_exon,'hsize',5,'sigma',1),'high',.99995,'low',.5));


%%
% analysisFolder = 'J:\20220426_E14_TwistRNA_seDNA\RNA_analysis\'
analysisFolder = [derekData,'20220426_E14_TwistRNA_seDNA\Analysis_CT4_v2\'];
pars.corrAlignFig = 10;
for h=1:nH
    for f=1:nFOV
        alignTableFile = [analysisFolder,'RNAref10_AlignTable_fov',num2str(f,'%03d'),'.csv'];
        if ~exist(alignTableFile,'file')
            currHyb = h;
            refHybFid = rna_ims{nH,f,3};  % Don't use hyb1 as a reference, we fixed the 647 after this hyb.  Using   
            curHybFid = rna_ims{h,f,3};

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

            writetable(alignTable,alignTableFile,'WriteMode','Append');
            % 
            % save image too. 
            if pars.corrAlignFig
                SetFigureSavePath([analysisFolder,'RNAref10_CorrAlign',filesep],'makeDir',true);
                imName = ['fov',num2str(f,'%03d'),'_H',num2str(currHyb,'%03d')];
                SaveFigure(f1,'name',imName,'formats',{'png'},'overwrite',true);
            end
        end
    end
end



%% align RNA to DNA
% we used this as the alignment reference for the DNA analysis
toeDax = FindFiles([derekData,'20220426_E14_TwistRNA_seDNA\DNA\Toe_002_020\','ConvZscan*dax']);
dnaRef = cell(length(toeDax),1);
for d=1:length(toeDax)
    dnaRef{d} = LoadDax(toeDax{d},'maxProject',true,'verbose',false);
end
dnaStk = cat(4,dnaRef{:});
[dnaRef_chn2,dnaBkd_chn2]= FlattenBackground(squeeze(dnaStk(:,:,2,:)),'showPlots',true,'backgroundCorrect','removeData');

rna_ims; % h,f,c
[rnaRef_chn2,rnaBkd_chn2]= FlattenBackground(rna_ims(end,:,2),'showPlots',true,'backgroundCorrect','removeData','strength',.5);

alignTableFile = [analysisFolder,'Align_RNAtoDNA_perFOV.csv'];
if ~exist(alignTableFile)
rnaAlignTable = cell(nFOV,1);
for f=1:nFOV  % f = 38
    rnaAlign2 = CorrAlignFast(dnaRef_chn2(:,:,f),rnaRef_chn2{f},'fineUpsample',1,'showplot',pars.corrAlignFig,'angles',-2:.5:0,'maxShift',100);  
    rnaAlignTable{f} = struct2table(rnaAlign2);
    % save image too. 
    if pars.corrAlignFig
        f1 = figure(pars.corrAlignFig); pause(.1);
        SetFigureSavePath([analysisFolder,'CorrAlign_RNAtoDNA2',filesep],'makeDir',true);
        imName = ['RNAtoDNA_fov',num2str(f,'%03d'),'_H',num2str(size(rna_ims,1),'%03d')];
        SaveFigure(f1,'name',imName,'formats',{'png'},'overwrite',true);
    end
end
rnaAlign_Table = cat(1,rnaAlignTable{:}) % 
writetable(rnaAlign_Table,alignTableFile);  
else
    rnaAlign_Table = readtable(alignTableFile);
end

%% RNA bursts 
% repeat analysis with an 0.2 threshold 
analysisFolder =[derekData, '20220426_E14_TwistRNA_seDNA\Analysis_CT4_v2\'];
chrTableHyb = readtable([analysisFolder,'chrTableHyb.csv']);
clusterThreshold = 0.6; % 0.6  0.3  cluster threshold in um
 maxDistToSE = 2 ;% um  max distance allowed between SE and transcript to be considered valid 

%  % find gene in chrTable
% selGeneID = chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Tbx3'));  %    Pvt1=Myc.  
% selSpotTable = spotTable(spotTable.genomeID==selGeneID,:);
% % find gene in RNA table
% selRnaHyb = StringFind(rnaTable.A1_Notes,'Tbx3');

selGeneIDs = {... % FOr the intron channel, Cy5 ( very not flat)
     chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Nanog')),... % h=1
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Sox2')),... % h=2
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Nanog')),... % h=3
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Tbx3')),... % h=4  % RNA is Tbx5, not on, not an SE, though near Tbx3
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Tbx3')),... % h=5
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Pou5f1')),... % h=6  % Fgf5, also not on, not an SE
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Pvt1')),... % h=7  % Myc
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Prdm14')),... % h=8
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Pou5f1')),... % h=9
         chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Nanog')),... % h=10
}
geneNames = {'Nanog','Sox2','Nanog','Tbx5','Tbx3','Fgf5','Pvt1','Prdm14','Pou5f1','Nanog'};%  Pou5f1 = Oct4

% Exon data (not used)
% selGeneIDs = {... % FOr the exon channel, 750, flattens nicely
%      chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Sox2')),... % h=1
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Pou5f1')),... % h=2
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Nanog')),... % h=3
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Tbx3')),... % h=4  %
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Tbx3')),... % h=5
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Pou5f1')),... % h=6  % 
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Pvt1')),... % h=7  % 
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Prdm14')),... % h=8
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Pou5f1')),... % h=9
%          chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Nanog')),... % h=10
% }

% nanogID = chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,'Nanog'));
% nanogSpots = spotTable(spotTable.genomeID==nanogID,:);

plotIms = false;
datTables = FindFiles([analysisFolder,'Fits\fov*masterTable.csv']);
nFOV = length(datTables);
rnaBurst_dnaClust = cell(nH,nFOV,2);
rna_ims; % h,f,c

mapsON = cell(10,1); % {h} = nan(84,84,1000); ccON = 0;
mapsOFF = cell(10,1); % {h} = nan(84,84,1000); ccOFF = 0;
for h = 1:10 % h=nH  % h=5
mapsON{h} = nan(88,88,5000); ccON = 0;
mapsOFF{h} = nan(88,88,5000); ccOFF = 0;

[rna_chn1,rnaBkd_chn1]= FlattenBackground(rna_ims(h,:,1),'showPlots',plotIms,'backgroundCorrect','median','strength',.5);  % exon channel, 750, low SNR, gen flat
[rna_chn2,rnaBkd_chn2]= FlattenBackground(rna_ims(h,:,2),'showPlots',plotIms); % cy5, better SNR, badly off center in h1, then fixed.
[rna_chn3,rnaBkd_chn3]= FlattenBackground(rna_ims(h,:,3),'showPlots',plotIms,'backgroundCorrect','median','strength',.5);
    disp(geneNames{h})
    for f=1:nFOV
        
        % move to DNA coordinates
           % first move this h to match the RNA ref, nH
           % then
           % f = 5;  % 8 5 is a nice cell / segment
           rnaDriftTableFile = [analysisFolder,'RNAref10_AlignTable_fov',num2str(f,'%03d'),'.csv'];
           rnaDriftTable = readtable(rnaDriftTableFile);
          %  rna_im = rna_chn1{f};
           rna_im = rna_chn2{f};
           rna_im = ApplyReg(rna_im,rnaDriftTable(h,:)); %    correct RNA drift to RNA ref hyb (last hyb)  
           rna_im = ApplyReg(rna_im,rnaAlign_Table(f,:),'invert',false);   %  convert RNA to DNA coordinates
           
           % % exon data  (not used)
           % rna_im = rna_chn1{f};
           % rna_im = ApplyReg(rna_im,rnaDriftTable(h,:)); %    correct RNA drift to RNA ref hyb (last hyb)  
           % rna_im = ApplyReg(rna_im,rnaAlign_Table(f,:),'invert',false);   %  convert RNA to DNA coordinates 
           % 
            cellIDs = imread([analysisFolder,'fov',num2str(f,'%02d'),'/cellpose/image_001_cp_masks.tif']); 
            if plotIms
                figure(1); clf; imagesc(IncreaseContrast(rna_im,'high',.9999,'low',.4));
                
                cellID_full = imresize(cellIDs,size(rna_im));
                mask = boundarymask(cellID_full); 
                im = labeloverlay(IncreaseContrast(rna_im,'high',.9999,'low',.4),mask,'Transparency',0);
                f4= figure(4); clf; imagesc(im); colormap(gray); colorbar;  hold on; pause(.1);
    
                buffer = 40;
                figure(1); clf;  imagesc(IncreaseContrast(rna_im(1+buffer:end-buffer,1+buffer:end-buffer),'high',.9999','low',.4)); hold on;
                xy = FindSpots(rna_im(1+buffer:end-buffer,1+buffer:end-buffer),'showPlot',true,'autoSelectThreshold',.995);
                figure(4); plot(xy(:,1)+buffer,xy(:,2)+buffer,'yo')
            end

            % find RNA peaks
            sptPixTable = SegmentSpotsPerNucleus(rna_im,'f',f,...
                                            'cellID',cellIDs,...
                                            'figShowResult',0);

            % this then took the brightest pixel in each 6-pixel box, as
            % the Gaussian fits sometimes lie a bit off peak. 
            imPropCrop = improps; 
            imPropCrop.nChns = 1;
            imPropCrop.zSteps = 1;
            imCrop = CT4_CropSpots(rna_im,imPropCrop,sptPixTable,...
                        'boxSize',6,'chns',1);
            nTraces = height(sptPixTable);
            maxBright = nan(nTraces,1);
            for t=1:length(imCrop)
                maxBright(t) = max(imCrop{t}(:));
            end
            sptPixTable.maxBright = maxBright;
            % figure(1); hold on; plot(sptPixTable.x_pix,sptPixTable.y_pix,'m+')
           %  figure(5); clf; hist(maxBright,nTraces);

           % load DNA
            spotTable = readtable([analysisFolder,'Fits\fov',num2str(f,'%03d'),'_masterTable.csv']);
            [cellNums,idx] = unique(spotTable.cellID);
           
%             % just plotting (DNA SE positions and cell numbers)
            if plotIms
                figure(1); hold on; 
                hold on; plot(spotTable.x/improps.xy2um,spotTable.y/improps.xy2um,'r+')
                text(spotTable.x_fid(idx)/improps.xy2um,spotTable.y_fid(idx)/improps.xy2um,cellstr(num2str(cellNums)),'Color',[0 1 1]);
                
                idx = ismember(spotTable.genomeID,selGeneIDs{h}) ; % mulitple SEs for the same gene
                geneSpots = spotTable(idx,:);
                hold on; plot(geneSpots.x/improps.xy2um,geneSpots.y/improps.xy2um,'g+')
            end
            nCells = max(cellNums);
            clusterSize = nan(nCells,2); % up to 2 alleles per cell 
            nascentRNA = nan(nCells,2); 
            xyz_all_bursts = [sptPixTable.x_pix*improps.xy2um, sptPixTable.y_pix*improps.xy2um,5*ones(height(sptPixTable),1),...
                              sptPixTable.maxBright];
            for c=1:nCells % c = 32
                cellTable = spotTable(spotTable.cellID==c,:); % SE DNA from this cell   
                idx = ismember(cellTable.genomeID,selGeneIDs{h}) ; % mulitple SEs for the same gene
                geneSpots = cellTable(idx,:); % SEs for the gene being tested   
                xyz_burst = xyz_all_bursts(sptPixTable.cellID==c,:);  % RNA for the gene 
                nBursts = size(xyz_burst,1);
                xyz_gene = [geneSpots.x,geneSpots.y];
                if ~isempty(geneSpots)
                    [idx,d] = knnsearch(xyz_gene, xyz_burst(:,1:2));
                    for b=1:nBursts
                        if d(b) > maxDistToSE % if it is a stray spot
                            xyz_burst(b,4) = 0; % no expression if brightest spot is a stray spot
                            xyz_burst(b,1:2) = xyz_gene(idx(b),:);  % use SE as a proxy for the gene if we don't have a nascent transcript to use. 
                        end
                    end
                end
                xyz = [cellTable.x, cellTable.y, cellTable.z]; %  SE locs in um
                isSelGene = ismember(cellTable.genomeID,selGeneIDs{h}); % sometimes multiple per 
                nPts = 88;% (# SEs in diploid data)
                geneIDs = find(isSelGene);
                if ~isempty(geneIDs)
                    for b=1:nBursts
                        nascentRNA(c,b) = xyz_burst(b,4);  % end up not using the actual burst xyz coordinates, just record its brightness, we'll use the SE coordinates for uniformity between the ONs and OFFs  
                         dmap = nan(nPts,nPts);
                        hasData = ~isnan(xyz(:,1));
                        bc = cellTable.genomeID + 44*cellTable.allele;
                        bc2 = bc(hasData);
                        if sum(hasData)>1
                            dmap(bc2,bc2) = squareform(pdist(xyz(hasData,1:3))); % 
                            % figure(3); clf; imagesc(dmap);
                            dmap(dmap==0) = nan; % 
                            clusterSize(c,b) = nansum(dmap(geneIDs(idx(b)),:) < clusterThreshold,2);   
                        end

                        if xyz_burst(b,4) > 1000
                            ccON = ccON+1;
                            mapsON{h}(:,:,ccON) = dmap;
                        else
                            ccOFF = ccOFF+1;
                            mapsOFF{h}(:,:,ccOFF) = dmap;
                        end
                    end
                end
            end
              rnaBurst_dnaClust{h,f,1} = nascentRNA(:);
              rnaBurst_dnaClust{h,f,2} = clusterSize(:);
            if plotIms
             figure(5); clf; PlotCorr(clusterSize(:),nascentRNA(:),'log',false);
            end
       
    end
end
%%
rnaBurst_dnaCluster = rnaBurst_dnaClust(dispGenes,:,:);

saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\'
% save([saveFolder,'RnaBurst_DnaCluster.mat'],'rnaBurst_dnaCluster');


% save([saveFolder,'RnaBurst_Maps.mat'],'mapsON','mapsOFF','-v7.3');
