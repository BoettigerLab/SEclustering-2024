%% Create Master Tables
% this step takes the fitted spot-tables, applies the chromatic correction
% computed in the previous step, and uses a brute-force clustering approach
% to assign alleles as well.  
% The resulting data-tables are complete with Cell_ID, allele_ID
% and other variables all embedded, and chromatic aberrations already
% corrected.  We call thest the master tables.
% These data still need to be demultiplexed by the barcodes, that happens
% in the final step before general analysis. 
%
% input "spotTable"
% ouptut "masterTable"
% Note - not actually the final table, we have yet to add the group labels
% from the demultiplexing, when we do this in the next step we will
% actually have the final tables ready for exploring. 
% ============================================================% 

%% Notes from Derek on the 220426 experiment
% From my notes, it looks like Rep_001_019 might be slightly more relevant than Hyb_001_019 to the rest of the images because I accidentally bleached the sample for ~10 frames with the wrong Pars file, which would mostly affect only the fiducial channel? Otherwise, the rest are control replicates


%% step 0 load data table and ref tables


expFolder = 'J:\20220426_E14_TwistRNA_seDNA\';
dnaFolder = [expFolder,'DNA\'];
analysisFolder = [expFolder,'Analysis_CT4_v2\'];
fitFolder = [analysisFolder,'Fits\'];

% ------- Load and resort hybTable as in original dax analysis
hybTable = readtable([dataFolder,'2022-04-26_E14-mES_TwistRNA_seDNA_eTable_v2.xlsx']);
hybFolders = strcat(dataFolder,eTable.FolderName,filesep);
hybFolderNames = hybTable.FolderName;

% the Reps, As and toes are at the beginning in this analysis
isH = strcmp(hybTable.DataType,'H');
hybTable = cat(1,hybTable(~isH,:),hybTable(isH,:))
hybFolders = cat(1,hybFolders(~isH),hybFolders(isH))
hybFolderNames = cat(1,hybFolderNames(~isH),hybFolderNames(isH))


% ----- match CT4 barcode IDs to actual barcodes
a1Table = table(hybTable.A1_Cy5,hybTable.A1_chr,hybTable.A1_L10_mm9_start,hybTable.A1_L10_mm9_end,...
    hybTable.A1_Cy5_gene,647*ones(height(hybTable),1),(2:2:2*height(hybTable))');
a1Table.Properties.VariableNames = {'read','chr','start','stop','Yo_gene','chn','barID'};

a2Table = table(hybTable.A2_750,hybTable.A2_chr,hybTable.A2_L10_mm9_start,hybTable.A2_L10_mm9_end,...
    hybTable.A2_750_gene,750*ones(height(hybTable),1),(1:2:2*height(hybTable))');
a2Table.Properties.VariableNames = {'read','chr','start','stop','Yo_gene','chn','barID'};

% it's useful to know what gene was co-stained, and who was what color, so
%   we'll keep columns 1:6 along for the ride
a1Table = cat(2,a1Table,hybTable(:,1:6)); 
a2Table = cat(2,a2Table,hybTable(:,1:6)); 
readTable = cat(1,a2Table,a1Table); % now we stack these tables

%%
% now we sort the readouts by chromosome, since we will cluster them by
% chromosome.  
cnt =0;
chrs = unique(readTable.chr);
chrTables = cell(length(chrs),1);
for c=1:length(chrs) % c=11
    isC = strcmp(readTable.chr,chrs{c});
    chrNum = str2double(regexprep(chrs{c},'chr',''));
    chrTable = readTable(isC,:);
    [~,i] = sort(chrTable.start);
    chrUID = ((cnt+1):cnt+sum(isC))';
    cnt = cnt+sum(isC);
    if isnan(chrNum)% handle X
        chrNum = 23; 
    end
    chrNum = repmat(chrNum,length(chrUID),1);
    chrTables{c} = cat(2,chrTable(i,:),table(chrUID,chrNum)); % apply sort
end
chrTable = cat(1,chrTables{:}); % combine again; 
[~,si] = sort(chrTable.chrNum);
chrTable = chrTable(si,:);
chrs = unique(chrTable.chr,'stable');
%% step 1.2 chromatic correction

load([analysisFolder,'tform3D_iter2_750to647.mat'],'tform3D_iter2');
disp(['loaded ',analysisFolder,'tform3D_iter2_750to647.mat']);




%% Build Polymer Lists and Distance Maps, separating chromosomes
% (this takes a little time, mostly spent on a small number of cells with
% >20 points per chromosome.  Could be accelerated). 
fitFiles = FindFiles([fitFolder,'spotTable_fov*.csv']);
nFOVs = length(fitFiles);

pars.figChrTer = 0; % optionall show figures
pars.figChrTer3D = 0;
spacer = 10;  % gap to leave in the map between 'mat' and 'pat' chrs.
maps = cell(nFOVs,1);
polys = cell(nFOVs,1);

dims = 1:3;
nDims = length(dims);


isHyb = contains(chrTable.FolderName,'Hyb');
idHyb = chrTable(isHyb,:).barID;
% repHyb = {'Hyb_035','Hyb_036','Hyb_037','Hyb_063','Hyb_066','Hyb_067'};
% isBad = contains(chrTable.FolderName,repHyb);
% idBad = chrTable(isBad,:).barID;
% useRep = {'Rep_035','Rep_036','Rep_037','Rep_063','Rep_066','Rep_067'};
% isRep = contains(chrTable.FolderName,useRep);
% idRep = chrTable(isRep,:).barID;
% [idBad,idRep]
% remove bad from ChrTable
chrTableHyb = chrTable;
% chrTableHyb(isBad,:) = chrTableHyb(isRep,:);
chrTableHyb(~isHyb,:) = [];
chrTableHyb.chrUID = (1:height(chrTableHyb))';

% writetable(chrTableHyb,[analysisFolder,'chrTableHyb.csv']);

nS = height(chrTableHyb);
nB = 2*nS+2*spacer;

% 

chrColorMap = GetColorMap('distColors',length(chrs));
masterTables = cell(nFOVs,1);
for f=1:nFOVs % f =9
    fovTableCSV = [fitFolder,'spotTable_fov',num2str(f-1,'%03d'),'.csv'];
    try
        % load data
        fovTabAll = readtable(fovTableCSV);
        isHybData = ismember(fovTabAll.barcodeID,idHyb);
        fovTab = fovTabAll(isHybData,:);

%         % Replace failed hybes with designated repeats 
%         isRep = ismember(fovTabAll.barcodeID,idRep);
%         isBad = ismember(fovTab.barcodeID,idBad);
%         fovTab(isBad,:) = []; % remove bad hybs
%         fovTab = cat(1,fovTab,fovTabAll(isRep,:)); % replace with replicate hybs 
        
        % convert to FOV coordinates
        fovTab.x  = fovTab.x_um + fovTab.x_fid; % a little easier to have these combined.
        fovTab.y  = fovTab.y_um + fovTab.y_fid;
        fovTab.z  = fovTab.z_um + fovTab.z_fid;
        % correct chromatic aberration (relative to fov)
        is750 = fovTab.dat_chn==1 & ~isnan(fovTab.x);
        new750 = tforminv(tform3D_iter2,fovTab.x(is750),fovTab.y(is750),fovTab.z(is750));
        fovTab.x(is750) = new750(:,1);
        fovTab.y(is750) = new750(:,2);
        fovTab.z(is750) = new750(:,3);
        
    catch
        continue
    end
    nCells = max(fovTab.cellID);
    disp(['processing fov ',num2str(f),' of ',num2str(nFOVs) ', cells=',num2str(nCells)]);
    mapF = cell(nCells,1); %  nan(nB,nB,nCells);
    polyF = cell(nCells,1); % nan(nB,size(fovTab,2)+6,nCells);
    newTables = cell(nCells,1);
    ct_barID = chrTableHyb.barID;
    ct_chr = chrTableHyb.chr;
    ct_chrUID = chrTableHyb.chrUID;
    parfor c=1:nCells % c = 10;  % should be able to parfor over cells
        disp(['cell = ',num2str(c),' of ',num2str(nCells)]);
        polyF{c} = nan(nB,size(fovTab,2)-2+3); % remove 2 non-x tables
        mapF{c} = nan(nB,nB);
        try
        %  disp(['processing fov ',num2str(f),' of ',num2str(nFOVs) ', cell ',num2str(c),' of ',num2str(nCells)]);
        cellTab = fovTab(fovTab.cellID==c ,:);
%         cellTab.x  = cellTab.x_um + cellTab.x_fid; % a little easier to have these combined.
%         cellTab.y  = cellTab.y_um + cellTab.y_fid;
%         cellTab.z  = cellTab.z_um + cellTab.z_fid;
        cellTab.genomeID = zeros(height(cellTab),1); % unique index
        cellTab.allele = zeros(height(cellTab),1);
        cellTab.chrID = zeros(height(cellTab),1);
        cellTabBarcodeID = cellTab.barcodeID; 
        for ch = 1:length(chrs)  % with some work we could par for this
             barIDs = ct_barID(strcmp(ct_chr,chrs{ch}));              %#ok<PFBNS>
             isChr = ismember(cellTabBarcodeID,barIDs) & ~isnan(cellTab.x);
             cellTab.chrID(isChr) = ch;
             [~,ids] = ismember(cellTab.barcodeID,ct_barID);
             cellTab.genomeID = ct_chrUID(ids);   %#ok<PFBNS>
             rawPts = [cellTab.x(isChr),cellTab.y(isChr),cellTab.z(isChr),cellTab.barcodeID(isChr)];
             if pars.figChrTer
                figure(pars.figChrTer); clf;
                plot(cellTab.x,cellTab.y,'k.','MarkerSize',1); hold on;
             end
             sortPts = AssignAllelesFromCluster(rawPts,...
                 'figScatter',0,'figMap',0,'verbose',false,...
                 'maxMatch',14);
             cellTab.allele(isChr) = sortPts(:,6)-1; % convert back to 1/0
        end
          
          %----------------
%           if pars.figChrTer3D
%                figure(2); clf;
%                alleleSym = {'o','>'};
%                alleleSz = [1.5,4];
%                for ch = 1:length(chrs)  % 3
%                    for a=1:2 
%                         i = cellTab.chrID==ch & cellTab.allele==a-1;
%                         xyz = 20*[cellTab.x(i),cellTab.y(i),cellTab.z(i)];
% %                         [~,o] = sort(cellTab.genomeID(i));
% %                         xyz = xyz(o,:);
% %                         PlotPolymerTube(xyz,'center',false,'sphereRadius',2,'tubeRadius',1,...
% %                             'lightOn',false,'colormap',chrColorMap(ch,:),'method','spline');
%                          PlotSpheres(xyz,'r',alleleSz(a),'color',chrColorMap(ch,:),'alpha',1,'lightingOn',false);
%                         hold on;
%                        plot3(xyz(:,1),xyz(:,2),xyz(:,3),'color',chrColorMap(ch,:)); hold on;
%                    end
%                end
%                material dull;
%                camlight left;
%                lighting gouraud;
%                set(gca,'color','w');
%                % should save figure for browsing, along with ID numbers
%           end          
          %---------------
         
          barcodeIndex = cellTab.genomeID + cellTab.allele*(nS+spacer) ; % 
          barcodeIndex(barcodeIndex==0) = 2*nS+2*spacer;
          nBs = length(barcodeIndex);
          remCols = StringFind(cellTab.Properties.VariableNames,{'spotID','hybName'});
          cellTab(:,remCols) = [];
          polyF{c}(barcodeIndex,:) = cellTab{:,:};
          xyz = [cellTab.x,cellTab.y,cellTab.z];
          mapF{c}(barcodeIndex,barcodeIndex) = squareform(pdist(xyz));     
          % figure(2); clf; imagesc(maps{f}(:,:,c));
          newTables{c} = cellTab;
        catch er 
            disp(['cell ',num2str(c)]);
            warning(er.message);
        end
    end
    polys{f} = cat(3,polyF{:});
    maps{f} = cat(3,mapF{:});
    
    masterTables{f} = cat(1,newTables{:});
    if ~isempty(masterTables{f})
    writetable(masterTables{f},[fitFolder,'fov',num2str(f,'%03d'),'_masterTable.csv']); 
    else
       disp(['no data in FOV',num2str(f)]) ;
    end
end
  
%%

allMaps = cat(3,maps{:});
allPols = cat(3,polys{:});

% show the full map, including the inter-chromosomal maps
mMap = nanmedian(allMaps,3);
figure(2); clf; imagesc(mMap); colorbar;

% combine both alleles of each chromosome
r2 =(nS+spacer+1):(2*nS+spacer); r1 = 1:nS;
catMaps = cat(3, allMaps(r1,r1,:), allMaps(r2,r2,:));
catPolys = cat(3, allPols(r1,:,:), allPols(r2,:,:));

% compute and show combined chr maps
[cMap1,n1] = ContactFrac(catMaps,'threshold',1.5); % in microns now
totCells = size(catMaps,3); % actually tot chromosomes
dMap = nanmedian(catMaps,3);
figure(3); clf; 
subplot(1,3,1); imagesc(log2(cMap1)); colorbar;
subplot(1,3,2); imagesc(dMap); colorbar;
subplot(1,3,3); imagesc(n1); colorbar;
set(gcf,'color','w');  % axis image;

size(allPols)

%%
dMap2 = dMap;
dMap2([22,43],:) = [];
dMap2(:,[22,43]) = [];
figure(10); clf; imagesc(dMap2); colorbar;  axis image;


