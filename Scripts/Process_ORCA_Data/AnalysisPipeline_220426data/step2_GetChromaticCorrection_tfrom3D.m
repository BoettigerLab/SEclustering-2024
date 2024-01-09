
%% 
% this time all datasets are in the same place
% v15 data

expFolder = 'J:\20220426_E14_TwistRNA_seDNA\';
dnaFolder = [expFolder,'DNA\'];
analysisFolder = [expFolder,'Analysis_CT4_v2\'];
fitFolder = [analysisFolder,'Fits\'];

%%

hybTable = readtable([dataFolder,'2022-04-26_E14-mES_TwistRNA_seDNA_eTable_v2.xlsx']);
hybFolders = strcat(dataFolder,eTable.FolderName,filesep);
hybFolderNames = hybTable.FolderName;

% the Reps, As and toes are at the beginning in this analysis
isH = strcmp(hybTable.DataType,'H');
hybTable = cat(1,hybTable(~isH,:),hybTable(isH,:))
hybFolders = cat(1,hybFolders(~isH),hybFolders(isH))
hybFolderNames = cat(1,hybFolderNames(~isH),hybFolderNames(isH))


%% match CT4 barcode IDs to actual barcodes

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

%% finalize the matched tables. 
isA = strcmp(a1Table.DataType,'A'); %  & nonhTable.chn==647;
repTable = a1Table(isA,:);  % this has the chromatic swap datasets
repHyb = repTable.FolderName{1};

nameParts = strsplit(repHyb,'_');
matchHyb = ['Hyb_',nameParts{3},'_',nameParts{2}]
% matchHybs = unique(regexprep(repTable.FolderName,'Chrom','Hyb'));

isMatch = StringFind(a1Table.FolderName, matchHyb,'exactly',true);
matchTable = a1Table(isMatch,:);  % this table just has the matched datasets
% now we have two tables, which here at least checkout as sorted in the
% same order, where one is the chormatic flips of the other. 


% Note, xy drift correction has already been applied. 

%%  Compute alignment
nB = height(matchTable);
fitFiles = FindFiles([fitFolder,'spotTable_fov*.csv']);
nFOV = length(fitFiles);
allPts = cell(nFOV,nB);
% load the spot data
for f = 1:nFOV % select an FOV % f=1
    fovTable = readtable(fitFiles{f});
    cntr = StringFind(fovTable.Properties.VariableNames,{'x_fid','y_fid','z_fid'},'exactly',true);
    fits = StringFind(fovTable.Properties.VariableNames,{'x_um','y_um','z_um'},'exactly',true);
    
    for b=1:nB % select one of the chormatic flips 

        % we do the two flipped genes separately, to avoid confusion in the
        % matching.  We start with match-table gene chn 1:
        isH750 = strcmp(fovTable.hybName,matchTable.FolderName{b}) & fovTable.dat_chn==1;  % in FOV, find the matchTable hyb, chn1
        isR647 = strcmp(fovTable.hybName,  repTable.FolderName{b}) & fovTable.dat_chn==2;  % find the rep table Hyb chn 2 (this is actually hte same gene)
        xyz750 = fovTable{isH750,cntr}+ fovTable{isH750,fits};
        xyz647 = fovTable{isR647,cntr}+ fovTable{isR647,fits};
    
        % (just test) -- these should largely overlap, modulo linear drift
        figure(1); clf;
        plot(xyz750(:,1),xyz750(:,2),'rs'); hold on;
        plot(xyz647(:,1),xyz647(:,2),'co');

        % simple nearest neighboring matching 
        [idx,dis] = knnsearch(xyz647,xyz750);
        hasMatch = dis<.5; % max 300nm separation
        pts647 = xyz647(idx,:);
        pts647 = pts647(hasMatch,:);
        pts750 = xyz750(hasMatch,:);

        % just test
        figure(1); clf;
        plot(pts750(:,1),pts750(:,2),'rs'); hold on;
        plot(pts647(:,1),pts647(:,2),'co');

        % Now we get match-table gene chn2 
        isR750 = strcmp(fovTable.hybName,  repTable.FolderName{b}) & fovTable.dat_chn==1;
        isH647 = strcmp(fovTable.hybName,matchTable.FolderName{b}) & fovTable.dat_chn==2;
        xyz750b = fovTable{isR750,cntr}+ fovTable{isR750,fits};
        xyz647b = fovTable{isH647,cntr}+ fovTable{isH647,fits};
        % simple nearest neighboring matching 
        [idx,dis] = knnsearch(xyz647b,xyz750b);
        hasMatch = dis<.5; % max 300nm separation
        pts647b = xyz647b(idx,:);
        pts647b = pts647b(hasMatch,:);
        pts750b = xyz750b(hasMatch,:);
        figure(2); clf;
        plot(pts750b(:,1),pts750b(:,2),'rs'); hold on;
        plot(pts647b(:,1),pts647b(:,2),'co');

        % these two should agree on chromatic correction but have opposite ideas of
        % the residual drift.  

        % correct extra drift
          % this does help, an extra 35 nm 3D, 40 nm xy  (2x better
          % resolution - 0.81 um vs 0.38 um final error)
          %  
        drift = mean(...
        [nanmean(pts647 - pts750)
        -nanmean(pts647b -pts750b)]);

        
        pts750A = pts750 + repmat(drift,size(pts750,1),1);
        pts647A = pts647;
        tform3Da = Polymap3D(pts647A,pts750A,'max2D',.300,'max3D',.400,'polyOrder',2,'bins',0:.01:.5,'units','um');

        pts750B = pts750b - repmat(drift,size(pts750b,1),1);
        pts647B = pts647b;
        tform3Db = Polymap3D(pts647B,pts750B,'max2D',.300,'max3D',.400,'polyOrder',2,'bins',0:.01:.5,'units','um');

        allPts{f,b,1} = cat(1,pts750A,pts750B);
        allPts{f,b,2} = cat(1,pts647A,pts647B);

        % new750 = tforminv(tform3Da,pts750B(:,1),pts750B(:,2),pts750B(:,3));
        new750 = tforminv(tform3Da,pts750B);
        disp([nanmedian(sqrt(sum((pts750B - pts647B).^2,2)))
        nanmedian(sqrt(sum((new750 - pts647B).^2,2)))])

        new750 = tforminv(tform3Db,pts750B(:,1),pts750B(:,2),pts750B(:,3));
        disp([nanmedian(sqrt(sum((pts750B - pts647B).^2,2)))
        nanmedian(sqrt(sum((new750 - pts647B).^2,2)))])
    end
end

%%
all647 = cat(1,allPts{:,:,1});
all750 = cat(1,allPts{:,:,2});
tform3D = Polymap3D(all647,all750,'max2D',.300,'max3D',.400,'polyOrder',2,'bins',0:.01:.5,'units','um','zscale',.5);

%% iterative correction
% only consider points as pairs that got reasonably close / improved on
% dirft correction.  This should filter out the erroneous pairs
xyz705Fix = tforminv(tform3D,all750(:,1),all750(:,2),all750(:,3));

% simple nearest neighboring matching 
[idx,dis] = knnsearch(all647(:,1:2),xyz705Fix(:,1:2));
hasMatch = dis<.04;
pts647 = all647(idx,:);
pts647 = pts647(hasMatch,:);
pts750 = all750(hasMatch,:);

tform3D_iter2 = Polymap3D(pts647,pts750,'max2D',.300,'max3D',.400,'polyOrder',3,'bins',0:.01:.5,'units','um','zscale',.5);

% subplot(2,2,1); cla; subplot(2,2,2); cla;
%% save
% 
% save([analysisFolder,'tform3D_750to647.mat'],'tform3D');
% disp(['wrote ',analysisFolder,'tform3D_750to647.mat']);
% 
% save([analysisFolder,'tform3D_iter2_750to647.mat'],'tform3D_iter2');
% disp(['wrote ',analysisFolder,'tform3D_iter2_750to647.mat']);




%% correct 
analysisFolder
for f=1:nFOV
    fovTable = readtable(fitFiles{f});
    cntr = StringFind(fovTable.Properties.VariableNames,{'x_fid','y_fid','z_fid'},'exactly',true);
    fits = StringFind(fovTable.Properties.VariableNames,{'x_um','y_um','z_um'},'exactly',true);
    
    isGood_isChn1 = fovTable.dat_chn==1 & ~isnan(fovTable.x_um);

    xyz750 = fovTable{isGood_isChn1,cntr}+ fovTable{isGood_isChn1,fits};
    xyz_fix = tforminv(tform3D_iter2,xyz750(:,1),xyz750(:,2),xyz750(:,3));
    fovTableCorrect = fovTable;
    fovTableCorrect{isGood_isChn1,fits} = xyz_fix;
    writetable(fovTableCorrect,[analysisFolder,'spotTableChromCorrect_fov',num2str(f,'%03d'),'.csv']);
end
