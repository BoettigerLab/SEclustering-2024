
%%  This step computes at 3D chromatic correction and saves the corresponding tform
% this chromatic correction will be loaded in the next step, and applied to
% the data before assigning alleles.
 

%% 
% this time all datasets are in the same place
% v15 data


expFolder = 'J:\Derek_TEMP\20210722_L10_SE_MultiplexedCells\';
dataFolder = [expFolder,'DNA_Expt\'];
analysisFolder = [expFolder,'Analysis_CT4_v15\'];
fitFolder = [analysisFolder,'Fits\'];

%%
% load Derek's reference table
refTable = readtable('J:\Derek_TEMP\20210722_L10_SE_MultiplexedCells\20210722 L10 multiplexed cells Experiment_2021-08-04.xlsx');
hybTable = refTable;

%% match CT4 barcode IDs to actual barcodes

a1Table = table(hybTable.A1_Cy5,hybTable.A1_chr,hybTable.A1_L10_mm9_start,hybTable.A1_L10_mm9_end,...
    hybTable.A1_Young_ID,hybTable.A1_mRNA_ID,hybTable.A1_Young_gene,647*ones(height(hybTable),1),(2:2:2*height(hybTable))');
a1Table.Properties.VariableNames = {'read','chr','start','stop','Yo_ID','mRNA_ID','Yo_gene','chn','barID'};

a2Table = table(hybTable.A2_750,hybTable.A2_chr,hybTable.A2_L10_mm9_start,hybTable.A2_L10_mm9_end,...
    hybTable.A2_Young_ID,hybTable.A2_mRNA_ID,hybTable.A2_Young_gene,750*ones(height(hybTable),1),(1:2:2*height(hybTable))');
a2Table.Properties.VariableNames = {'read','chr','start','stop','Yo_ID','mRNA_ID','Yo_gene','chn','barID'};

% it's useful to know what gene was co-stained, and who was what color, so
%   we'll keep columns 1:6 along for the ride
a1Table = cat(2,a1Table,hybTable(:,1:6)); 
a2Table = cat(2,a2Table,hybTable(:,1:6)); 
readTable = cat(1,a2Table,a1Table); % now we stack these tables

%% finalize the matched tables. 
isA = strcmp(a1Table.DataType,'A'); %  & nonhTable.chn==647;
repTable = a1Table(isA,:);  % this has the chromatic swap datasets
matchHybs = unique(regexprep(repTable.FolderName,'Chrom','Hyb'));
isMatch = StringFind(a1Table.FolderName, matchHybs,'exactly',true);
matchTable = a1Table(isMatch,:);  % this table just has the matched datasets
% now we have two tables, which here at least checkout as sorted in the
% same order, where one is the chormatic flips of the other. 
%%  Test an example
nB = height(matchTable);
fitFiles = FindFiles([fitFolder,'spotTable_fov*.csv']);
nFOV = length(fitFiles);
allPts = cell(nFOV,nB);
% load the spot data
for f = 1:nFOV % select an FOV % f=1
    fovTable = readtable(fitFiles{f});
    cntr = StringFind(fovTable.Properties.VariableNames,{'x_fid','y_fid','z_fid'},'exactly',true);
    fits = StringFind(fovTable.Properties.VariableNames,{'x_um','y_um','z_um'},'exactly',true);
    
    for b=1:nB % select one of the 4 chormatic flips to work with first
        isH750 = strcmp(fovTable.hybName,matchTable.FolderName{b}) & fovTable.dat_chn==1;
        isR647 = strcmp(fovTable.hybName,  repTable.FolderName{b}) & fovTable.dat_chn==2;
        xyz750 = fovTable{isH750,cntr}+ fovTable{isH750,fits};
        xyz647 = fovTable{isR647,cntr}+ fovTable{isR647,fits};

        % (just test)
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
        figure(2); clf;
        plot(pts750(:,1),pts750(:,2),'rs'); hold on;
        plot(pts647(:,1),pts647(:,2),'co');

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

        % correct extrac drift
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
tform3D = Polymap3D(all647,all750,'max2D',.300,'max3D',.400,'polyOrder',2,'bins',0:.01:.5,'units','um');
%% save

save([analysisFolder,'tform3D_750to647.mat'],'tform3D');
disp(['wrote ',analysisFolder,'tform3D_750to647.mat']);


