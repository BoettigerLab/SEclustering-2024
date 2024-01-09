
% Huang Ren 2021 data downloaded from 4DN portal
traceT = readtable('U:\GenomeData\ByPublication\HuangRen2021\4DNFI2ZSFLHH.csv','NumHeaderLines',16);
spotT = readtable('U:\GenomeData\ByPublication\HuangRen2021\4DNFI9KE6AII.csv');

% Paired on data portal as biorep 1 tech rep 1,   a name "4DNESNEKOCAP" (which also is not preserved anywhere) 
traceT2 = readtable('U:\GenomeData\ByPublication\HuangRen2021\4DNFIFBN2VO9.csv','NumHeaderLines',16);
spotT2 = readtable('U:\GenomeData\ByPublication\HuangRen2021\4DNFIKPGMZJ8.csv');

% GRCm38
IDs_is129 = traceT2.x_columns__Trace_ID(traceT2.Allele==129);
IDs_isCAST = traceT2.x_columns__Trace_ID( ~(traceT2.Allele==129));
spotT2.Trace_ID

idCAST = contains(spotT.Trace_ID,'CAST');
spotT_129 = spotT(~idCAST,:) ;
spotT_CAST = spotT(idCAST,:) ;

%% separate out the polymer and distance maps for each allele
curSpotT = spotT_129;
bars = unique(curSpotT.Chrom_Start);
nB = length(bars);
cellIDs = unique(curSpotT.Cell_ID_);
nC = length(cellIDs);
polys = nan(nB,3,nC);
for c=1:nC
    cellTable = curSpotT(curSpotT.Cell_ID_==cellIDs(c),:);
    cellTable.Chrom_Start;

    [~,ia,ib] = intersect(bars,cellTable.Chrom_Start);
    xyz = [cellTable.X(ib),cellTable.Y(ib),cellTable.Z(ib)];
    polys(:,:,c) = xyz;
end
polys_129 = polys;

curSpotT = spotT_CAST;
bars = unique(curSpotT.Chrom_Start);
nB = length(bars);
cellIDs = unique(curSpotT.Cell_ID_);
nC = length(cellIDs);
polys = nan(nB,3,nC);
for c=1:nC
    cellTable = curSpotT(curSpotT.Cell_ID_==cellIDs(c),:);
    cellTable.Chrom_Start;

    [~,ia,ib] = intersect(bars,cellTable.Chrom_Start);
    xyz = [cellTable.X(ib),cellTable.Y(ib),cellTable.Z(ib)];
    polys(:,:,c) = xyz;
end
polys_CAST = polys;

dmap_129 = PolyToDistMap(polys_129);
medMap_129 = nanmedian(dmap_129,3);
medMap_129(26,:) = nan;
medMap_129(:,26) = nan;

dmap_CAST = PolyToDistMap(polys_CAST);
medMap_CAST = nanmedian(dmap_CAST,3);
figure(1); clf; 
subplot(1,2,1); imagesc(medMap_129); colorbar; clim([0,350])
subplot(1,2,2); imagesc(medMap_CAST); colorbar; clim([0,350])
GetColorMap('redToWhiteSat')

scr = 30:32; % barcodes corresponding to the SCR in this dataset
medMap_CAST(scr,scr)
medMap_129(scr,scr)


%% Same as above but on the replicate data

idCAST = contains(spotT2.Trace_ID,'CAST');
spotT_129 = spotT2(~idCAST,:) ;
spotT_CAST = spotT2(idCAST,:) ;

curSpotT = spotT_129;
bars = unique(curSpotT.Chrom_Start);
nB = length(bars);
cellIDs = unique(curSpotT.Cell_ID_);
nC = length(cellIDs);
polys = nan(nB,3,nC);
for c=1:nC
    cellTable = curSpotT(curSpotT.Cell_ID_==cellIDs(c),:);
    cellTable.Chrom_Start;

    [~,ia,ib] = intersect(bars,cellTable.Chrom_Start);
    xyz = [cellTable.X(ib),cellTable.Y(ib),cellTable.Z(ib)];
    polys(:,:,c) = xyz;
end
polys_129 = polys;



curSpotT = spotT_CAST;
bars = unique(curSpotT.Chrom_Start);
nB = length(bars);
cellIDs = unique(curSpotT.Cell_ID_);
nC = length(cellIDs);
polys = nan(nB,3,nC);
for c=1:nC
    cellTable = curSpotT(curSpotT.Cell_ID_==cellIDs(c),:);
    cellTable.Chrom_Start;

    [~,ia,ib] = intersect(bars,cellTable.Chrom_Start);
    xyz = [cellTable.X(ib),cellTable.Y(ib),cellTable.Z(ib)];
    polys(:,:,c) = xyz;
end
polys_CAST = polys;

dmap_129 = PolyToDistMap(polys_129);
medMap_129_2 = nanmedian(dmap_129,3);
medMap_129_2(26,:) = nan;
medMap_129_2(:,26) = nan;

dmap_CAST = PolyToDistMap(polys_CAST);
medMap_CAST_2 = nanmedian(dmap_CAST,3);
figure(1); clf; 
subplot(1,2,1); imagesc(medMap_129_2); colorbar; clim([100,350]);
subplot(1,2,2); imagesc(medMap_CAST_2); colorbar; clim([100,350]);
GetColorMap('redToWhiteSat')

medMap_CAST_2(scr,scr)
medMap_129_2(scr,scr)

nanmean(dmap_129(scr,scr,:),3)

%% load Tonia data
 analysisFolder = 'Z:\Tonia\2018-05-24-mES_Bo_Sox2\Analysis\';
[polysT,mapsT] = CombineAllFits(analysisFolder,'dims','xyz','parallel',1); % ,'byFOV',true   

scrT = 29:31; % barcodes corresponding to the SCR

medMapT = nanmedian(mapsT,3);
figure(3); clf; imagesc(medMapT); colorbar; clim([100,350]); 
GetColorMap('redToWhiteSat');
medMapT(scrT,scrT)
%%
dax1 = ['Z:\Tonia\2018-05-24-mES_Bo_Sox2\Hyb_037\ConvZscan_',num2str(f,'%02d'),'.dax'];
im1 = ReadDax(dax1);
im1_d = im1(:,:,1:2:end);
im1_f = im1(:,:,2:2:end);

dax2 = ['Z:\Tonia\2018-05-24-mES_Bo_Sox2\Hyb_038\ConvZscan_',num2str(f,'%02d'),'.dax'];
im2 = ReadDax(dax2);
im2_d = im2(:,:,1:2:end);
im2_f = im2(:,:,2:2:end);

imO = IncreaseContrast(cat(3,max(im1_d,[],3),max(im2_d,[],3)),'high',.9999,'low',.75);
figure(1); clf; Ncolor(imO);
imO = IncreaseContrast(cat(3,max(im1_f,[],3),max(im2_f,[],3)),'high',.9999,'low',.75);
figure(1); clf; Ncolor(imO);

%% another data set (this data is a bit different, the nearest probes are quite 
analysisFolder ='Z:\Tonia\2020-06-10_mES_multi_TwistRNA_Sox2_5kb_DNA\DNA\ChrTracer3_Out\';
[polysT,mapsT2] = CombineAllFits(analysisFolder,'dims','xyz','parallel',1,'bins',52); % ,'byFOV',true   
whos mapsT
scrT2 = 38:41;


medMapT2 = nanmedian(mapsT2,3);
figure(3); clf; imagesc(medMapT2); colorbar; clim([100,450]); 
GetColorMap('redToWhiteSat');
medMapT2(scrT2,scrT2)

figure(3); clf; 
subplot(1,2,1); imagesc(medMapT(1:33,1:33)); colorbar; clim([100,450]); 
subplot(1,2,2); imagesc(medMapT2(8:40,8:40)); colorbar; clim([100,450]); 

figure(4); clf; imagesc(medMapT2(8:40,8:40) - medMapT(1:33,1:33)); colorbar; caxis([-300,300])
GetColorMap('RedWhiteBlueSat')

% analysisFolder ='Z:\Tonia\2020-06-10_mES_multi_TwistRNA_Sox2_5kb_DNA\DNA\ChrTracer3_Out\';
% [~,mapsXY] = CombineAllFits(analysisFolder,'dims','xy','parallel',1,'bins',52); % ,'byFOV',true   
% [~,mapsXZ] = CombineAllFits(analysisFolder,'dims','xz','parallel',1,'bins',52); % ,'byFOV',true   
% scrT = 29:32;
% medMapXY = nanmedian(mapsXY,3);
% medMapXZ = nanmedian(mapsXZ,3);
% figure(3); clf; 
% subplot(1,2,1); imagesc(medMapXY); colorbar; clim([100,450]); 
% subplot(1,2,2); imagesc(medMapXZ); colorbar; clim([100,450]); 
%% Mateo et al 
% selectReads = 1:38; hs= selectReads;
% sox2 = CombineFOVdataFromFolder('R:\2018-05-24-mES-Sox2\2018-05-26-mES-Sox2-CT2out\','selectReads',selectReads);
% bH = 25 ;% missed hybe.


%%
eDists = {  squeeze(dmap_129(30,31,:)),squeeze(mapsT(29,30,:))};
figure(1); clf;   BoxPlotCell(eDists);
figure(2); clf; violin(eDists,'bandwidth',.02,'variableNames',{'Huang 2021','This study'}); ylim([0,1000])

cellfun(@nanmedian, eDists)
cellfun(@nanmean, eDists)

%%
eNames = {'adjacent 5kb (Huang 2021)','adjacent 5kb (This study)','enhancer ends (Huang 2021)','enhancer ends 5kb (This study)'};
eDists = {  squeeze(dmap_129(30,31,:)),squeeze(mapsT(29,30,:)),  squeeze(dmap_129(30,32,:)),squeeze(mapsT(29,31,:)) };
figure(2); clf; violin(eDists,'bandwidth',.02,'variableNames',eNames,'plotMean',false,'medianStyle',{'r.'});
ylim([0,1000]); ylabel('center-to-center distance (nm)')
cellfun(@(x) sum(~isnan(x)),eDists)
cellfun(@nanmedian, eDists)
cellfun(@nanmean, eDists)

%% 
'Z:\Tonia\2018-05-24-mES_Bo_Sox2\'


im =cat(4,da1.data.im{[1,2]});
im(:,:,:,1) = 1.5*im(:,:,:,1);
figure(1); clf; ProjectIm4D(im(585:615,860:890,40:70,:));

im =cat(4,da1.data.im{[2,3]});
im(:,:,:,2) = 1.5*im(:,:,:,2);
figure(1); clf; ProjectIm4D(im(585:615,860:890,40:70,:));


