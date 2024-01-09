%% Figure 5 plots

saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\';
load([saveFolder,'RnaBurst_DnaCluster.mat'],'rnaBurst_dnaCluster','rnaGeneNames');
load([saveFolder,'RnaBurst_Maps.mat'],'mapsON','mapsOFF'); % these are very large files, will be slow
geneNames = rnaGeneNames;
%%  Violin plots per gene
dispGenes =1:6;  % [1,4,6,10];%  
hh =0;
cMax = 10; % max cluster size to test
figure(5); clf; figure(4); clf;
or = zeros(6,1); ci=zeros(6,2);
clustSizeAll = cell(nH,1);
rnaLevelsAll = cell(nH,1);
for h=dispGenes % 1:10 % h=2
    hh=hh+1;
    maxNascentRNA = cat(1,rnaBurst_dnaClust{h,:,1});
    maxClusterSize = cat(1,rnaBurst_dnaClust{h,:,2});

    clustSizeAll{h} = maxClusterSize;
    rnaLevelsAll{h} = maxNascentRNA;
    ndata = isnan(maxNascentRNA) | isnan(maxClusterSize);
    maxNascentRNA(ndata) = [];
    maxClusterSize(ndata) = [];
    rnaCnt = cell(1,cMax);
    for c=1:cMax
        rnaCnt{c}  = maxNascentRNA(maxClusterSize==c-1);
    end
    clustSize = cellfun(@length,rnaCnt);

    clusterTheta =6; %  find(clustSize>2,1,'last')-2
    clusterRNA = cat(1,rnaCnt{clusterTheta:end});
    isolateRNA = cat(1,rnaCnt{1});
    isClustered = [true(length(clusterRNA),1); false(length(isolateRNA),1)];
    isTranscribing = [clusterRNA>1000; isolateRNA>1000];
    [or(hh),ci(hh,:)] = OddsRatioCI(isTranscribing,isClustered,'cI',.75);  

    % a = sum(isClustered & isTranscribing)
    % b = sum(isClustered & ~isTranscribing)
    % c= sum(~isClustered & isTranscribing)
    % d= sum(~isClustered & ~isTranscribing)
    % (a/(a+b)) / (c/(c+d))
    for c=1:length(rnaCnt)
        rnaDis = rnaCnt{c} ;
        if sum(rnaDis)== 0
            rnaDis = cat(1,rnaDis,1);
        end
        if length(rnaDis) <=2
            rnaDis = [0,1]';
        end
        rnaCnt{c} = rnaDis;
    end
    figure(4); subplot(3,2,hh); 
    semilogy(clustSize,'.','color',[.7 .2 .2]);    ylabel('n alleles');
    hold on; bar(clustSize);
    title(geneNames{h});
    figure(5); subplot(3,2,hh); cla; %  yyaxis left; 
    try
    violin(rnaCnt,'plotMean',false,'lineColor','none','bandwidth',.05,'faceColor',[1 .6 .6]); ylim([0,6000]);
    catch er
        warning(['error plotting data from ',num2str(h)]);
        warning(er.message)
    end
    ylabel('nascent transcription (au)')
    xlabel('SE community size (600 nm)')
    [c1,p1] = corr(maxClusterSize,maxNascentRNA);
    cf = polyfit(maxClusterSize,maxNascentRNA,1);
    xdat = maxClusterSize;
     yy = polyval(cf,xdat);
     hold on; plot(xdat,yy,'-','color',[1 .4 .4]);
    title([geneNames{h},'  r=',num2str(c1,2),' p=',num2str(p1,2)]);
    xlim([0,10]);
end
figure(6); clf; semilogy(or,'.','MarkerSize',30); ylim([1/10 10]);
hold on; plot(ci(:,1),'k.'); hold on; plot(ci(:,2),'k.');
plot([1,6],[1,1],'k--'); ylabel('Odds Ratio')




%% combine all in same violin plot
rnaCnt = cell(1,cMax);
maxNascentRNA = cat(1,rnaLevelsAll{:});
maxClusterSize = cat(1,clustSizeAll{:})
for c=1:cMax
    rnaCnt{c}  = maxNascentRNA(maxClusterSize==c-1);
end
skip = isnan(maxNascentRNA) | isnan(maxClusterSize); 
maxNascentRNA(skip)=[];
maxClusterSize(skip)=[];
figure(5); clf;
clustSize = cellfun(@length,rnaCnt);
semilogy(clustSize,'.','color',[.7 .2 .2]);    ylabel('n alleles');
hold on; bar(clustSize);
% 
figure(4); clf;
violin(rnaCnt(1:10),'plotMean',false,'lineColor','none','bandwidth',.05,'faceColor',[1 .6 .6]); ylim([0,6000]);
ylabel('nascent transcription (au)')
xlabel('SE community size (600 nm)')
[c1,p1] = corr(maxClusterSize,maxNascentRNA);
 cf = polyfit(maxClusterSize,maxNascentRNA,1);
xdat = maxClusterSize;
yy = polyval(cf,xdat);
hold on; plot(xdat,yy,'-','color',[1 .4 .4]);
xlim([0,10]);
title(['All Genes','  r=',num2str(c1,2),' p=',num2str(p1,2)]);

clusterTheta =6; %  find(clustSize>2,1,'last')-2
clusterRNA = cat(1,rnaCnt{clusterTheta:end});
isolateRNA = cat(1,rnaCnt{1:5});
isClustered = [true(length(clusterRNA),1); false(length(isolateRNA),1)];
isTranscribing = [clusterRNA>1000; isolateRNA>1000];
[or_all,ci_all] = OddsRatioCI(isTranscribing,isClustered,'cI',.75)
[rr_all,rci_all] = RelativeRisk(isTranscribing,isClustered,'cI',.75);  
title(['All Genes','  r=',num2str(c1,2),' p=',num2str(p1,2),...
    '   OR=',num2str(or_all,2),' (',num2str(ci_all(1),2),':',num2str(ci_all(2),2),')']);


%% effect of cluster size threshold on Odds Ratio
%% combine all
cTheta = 1:12
or_t = zeros(12,1);
ci_t = zeros(12,2);
for t=1:12
    rnaCnt = cell(1,12);
    maxNascentRNA = cat(1,rnaLevelsAll{:});
    maxClusterSize = cat(1,clustSizeAll{:});
    for c=1:12
        rnaCnt{c}  = maxNascentRNA(maxClusterSize==c-1);
    end
    clusterRNA = cat(1,rnaCnt{cTheta(t):end});
    isolateRNA = cat(1,rnaCnt{1}); % :clusterTheta-1
    isClustered = [true(length(clusterRNA),1); false(length(isolateRNA),1)];
    isTranscribing = [clusterRNA>1000; isolateRNA>1000];
    [or_t(t),ci_t(t,:)] = OddsRatioCI(isTranscribing,isClustered,'cI',.75,'iters',100);
end
or_t(isnan(or_t)) = 1;
figure(2); clf; plot(cTheta',or_t,'linewidth',2,'color',lines(1)); hold on;
plot(cTheta',ci_t(:,1),'linewidth',1,'color',[.5 .5 .5]); 
plot(cTheta',ci_t(:,2),'linewidth',1,'color',[.5 .5 .5]); 
xlabel('clusters size'); ylabel('Odds Ratio for burst if clustered')

%% correlations in bursting between genes
% this is fun but not too relevant to SE clustering
% Sox2 and Myc show correlated bursts. 
% Otherwise most genes are pretty independent (ever so slightly correlated)

burstSize = rnaBurst_dnaClust(:,:,1); % burst size hyb x fov 
rnaBurst_dnaClust(:,:,2); % cluster size

burstFov = cell(nFOV,1);
for f=1:nFOV
burstFov{f} = cat(2,burstSize{:,f});
end
burstFovs = cat(1,burstFov{:});
figure(1); clf; imagesc(burstFovs);
burstFovs = burstFovs(:,dispGenes);
nG = length(dispGenes);
bC = nan(nG,nG);
for h=1:nG
    for g=1:nG
        b1 = burstFovs(:,h);
        b2 = burstFovs(:,g);
        skip = isnan(b1) | isnan(b2);
        bC(h,g) = corr(b1(~skip),b2(~skip));
    end
end
figure(1); clf; imagesc(bC); colorbar; clim([-1 1]);
GetColorMap('RedWhiteBlue');
set(gca,'xticklabel',geneNames(dispGenes),'yTickLabel',geneNames(dispGenes))
title('Co-burst Correlation')


%% maps  ON vs. Off distance maps
figure(41); clf;
figure(31); clf; 
for h=1:10
    % h=2;
on = nanmedian(mapsON{h},3);
off = nanmedian(mapsOFF{h},3);

figure(31);
subplot(2,10,h); imagesc(on); colorbar;  clim([.5 5.5]); title([geneNames{h},'  on n=',num2str(sum(~isnan(mapsON{h}(10,11,:))))]);
subplot(2,10,10+h); imagesc(off); colorbar; clim([.5 5.5]); title([geneNames{h},'  off n=',num2str(sum(~isnan(mapsOFF{h}(10,11,:))))]);

figure(2); clf; PlotCorr(on(:),off(:),'log',false);
idx = find((on./off) <.9);
[cy1,cx1] = ind2sub(size(on),idx );
cy1(cy1>44) = cy1(cy1>44)-44;
cx1(cx1>44) = cx1(cx1>44)-44;
[chrTableHyb.Yo_gene(cy1), chrTableHyb.Yo_gene(cx1)]

idx = find((on - off) >.4);
[cy1,cx1] = ind2sub(size(on),idx );
cy1(cy1>44) = cy1(cy1>44)-44;
cx1(cx1>44) = cx1(cx1>44)-44;
[chrTableHyb.Yo_gene(cy1), chrTableHyb.Yo_gene(cx1)]


figure(41); 
subplot(2,5,h); imagesc(on-off); colorbar; caxis([-1 1]);
idx = chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,geneNames{h}));
title([geneNames{h}, '  ',num2str(idx')]);
end

figure(41); colorbar; GetColorMap('RedWhiteBlueSat');
figure(31); GetColorMap('redToWhiteSat')

%% contact maps (on vs off)
figure(40); clf;
figure(30); clf; 
for h=1:10
    % h=2;
on = ContactFrac(mapsON{h},'threshold',.6);
off =  ContactFrac(mapsOFF{h},'threshold',.6);
figure(30);
subplot(10,2,h); imagesc(on); colorbar;  clim([0 .3]); title(['on n=',num2str(sum(isnan(mapsON{h}(1,1,:))))]);
subplot(10,2,10+h); imagesc(off); colorbar;   clim([0 .3]); title(['off n=',num2str(sum(isnan(mapsOFF{h}(1,1,:))))]);

figure(40); 
subplot(2,5,h); imagesc(log2(on./off)); colorbar; caxis([-6 6]);
idx = chrTableHyb.chrUID(StringFind(chrTableHyb.Yo_gene,geneNames{h}));
title([geneNames{h}, '  ',num2str(idx')]);
end

figure(40); colorbar; GetColorMap('RedWhiteBlueSat');
figure(30); colormap(flipud(GetColorMap('redToWhiteSat')))





%% Coburst burst-size analysis in clusters

load([saveFolder,'RNA_DNA_coburst.mat'],'loci_RNA_DNA_coburst');

nG = length(dispGenes)
cmap = hsv(nG);
nCells = 5000;
txLocs = nan(nCells*2,nG,4);

cc = 0;
for f=1:nFOV
    figure(1); clf;
    sz= cellfun(@size,loci_RNA_DNA(:,f),'UniformOutput',false);
    sz = cat(1,sz{:});
    nCells_in_FOV = max(sz(:,1))
    for c=1:nCells_in_FOV
        cc=cc+1;
        g=0;
        for h=1:nH     
            if ~isempty(loci_RNA_DNA{h,f})
                g=g+1;
                rna_xy = loci_RNA_DNA{h,f}{c,1};
                dna_xy = loci_RNA_DNA{h,f}{c,2};
                if ~isempty(rna_xy)
                    txLocs([cc,nCells+cc],g,:) = rna_xy(:,1:4);
                end
                if ~isempty(rna_xy)
                    plot(rna_xy(:,1),rna_xy(:,2),'x','color',cmap(g,:),'MarkerSize',15); hold on;
                end
                plot(dna_xy(:,1),dna_xy(:,2),'o','color',cmap(g,:)); hold on;
            end
        end
    end
end

%%
nCells = cc;
contactsPerRNA = nan(nCells,nG);
txBurstPerRNA  = nan(nCells,nG);

for c=1:nCells % note n-cells has been updated
dMap = squareform(pdist([txLocs(c,:,1)',txLocs(c,:,2)']));
% figure(1); clf; imagesc(dMap); colormap(jet);
txBurstPerRNA(c,:) = txLocs(c,:,4);
hasBurst =  txLocs(c,:,4) > 1000;
dMap(~hasBurst,:) = NaN;
dMap(:,~hasBurst) = NaN;
contactsPerRNA(c,:) = nansum(dMap<.6);
end

txBurstPerRNA(  ~(txBurstPerRNA>1000) ) = NaN;
figure(1); clf; PlotCorr(contactsPerRNA(:),txBurstPerRNA(:),'log',false)
v1 = contactsPerRNA(:);
v2 = txBurstPerRNA(:);
skip = isnan(v1) | isnan(v2); 
[r,p] = corr(v1(~skip),v2(~skip));
cf = polyfit(v1(~skip),v2(~skip),1);
yy = polyval(cf,v1(~skip));
hold on; plot(v1(~skip),yy,'-','color',[1 .4 .4]);

%%
txBurstPerRNAall = txBurstPerRNA(:);
contactsPerRNAall = contactsPerRNA(:);
skip = isnan(txBurstPerRNAall) | isnan(contactsPerRNAall);
txBurstPerRNAall(skip) = [];
contactsPerRNAall(skip) =[];
rnaCnt = {};
for c=1:3
    rnaCnt{c}  = txBurstPerRNAall(contactsPerRNAall==c);
end

clusterRNA = cat(1,rnaCnt{2:end});
isolateRNA = cat(1,rnaCnt{1});
isClustered = [true(length(clusterRNA),1); false(length(isolateRNA),1)];
isTranscribing = [clusterRNA>3500; isolateRNA>3500];  % 4000
sum(isTranscribing)
[or_all,ci_all] = OddsRatioCI(isTranscribing,isClustered,'cI',.75)
[rr_all,rci_all] = RelativeRisk(isTranscribing,isClustered,'cI',.75); 

figure(1); clf; 
 clustSize = cellfun(@length,rnaCnt);
semilogy(clustSize,'r.');    ylabel('log10(n alleles)');
hold on; bar(clustSize);
% 
figure(2); clf; 
v1 = contactsPerRNA(:);
v2 = txBurstPerRNA(:);
skip = isnan(v1) | isnan(v2); 
[r,p] = corr(v1(~skip),v2(~skip));
cf = polyfit(v1(~skip),v2(~skip),1);
yy = polyval(cf,v1(~skip));

     
violin(rnaCnt(~cellfun(@isempty,rnaCnt)),'plotMean',false,'lineColor','none','bandwidth',.04,'faceColor',[1 .6 .6]);%

     hold on; plot(v1(~skip),yy,'-','color',[1 .4 .4]);
title(['All Genes','  r=',num2str(r,2),' p=',num2str(p,2),...
    '   OR=',num2str(or_all,2),' (',num2str(ci_all(1),2),':',num2str(ci_all(2),2),')']);
ylabel('nascent RNA');
  ylim([0,6000]); xlim([.5,4.5])

%%
thetas = 0:100:5000;
nT = length(thetas)
or_t = zeros(nT,1);
ci_t = zeros(nT,2);
for t=1:nT
    clusterRNA = cat(1,rnaCnt{2:end});
    isolateRNA = cat(1,rnaCnt{1});
    isClustered = [true(length(clusterRNA),1); false(length(isolateRNA),1)];
    isTranscribing = [clusterRNA>thetas(t); isolateRNA>thetas(t)];  % 4000
    sum(isTranscribing)
    [or_t(t),ci_t(t,:)] = OddsRatioCI(isTranscribing,isClustered,'cI',.75,'iters',100);
end
or_t(isnan(or_t)) = 1;
figure(3); clf; plot(thetas',or_t,'linewidth',2,'color',lines(1)); hold on;
plot(thetas',ci_t(:,1),'linewidth',1,'color',[.5 .5 .5]); 
plot(thetas',ci_t(:,2),'linewidth',1,'color',[.5 .5 .5]); 
xlabel('"large" burst size'); ylabel('Odds Ratio for "large burst" if clustered')

