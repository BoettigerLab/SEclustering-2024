%% Fig 2 panels

% Run Process_Data_Tables_To_Matrices first 


%% gene names
nearestGene = chrTableHyb.Yo_gene;
nearestAllele = [nearestGene; nearestGene];
find(StringFind(nearestGene,'Sox2','boolean',true))


% idx_pp = contains(nearestGene,{'Sox2','Mycn','Tbx3','Pou5f1','Pvt1','Nanog','Prdm14'});
%   gene names to display
show_genes = [74:76,133:137,242:245] ; % Sox2, Nanog, Mycn SEs + 1 downstream control
idx_pp = false(length(nearestGene),1);
idx_pp(show_genes) = true;

%%
cFrac = cFrac2; % 200 nm cutoff 
meanObs = nanmean(nObs{2});
meanObs = meanObs(1:376);  % validated already the two halves are identical in badhybes  
bH = [209:212 ,find(meanObs < .5*max(meanObs))];  % 209:212 had an error in the barcode updates  
bH = [bH, 48, 337]; % this pair is also stripey,
bH = [bH,bH+376];
% these bad hybes have low detection efficiency and often as not record
% noise instead of true signal. To avoid the resulting confusion, we simply
% exclude these SE from furtehr analysis. 

% some cross-reactive hybes - these ones don't give stripes, just a
% peculiar off-diagonal network of dots. We'll drop these from the analysis. 
cxr = [91,213,231,257,291,296,333]; 
cxr = [cxr,cxr+376]; % both alleles;

cleanMap = cFrac2{2};
cleanMap(cxr,cxr) = nan;  % does all combos of cxr with cxr
cleanMap = InterpMapNans(cleanMap,'badHybes',bH);
colormap(flipud(GetColorMap('redToWhiteK')))
cf200 = (cleanMap(1:376,1:376)+cleanMap( (376+1):2*376,(376+1):2*376))/2;
figure(100); clf; imagesc(cf200); caxis([0,.02]); colorbar;
colormap(flipud(GetColorMap('redToWhiteK'))); axis image;

cfN = InterpMapNans(cFrac2{3},'badHybes',bH,'window',4);
cfN = (cfN(1:376,1:376)+cfN( (376+1):2*376,(376+1):2*376))/2;
figure(101); clf; imagesc(cfN); caxis([0,.02]); colorbar;
colormap(flipud(GetColorMap('redToWhiteK'))); axis image;
title('NPC contact map 200 nm');

% log scale as well 
figure(2); clf; imagesc(log10(InterpMapNans(cFrac2{2},'badHybes',bH,'window',4))); colorbar;
colormap(flipud(GetColorMap('redToWhiteK'))); caxis([-3,0])

%% zoom in on some known contacts as controls
cFrac = cFrac2; 

figure(1); clf;
k=0; N = 3;  % 74 and 75 are sox2
k=k+1; subplot(N,2,k); imagesc(cFrac{2}(70:80,70:80)); caxis([0.01,.2]); colorbar; axis image; title('sox2'); axis off;
k=k+1;subplot(N,2,k); imagesc(cFrac{3}(70:80,70:80)); caxis([0.01,.2]); colorbar; axis image; axis off;
colormap(flipud(GetColorMap('redToWhiteK'))); 

% 135 is nanog
k=k+1; subplot(N,2,k); imagesc(cFrac{2}(130:140,130:140)); caxis([0.01,.2]); colorbar; axis image; title('nanog'); axis off;
k=k+1; subplot(N,2,k); imagesc(cFrac{3}(130:140,130:140)); caxis([0.01,.2]); colorbar; axis image; axis off;
colormap(flipud(GetColorMap('redToWhiteK'))); 

% 242 - 244 are Mycn
colormap(flipud(GetColorMap('redToWhiteK'))); 
k=k+1; subplot(N,2,k); imagesc(cFrac{2}(238:248,238:248)); caxis([0.01,.2]); colorbar; axis image; title('mycn'); axis off;
k=k+1;  subplot(N,2,k); imagesc(cFrac{3}(238:248,238:248)); caxis([0.01,.2]); colorbar; axis image; axis off;
colormap(flipud(GetColorMap('redToWhiteK'))); 

% report the actual values as well, so we can quote them in text.  
ig = find(StringFind(nearestGene,'Nanog','boolean',true))
cFrac{2}(ig-2:ig+2,ig-2:ig+2)
cFrac{3}(ig-2:ig+2,ig-2:ig+2)
figure(10); clf; imagesc( log2( cFrac{2}(ig-2:ig+2,ig-2:ig+2)./cFrac{3}(ig-2:ig+2,ig-2:ig+2)) );
GetColorMap('RedWhiteBlueSat'); colorbar; caxis([-3,3]);

ig = find(StringFind(nearestGene,'Sox2','boolean',true))
cFrac{2}(ig,ig)
cFrac{3}(ig,ig)
figure(10); clf; imagesc( log2( cFrac{2}(ig,ig)./cFrac{3}(ig,ig)) );
GetColorMap('RedWhiteBlueSat'); colorbar; caxis([-3,3]);

ig = find(StringFind(nearestGene,'Mycn','boolean',true))
cFrac{2}(ig,ig)
cFrac{3}(ig,ig)
figure(10); clf; imagesc( log2( cFrac{2}(ig,ig)./cFrac{3}(ig,ig)) );
GetColorMap('RedWhiteBlueSat'); colorbar; caxis([-3,3]);

ig = find(StringFind(nearestGene,'Tbx3','boolean',true))
cFrac{2}(ig,ig)
cFrac{3}(ig,ig)
figure(10); clf; imagesc( log2( cFrac{2}(ig,ig)./cFrac{3}(ig,ig)) );
GetColorMap('RedWhiteBlueSat'); colorbar; caxis([-3,3]);

ig = find(StringFind(nearestGene,'Hist','boolean',true));
figure(10); clf; 
imagesc( log2( cFrac{2}(ig,ig)./cFrac{3}(ig,ig)) );
GetColorMap('RedWhiteBlueSat'); colorbar; caxis([-3,3]);

%% Show that even cis-cis contacts are generally rare 
nR = 376;
% mask cis
cisID = false(nR,nR);
cisCoords = cell(23,1);
for c=1:23
idC = chrTableHyb.chrUID(chrTableHyb.chrNum==c);
cisID(idC,idC) = true;
cisCoords{c} = idC;
end
figure(4); clf; imagesc(cisID); colorbar;
cis200 = cf200(cisID);
cis200(cis200==1) = nan;

% show in log scale for readability
figure(2); clf;  loglog(1,1,'w.'); hold on;
[b,x] = hist(cis200(:),0:.005:1);
b = b./sum(b);
bar(x,b,1);
colormap('default'); xlim([.001,1])
xlabel('SE cis-contact frequency (300 nm)')
ylabel('freq')

figure(2); clf;  % also show in linear scale
[b,x] = hist(cis200(:),0:.0005:1);
b = b./sum(b);
bar(x,b,1);
colormap('default'); xlim([.001,1])
xlabel('SE cis-contact frequency (200 nm)')
ylabel('freq')

% the median cis-contact frequency, to quote in text
nanmedian(cis200(:))



%% compare fraction alone, vs. 1, 2,3+ contacts
scMap2= single(esMap(:,:,:) < .2); % .3 to 1
scObs2 = single(esMap(:,:,:) < inf); % .3 to 1

 f_alone =  nansum(nansum(scMap2,2)==1,3)./(nansum(nansum(scObs2,2)>0,3)); %#ok<*NANSUM>
 f_1 =  nansum(nansum(scMap2,2)==2,3)./(nansum(nansum(scObs2,2)>0,3));
 f_2 =  nansum(nansum(scMap2,2)==3,3)./(nansum(nansum(scObs2,2)>0,3));
 f_3 =  nansum(nansum(scMap2,2)>=4,3)./(nansum(nansum(scObs2,2)>0,3));


figure(30); clf; 
plot(f_alone(1:376),'o','color',[0 0 0]); hold on;  plot(f_alone(376+1:end),'+','color',[.4 .4 .4]); 
plot(f_1(1:376),'o','color',[1,0,0]);  hold on; plot(f_1(376+1:end),'+','color',[1,.4,.4]); 
plot(f_2(1:376),'o','color',[0,1,0]);  hold on; plot(f_2(376+1:end),'+','color',[.4,1,.4]); 
plot(f_3(1:376),'o','color',[0,0,1]);  hold on; plot(f_3(376+1:end),'+','color',[.4,.4,1]); 
xlabel('SE ID')
ylabel('fraction of cells')

text(find(bol_sel), f_alone(bol_sel),nearestGene(bol_sel),'color','b')

%% simulate uniform random distribution of spots as a control
nPts = 2*376;
sz = 12; % just for plotting
rm = 4; % cell radius in um, estimated from data
nEx = 4000; % number of examples for comparison

theta = rand(nPts,nEx)*2*pi;
phi = acos(1-2*rand(nPts,1));%  rand(nPts,1)*pi;
r = rm*rand(nPts,nEx).^(1/3);
x = r.*sin(phi).*sin(theta);
y = r.*sin(phi).*cos(theta);
z = r.*cos(phi);
c=1:nPts/2; % 1:nPts; % sel
c= ones(1,nPts/2)
c = [c,c];
figure(1); clf;
e=1;
scatter3(x(:,e),y(:,e),z(:,e),sz,c,'filled'); hsv(nPts); 
axis image; colorbar;
xlim([-rm ,rm]); ylim([-rm ,rm]); 
randMaps = zeros(nPts,nPts,nEx);
for e=1:nEx
    randMaps(:,:,e) = squareform(pdist([x(:,e),y(:,e),z(:,e)]));
end

%% 200 nm cut-off view cluster sizes 
contactCnt = squeeze(nansum(scMap2,1)); % genes x cells
obsCnt = squeeze(nansum(scObs2,1)); % genes x cells
cpc2 = nanmean(contactCnt,2)./(nansum(obsCnt>0,2)/size(obsCnt,2));  % Average cluster size per rSE   ('counts per cluster' or cpc)
cout = quantile(contactCnt,.95,2); % genes x cells
cout2 = cout./ (nansum(obsCnt>0,2)/size(obsCnt,2));  % Top (95-percentile) cluster sizes 


% average cluster size random 
scMap_rand= single(randMaps(:,:,:) < .2); % .3 to 1
scObs_rand = single(randMaps(:,:,:) < inf); % 
contactCnt = squeeze(nansum(scMap_rand,1)); % genes x cells
obsCnt = squeeze(nansum(scObs_rand,1)); % genes x cells
cpc_rand = nanmean(contactCnt,2)./(nansum(obsCnt>0,2)/size(obsCnt,2));
cout_rand = quantile(contactCnt,.95,2); % genes x cells
cout_rand = cout_rand./ (nansum(obsCnt>0,2)/size(obsCnt,2)); 


figure(7); clf; plot(cpc2(1:376),'.'); hold on; 
plot(cpc2(376+1:end),'+','MarkerSize',3); hold on; 
plot(cpc_rand(1:376));
xlabel('SE ID'); ylabel('ave. cluster size');
cpcP = cpc2(1:376);
text( find(cpcP>quantile(cpcP,.95)),cpcP(cpcP>quantile(cpcP,.95)), nearestGene(cpcP>quantile(cpcP,.95)),'color','r' );
text(find(idx_pp), cpcP(idx_pp),nearestGene(idx_pp),'color','b')
ylim([1,3]);


figure(8); clf; plot(cout2(1:376),'.'); hold on;   % 
plot(cout2(376+1:end),'+','MarkerSize',3); hold on;  % 
 plot(cout_rand(1:376));
xlabel('SE ID'); ylabel('95-percentile cluster size');
cout2b = cout2(1:376);
text( find(cout2b>quantile(cout2b,.95)),cout2b(cout2b>quantile(cout2b,.95)), nearestGene(cout2b>quantile(cout2b,.95)),'color','r' );
text(find(idx_pp), cout2b(idx_pp),nearestGene(idx_pp),'color','b')



%% same data on sorted plots

%% compare fraction alone, vs. 1, 2,3+ contacts
scMap2= single(esMap(:,:,:) < .2); % .3 to 1
scObs2 = single(esMap(:,:,:) < inf); % .3 to 1

 f_alone =  nansum(nansum(scMap2,2)==1,3)./(nansum(nansum(scObs2,2)>0,3)); %#ok<*NANSUM>
 f_1 =  nansum(nansum(scMap2,2)==2,3)./(nansum(nansum(scObs2,2)>0,3));
 f_2 =  nansum(nansum(scMap2,2)==3,3)./(nansum(nansum(scObs2,2)>0,3));
 f_3 =  nansum(nansum(scMap2,2)>=4,3)./(nansum(nansum(scObs2,2)>0,3));

[~,fsort] = sort(f_alone(1:376));
figure(3); clf; m=2;
plot(f_alone(fsort),'o','color',[0 .6 .9],'MarkerSize',m); hold on; f0b = f_alone(376+1:end); plot(f0b(fsort),'>','color',[0 .6 .9],'MarkerSize',m); 
plot(f_1(fsort),'o','color',[.7,.6,0],'MarkerSize',m);  hold on; f1b = f_1(376+1:end); plot(f1b(fsort),'>','color',[.87,.6,0],'MarkerSize',m); 
plot(f_2(fsort),'o','color',[.9,.1,.1],'MarkerSize',m);  hold on; f2b = f_2(376+1:end); plot(f2b(fsort),'>','color',[.9,.1,.1],'MarkerSize',m); 
plot(f_3(fsort),'o','color',[1,0,1],'MarkerSize',m);  hold on; f3b = f_3(376+1:end); plot(f3b(fsort),'>','color',[1,0,1],'MarkerSize',m); 
xlabel('SE ID')
ylabel('fraction of cells')
% label sel genes
idx_sel = [74:76,134:138,241:245]; % sox, nanog, mycn context
bol_sel = false(1,376);
bol_sel(idx_sel) = true;
bol_sort = bol_sel(fsort);
fa = f_alone(fsort);
ng = nearestGene(fsort);
% plot(find(bol_sort),fa(bol_sort)-.05,'+','color',.034*ones(1,3));
% text(2+find(bol_sort),fa(bol_sort)-.05 ,ng(bol_sort),'color','b');
plot(find(bol_sort),ones(sum(bol_sort),1),'+','color','b');  % 
text(2+find(bol_sort),.95*ones(sum(bol_sort),1),ng(bol_sort),'color','b')
ng(1:10)
ng(end-10:end)
N=sum(~isnan(f_alone(1:376)))
plot(1:5,ones(1,5),'+','color','b');  % 
text(1:5,.92*ones(1,5),ng(1:5),'color','b');  % 
plot(N-4:N,ones(1,5),'+','color','b');  % 
text(N-4:N,.97*ones(1,5),ng(N-4:N)','color','b');  % 


[~,fsort] = sort(cpc2(1:376));
figure(7); clf; plot(cpc2(fsort),'o','MarkerSize',m); hold on; 
cpc2b = cpc2(376+1:end);
plot(cpc2b(fsort),'>','MarkerSize',m); hold on; 
plot(sort(cpc_rand(1:376)),'+','MarkerSize',m);
xlabel('SE ID'); ylabel('ave. cluster size');
cpcP = cpc2(fsort);
% text( find(cpcP>quantile(cpcP,.95)),cpcP(cpcP>quantile(cpcP,.95)), nearestGene(cpcP>quantile(cpcP,.95)),'color','r' );
% text(find(idx_pp), cpcP(idx_pp),nearestGene(idx_pp),'color','b')
ylim([1,2]);
fa = cpc2b(fsort);
ng = nearestGene(fsort);
bol_sort = bol_sel(fsort);
plot(find(bol_sort),fa(bol_sort)-.05,'+','color',.034*ones(1,3));
text(2+find(bol_sort),fa(bol_sort)-.05 ,ng(bol_sort),'color','b');
plot(1:5,fa(1:5),'+','color','r');  % 
text(1:5,.92*fa(1:5),ng(1:5),'color','r');  % 
plot(N-4:N,fa(N-4:N),'+','color','r');  % 
text(N-4:N,fa(N-4:N),ng(N-4:N),'color','r');  % 


[~,fsort] = sort(cout2(1:376));
figure(8); clf; plot(cout2(fsort),'o','MarkerSize',m); hold on;   % 
coutb = cout2(376+1:end);
plot(coutb(fsort),'>','MarkerSize',m); hold on;  % 
 plot(sort(cout_rand(1:376)),'+','MarkerSize',m);
xlabel('SE ID'); ylabel('95-percentile cluster size');
cout2b = cout2(fsort);
% text( find(cout2b>quantile(cout2b,.95)),cout2b(cout2b>quantile(cout2b,.95)), nearestGene(cout2b>quantile(cout2b,.95)),'color','r' );
% text(find(idx_pp), cout2b(idx_pp),nearestGene(idx_pp),'color','b')
fa = cout2(fsort);
ng = nearestGene(fsort);
bol_sort = bol_sel(fsort);
plot(find(bol_sort),fa(bol_sort)-.05,'+','color',.034*ones(1,3));
text(2+find(bol_sort),fa(bol_sort)-.05 ,ng(bol_sort),'color','b');
plot(1:5,fa(1:5),'+','color','r');  % 
text(1:5,.92*fa(1:5),ng(1:5),'color','r');  % 
plot(N-4:N,fa(N-4:N),'+','color','r');  % 
text(N-4:N,fa(N-4:N),ng(N-4:N),'color','r');  % 
%% distance correlation
dis = chrTableHyb.chrNum*1e9 + chrTableHyb.start;
disM = squareform(pdist([dis,dis]));
disM(disM>1e9) = nan;
figure(4); clf; imagesc(disM); colorbar;

cntFreq = ContactFrac(esMap(1:376,1:376,:),'threshold',.2);

figure(4); clf; PlotCorr(disM(:),cntFreq(:));
xlabel('bp distance (cis)');
ylabel('med. contact freq');



v1 = log10(disM(:));
v2 = log10(cntFreq(:));
skip = isnan(v1) | isnan(v2);
skip = skip | isinf(v1) | isinf(v2);
v1(skip) = []; v2(skip) = []; 
figure(4); clf; hexscatter(v1,v2,'res',30);
cp = corr(v1,v2);
title(['r=',num2str(cp,3)]);
colormap(flipud(GetColorMap('magma'))); colorbar;
xlabel('log10  bp distance (cis)');
ylabel('log10  med. contact freq (<200 nm)');
set(gca,'color','w')

