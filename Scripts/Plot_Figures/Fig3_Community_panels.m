%% Fig 3 Community analysis

% Run Process_Data_Tables_To_Matrices to load the data first. 

%% gene names
nearestGene = chrTableHyb.Yo_gene;
nearestAllele = [nearestGene; nearestGene];
show_genes = [74:76,133:137,242:245] ; % Sox2, Nanog, Mycn SEs + 1 downstream control
idx_pp = false(length(nearestGene),1);
idx_pp(show_genes) = true;


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
c = [c,c];
figure(1); clf;
e=2;
scatter3(x(:,e),y(:,e),z(:,e),sz,c,'filled'); hsv(nPts); 
axis image; colorbar;
xlim([-rm ,rm]); ylim([-rm ,rm]); 
randMaps = zeros(nPts,nPts,nEx);
for e=1:nEx
    randMaps(:,:,e) = squareform(pdist([x(:,e),y(:,e),z(:,e)]));
end

% averag cluster size random 
scMapRand= single(randMaps(:,:,:) < .6); % .3 to 1
scObsRand = single(randMaps(:,:,:) < inf); % 
contactCntRnd = squeeze(nansum(scMapRand,1)); % genes x cells
obsCnt = squeeze(nansum(scObsRand,1)); % genes x cells
cpc_rand = nanmean(contactCntRnd,2)./(nansum(obsCnt>0,2)/size(obsCnt,2));
cout_rand = quantile(contactCntRnd,.95,2); % genes x cells
cout_rand = cout_rand./ (nansum(obsCnt>0,2)/size(obsCnt,2)); 


%% 600

scMap= single(esMap(:,:,:) < .6); % .3 to 1
scObs = single(esMap(:,:,:) < inf); % 
contactCnt = squeeze(nansum(scMap,1)); % genes x cells
obsCnt = squeeze(nansum(scObs,1)); % genes x cells
cpc6 = nanmean(contactCnt,2)./(nansum(obsCnt>0,2)/size(obsCnt,2));
cout = quantile(contactCnt,.95,2); % genes x cells
cout = cout./ (nansum(obsCnt>0,2)/size(obsCnt,2)); 
cout6 = cout;

%% frac alone compare fraction alone, vs. 1, 2,3+ contacts

 f_alone =  nansum(nansum(scMap,2)==1,3)./(nansum(nansum(scObs,2)>0,3)); %#ok<*NANSUM>
 f_1 =  nansum(nansum(scMap,2)==2,3)./(nansum(nansum(scObs,2)>0,3));
 f_2 =  nansum(nansum(scMap,2)==3,3)./(nansum(nansum(scObs,2)>0,3));
 f_3 =  nansum(nansum(scMap,2)>=4,3)./(nansum(nansum(scObs,2)>0,3));


figure(3); clf; 
plot(f_alone(1:376),'o','color',[0 0 0]); hold on;  plot(f_alone(376+1:end),'+','color',[.4 .4 .4]); 
plot(f_1(1:376),'o','color',[1,0,0]);  hold on; plot(f_1(376+1:end),'+','color',[1,.4,.4]); 
plot(f_2(1:376),'o','color',[0,1,0]);  hold on; plot(f_2(376+1:end),'+','color',[.4,1,.4]); 
plot(f_3(1:376),'o','color',[0,0,1]);  hold on; plot(f_3(376+1:end),'+','color',[.4,.4,1]); 
xlabel('SE ID')
ylabel('fraction of cells')

%%




[cpcS,idx] = sort(cpc6(1:376),'ascend')
 figure(7); clf; plot(cpcS,'.'); hold on; 
xlabel('SE ID'); ylabel('ave. community size');
nG = nearestGene(idx); 
idx_ppS = idx_pp(idx);
text( find(cpcS>quantile(cpcS,.975)),cpcS(cpcS>quantile(cpcS,.975)), nG(cpcS>quantile(cpcS,.975)),'color','r' );
text(find(idx_ppS), cpcS(idx_ppS),nG(idx_ppS),'color',[.5 0 1]);
numS = cellstr(num2str(find(idx_ppS)));
text(find(idx_ppS), cpcS(idx_ppS)-.5, numS,'color',[.5 0 1]);
nG(cpcS<quantile(cpcS,.025));
text( find(cpcS<quantile(cpcS,.025)),cpcS(cpcS<quantile(cpcS,.025)), nG(cpcS<quantile(cpcS,.025)),'color','b' );
text( find(cpcS<quantile(cpcS,.025)),cpcS(cpcS<quantile(cpcS,.025))-.5, numS,'color','b' );
plot(sort(cpc_rand(1:376)),'.');



[cpcS,idx] = sort(cout6(1:376),'ascend')
 figure(8); clf; plot(cpcS,'.'); hold on; 
xlabel('SE ID'); ylabel('95-percentile cluster size');
nG = nearestGene(idx); idx_ppS = idx_pp(idx);
text( find(cpcS>quantile(cpcS,.975)),cpcS(cpcS>quantile(cpcS,.975)), nG(cpcS>quantile(cpcS,.975)),'color','r' );
text(find(idx_ppS), cpcS(idx_ppS),nG(idx_ppS),'color',[.5 0 1])
text( find(cpcS<quantile(cpcS,.025)),cpcS(cpcS<quantile(cpcS,.025)), nG(cpcS<quantile(cpcS,.025)),'color','b' );
 plot(sort(cout_rand(1:376)),'.');





%% entropy

N = size(esMap,1);
theta = 0.6;
pmap = nansum(esMap<theta,3)./nansum(esMap<inf,3);
entrop = nan(N,1);
for e = 1:N
    idx = 1:N;
    idx(e) = [];
    pe = pmap(e,idx);
    pe = pe./nansum(pe);
    entrop(e) = -nansum(pe.*log2(pe));
end
entrop(entrop==0) = nan;
[esort,idx] = sort(entrop(1:376),'ascend');
figure(6); clf; plot(esort);
nG = nearestGene(idx); idx_ppS = idx_pp(idx);
text( find(esort>quantile(esort,.975)),esort(esort>quantile(esort,.975)), nG(esort>quantile(esort,.975)),'color','r' );
text(find(idx_ppS), esort(idx_ppS),nG(idx_ppS),'color',[.5 0 1]);
numS = cellstr(num2str(find(idx_ppS)));
text(find(idx_ppS), esort(idx_ppS)-.5, numS,'color',[.5 0 1]);
text( find(esort<quantile(esort,.025)),esort(esort<quantile(esort,.025)), nG(esort<quantile(esort,.025)),'color','b' );
ylabel('community entropy')



%% 3 way absolute interaction

% % loading is much faster than re-running
saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\';
prevCalc = [saveFolder,'trip600nm_alltrans.mat'];
if exist(prevCalc,'file')
	load(prevCalc,'sort_triples','sort_labABC')
else  
    theta = 0.6; 
    nB = size(esMap,1); 
    isABC = nan(nB,nB,nB);
    hasABC = nan(nB,nB,nB);
    labelABC= nan(nB,nB,nB,3);
    for a=1:nB
        for b=a+1:nB
            for c=b+1:nB
                if c~=a && c~= b && a~=b
                    isABC(a,b,c) =sum( esMap(a,b,:) < theta & esMap(b,c,:) < theta);
                    hasABC(a,b,c) =sum( esMap(a,b,:) < inf & esMap(b,c,:) < inf);
                    labelABC(a,b,c,:) = [a,b,c];
                end
            end
        end
        disp(a/nB);
    end
    
    
    trip = isABC./hasABC;
    keep = ~isnan(trip);
    labA = labelABC(:,:,:,1); labB = labelABC(:,:,:,2); labC = labelABC(:,:,:,3);
    labA = labA(keep); labB = labB(keep); labC = labC(keep);
    
    
    trip = trip(keep);
    [sort_triples,trip_idx] = sort(trip);
    labABC = [labA(:),labB(:),labC(:)];
    sort_labABC = labABC(trip_idx,:);
    figure(5); clf; plot(sort_triples);
    sort_triples(end-10:end)
    sort_labABC(end-10:end,:)
end
% saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\';
% save([saveFolder,'trip600nm_alltrans.mat'],'sort_triples','sort_labABC')

nearestGene2 = cat(1,nearestGene,nearestGene);
chrTable2x = cat(1,chrTableHyb,chrTableHyb);

figure(7); clf; semilogy(sort_triples);
xlabel('triplet ID (sorted)')
ylabel('frequency')

[nTrip,~] = size(sort_labABC)
idx_Nanog = find(ismember(sort_labABC,[134,135,136],'rows')); % Nanog
hold on; text(idx_Nanog,sort_triples(idx_Nanog),'Nanog');
idx_mycN = find(ismember(sort_labABC,[242,243,244],'rows')); % mycn
hold on; text(idx_mycN,sort_triples(idx_mycN),'MycN');
idx_sox2 = find(ismember(sort_labABC,[74,75,76],'rows')); % sox2
hold on; text(idx_sox2,sort_triples(idx_sox2),'sox2');



%% cooperativity
% this is slow to recompute, much faster to load
% therfore it is better to load a pre-computed version if it exists
try
    saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\';
    load([saveFolder,'Cooperativity_600nm.mat'],'obsVexp_es','obsVexp_np','labelABC');
catch
    theta = 0.6; 
    nB = 376;%   size(distMapCell{2},1);   % the full 750 is impossible (would be ~10x longer, and this is already ~1 week?)
    %    I suppose a small parfor with a 3x parpool might have been alright? 
    obsVexp_es = nan(nB,nB,nB);
    obsVexp_np = nan(nB,nB,nB);
    labelABC= nan(nB,nB,nB,3);
    for a=1:nB
        for b=a+1:nB
            for c=b+1:nB
                if c~=a && c~= b && a~=b
                    isABC =sum( esMap(a,b,:) < theta & esMap(b,c,:) < theta);
                    hasABC =sum( esMap(a,b,:) < inf & esMap(b,c,:) < inf);
                    obs = isABC./hasABC;
                    ab = sum( esMap(a,b,:) < theta)./sum( esMap(a,b,:) < inf);
                    bc = sum( esMap(b,c,:) < theta)./sum( esMap(b,c,:) < inf);
                    exp = ab.*bc;
                    if obs == 0  || exp==0
                        exp = exp+1;
                        obs = obs+1;
                    end
                    obsVexp_es(a,b,c) = obs./exp;
    
                    isABC =sum( npMap(a,b,:) < theta & npMap(b,c,:) < theta);
                    hasABC =sum( npMap(a,b,:) < inf & npMap(b,c,:) < inf);
                    obs = isABC./hasABC;
                    ab = sum( npMap(a,b,:) < theta)./sum( npMap(a,b,:) < inf);
                    bc = sum( npMap(b,c,:) < theta)./sum( npMap(b,c,:) < inf);
                    exp = ab.*bc;
                    if obs == 0  || exp==0
                        exp = exp+1;
                        obs = obs+1;
                    end
                    obsVexp_np(a,b,c) = obs./exp;
                    labelABC(a,b,c,:) = [a,b,c];
                end
            end
        end
        disp(a/nB);
    end
end

figure(1); clf;
subplot(2,1,1);
r = log2(obsVexp_es(:)); r(isnan(r)) = []; r(r==0) = []; 
hist( r , -30:.1:30);
title('obs/exp ESC 3-way SE <600nm');
xlim([-8,8]); xlabel('log2(obs/exp) 3-way cluster');
ylim([0,1e5])

r = log2(obsVexp_np(:));  r(isnan(r)) = []; r(r==0) = [];
subplot(2,1,2); 
hist(  r, -30:.1:30);
title('obs/exp NPC 3-way SE <600nm');
xlim([-8,8]); xlabel('log2(obs/exp) 3-way cluster')
ylim([0,1e5])

%%
figure(1); clf;
r = log2(obsVexp_es(:)); r(isnan(r)) = []; r(r==0) = []; 
hist( r , -30:.1:30);
title('obs/exp ESC 3-way SE <600nm');
xlim([-8,8]); xlabel('log2(obs/exp) 3-way cluster');
ylim([0,1e5])

%%
aveCoop = squeeze(nanmean(obsVexp_es,1));
aveCoop(isnan(aveCoop)) = 0;
aveCoop = aveCoop + triu(aveCoop,1)';
aveCoop(aveCoop==0)=nan;
figure(2); clf; imagesc(log2(aveCoop)); colorbar;
GetColorMap('RedWhiteBlueSat'); colorbar; clim([-6 6]);
%%
idx_sel = [74:76,134:138,241:245]; % sox, nanog, mycn +nearby SEs
bol_sel = false(1,376);
bol_sel(idx_sel) = true;
seCoop = nanmean((aveCoop));
[sortCoop,fsort] = sort(seCoop);
st = 0;
N = sum(~isnan(sortCoop));
figure(3); clf; plot(sortCoop,'k.'); hold on;
fa = seCoop(fsort);
ng = nearestGene(fsort);
bol_sort = bol_sel(fsort);
plot(find(bol_sort),fa(bol_sort)-.05,'+','color',.034*ones(1,3));
text(2+find(bol_sort),fa(bol_sort)-.05 ,ng(bol_sort),'color','b');
plot(st+1:st+5,fa(st+1:st+5),'+','color','r');  % 
text(st+1:st+5,.92*fa(st+1:st+5),ng(st+1:st+5),'color','r');  % 
plot(N-4:N,fa(N-4:N),'+','color','r');  % 
text(N-4:N,fa(N-4:N),ng(N-4:N),'color','r');  % 
ylabel('average cooperativity');
xlabel('SE ID (sorted)')
ng(st+1:st+10)
ng(N-10:N)

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


figure(2); clf; imagesc(log2(aveCoop)); colorbar;
GetColorMap('RedWhiteBlueSat'); colorbar; clim([-6 6])
cmap = jet(23);
hold on;
for c=1:23
    try
    y=cisCoords{c}(1);
    x = 5;
    w = 10;
    h = cisCoords{c}(end)-cisCoords{c}(1)+1;
    rectangle('Position',[x,y,w,h],'Curvature',1,'FaceColor',cmap(c,:));
    catch
    end
end


%% cis vs trans cooperativity

cisID3 = repmat(cisID,1,1,376);
cisCoop = obsVexp_es(cisID3);
transCoop = obsVexp_es(~cisID3);

figure(1); clf;
subplot(3,1,1);
r = log2(obsVexp_es(:)); r(isnan(r)) = []; r(r==0) = []; 
hist( r , -30:.05:30);
title('obs/exp 3-way SE <600nm');
xlim([-8,8]); xlabel('log2(obs/exp) 3-way cluster');
ylim([0,.5e5])

subplot(3,1,2);
r = cisCoop(:);  r(isnan(r)) = [];
r = log2(r); r(isnan(r)) = []; r(r==0) = []; 
hist( r , -30:.2:30);
title('obs/exp 3-way cis-only');
xlim([-8,8]); xlabel('log2(obs/exp) 3-way cluster');
ylim([0,.5e5])

subplot(3,1,3);
r = transCoop(:);  r(isnan(r)) = [];
r = log2(r); r(isnan(r)) = []; r(r==0) = []; 
hist( r , -30:.05:30);
title('obs/exp ESC 3-way trans-only');
xlim([-8,8]); xlabel('log2(obs/exp) 3-way cluster');
ylim([0,.5e5])


