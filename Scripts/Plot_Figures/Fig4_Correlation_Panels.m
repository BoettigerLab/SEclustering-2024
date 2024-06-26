%% Fig4 Correlation Analysis

% load correlates
saveFolder = 'U:\Manuscripts\SE Clustering Paper\Data\';
chrTableHyb = readtable([saveFolder,'chrTableHyb_mm10_v2.xlsx']) ; % mm10  
seqScores =  readtable('U:\Manuscripts\SE Clustering Paper\Data\ForZenodo\Processed_Seq_Data\SE_TF-ChIPseq_SPRITE-speckle_DamID-LAD_scores.xlsx');
load([saveFolder,'cpc6.mat'],'cpc6');



%% compute pairwise distance in kilobases (linear) between all SEs
dis = chrTableHyb.chrNum*1e9 + chrTableHyb.start;
disM = squareform(pdist([dis,dis]));
disM(disM>1e9) = nan;
figure(4); clf; imagesc(disM); colorbar;

eclose = nansum(disM<1e6,2);  % number of enhancers within 1 Mb
disI = disM; disI(disM==0) = nan;
eWeight = nansum(1./disI,2);  % sum the 1/linear-distance to all other Es. 

figure(4); clf; PlotCorr(eWeight,eclose);
figure(4); clf; PlotCorr(eclose,cpc6(1:376));
xlabel('SEs within 1 Mb');
ylabel('Ave SE cluster size')

figure(4); clf; PlotCorr(eWeight,cpc6(1:376))
xlabel('sum of 1/distance to all SEs');
ylabel('Ave SE cluster size')

cMap = ContactFrac(esMap(1:376,1:376,:),'threshold',.6);
figure(5); clf; PlotCorr(cMap(:),disM(:));
xlabel('contact freq (<600nm)')
ylabel('distance (kb)')


v1 = log10(eWeight);
v2 = cpc6(1:376);
skip = isnan(v1) | isnan(v2);
v1(skip) = []; v2(skip) = [];
figure(4); clf; hexscatter(v1,v2,'res',30);
cp = corr(v1,v2);
title(['r=',num2str(cp,3)]);
colormap(flipud(GetColorMap('magma'))); colorbar;
xlabel('weighted distance 1/\Sigma_i(d_i)');
ylabel('community size (<600nm)');
set(gca,'color','w')

colormap(flipud(GetColorMap('blue'))); colorbar;
set(gca,'color',[.9,.9,.9])
%% combine in table 
dTable1 = seqScores;
dTable1.eWeight = eWeight;
dTable1.eclose = eclose;
dTable1.speckle = speckleSE;




%% plot correlates
pol2 = mean([seqScores.PolII_r1_normalized,seqScores.PolII_r1_normalized],2);
med1 = mean([seqScores.med1_r1_normalized,seqScores.med1_r1_normalized],2);
figure(4); clf;
xVars = {eclose, eWeight, seqScores.NANOG_normalized,seqScores.POU5F1_normalized,seqScores.SOX2_normalized,seqScores.Bingren_mESC_H3K27ac, pol2,med1,seqScores.LAD_score,seqScores.SPRITE_Speckle_score};
xName = {'# SE <1Mb','weighted distance 1/\Sigma_i(d_i)','Nanog','Oct4','Sox2','H3K27ac','PolII','Med1','LAD','Speckle'}

cmap = flipud(gray(10))*.8; % flipud(GetColorMap('redToWhiteK',10));
for x=1:length(xVars) % x=6
    subplot(2,5,x);
    v1 =xVars{x}; %  log10(xVars{x});
    v2 = cpc6(1:376);
    if (max(v1) > min(v1)*10) && min(v1) > 0
        v1 = log10(v1);
    end
    skip = isnan(v1) | isnan(v2) | isinf(v1) | isinf(v2);
    v1(skip) = []; v2(skip) = [];
    hexscatter(v1,v2,'res',40);
    colormap(flipud(GetColorMap('magma'))); colorbar;
    xlabel(xName{x});
    ylabel('community size (<600nm)');
    % set(gca,'color',[.3 .3 .3]);
    clim([0,5]);
    [c1,p1] = corr(v1,v2);
    cf = polyfit(v1,v2,1);
    yy = polyval(cf,v1);
    i = -log10(p1);
    i = round(min(i,10));
    hold on; plot(v1,yy,'-','color',cmap(i,:),'linewidth',2);
    title(['  r=',num2str(c1,2),' p=',num2str(p1,2)]);
    pause(.1);
end

%% 
clrs = [.5 0 1; 
        .5 0 1;
        0 0 .5;
        0 0 .5;
        0 0 .5;
        .8 .5 0;
        0 .5 .8;
        0 .5 .8;
        0 .2 .4;
        0 .5 .8;
        0 .5 .8];

figure(5); clf;
figure(4); clf;
xVars = {eclose, eWeight, seqScores.NANOG_normalized,seqScores.POU5F1_normalized,seqScores.SOX2_normalized,seqScores.Bingren_mESC_H3K27ac, pol2,med1,seqScores.LAD_score,seqScores.SPRITE_Speckle_score};
xName = {'# SE <1Mb','weighted distance 1/\Sigma_i(d_i)','Nanog','Oct4','Sox2','H3K27ac','PolII','Med1','LAD','Speckle'}

rp = zeros(12,2);
for x=1:length(xVars) % x=6
   
    v1 =xVars{x}; %  log10(xVars{x});
    v2 = cpc6(1:376);
    skip = isnan(v1) | isnan(v2) | isinf(v1) | isinf(v2);
    v1(skip) = []; v2(skip) = [];
    figure(4);  subplot(4,3,x);
    cp = PlotCorr(v1,v2,'res',40,'hex',true,'log',true);
    rp(x,1) = cp.log10rho;
    rp(x,2) = cp.log10pvalue;
    colormap(flipud(GetColorMap('magma'))); colorbar;
    xlabel(xName{x});
    ylabel('community size (<600nm)');
    % set(gca,'color',[.3 .3 .3]);
    clim([0,5]);

    figure(5); 
    if rp(x,2)<0.05
        plot(rp(x,1),-log10(rp(x,2)),'.','color',clrs(x,:),'MarkerSize',25); hold on;
        text(rp(x,1)-.4,-log10(rp(x,2)),xName{x},'color',clrs(x,:));
    else
        plot(rp(x,1),-log10(rp(x,2)),'o','color',clrs(x,:),'MarkerSize',5); hold on;
        text(rp(x,1)-.4,-log10(rp(x,2)),xName{x},'color',clrs(x,:));
    end
end
figure(5); xlim([-.7 .7]); ylim([0,40])
ylabel('-log10(p)'); xlabel('Pearsons R')



%% show correlates as a graph in sub mat
figure(1); clf;
 v1 =xVars{7}; %  log10(xVars{x});
 [~,i] = sort(v1);
 seqChns = xVars([1:7,9,10]);
 seqNames = xName([1:7,9,10]);
for x=1:length(seqChns) %
    v1 =seqChns{x}; %  log10(xVars{x});
    figure(1);
    subplot(length(seqChns),1,x);
    semilogy(v1,'.-');
    box off;
    ylabel(seqNames{x});
    clim([0,5]);
end
xlabel('SE ID')
