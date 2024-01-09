
%% Compare ORCA, Hi-C, SPITE and GAM data 


%%  GAM test load
% Note, GAM data is mm9 not mm10.  we will map it to the mm9 SEs. 
GenomeFolder = 'U:\';
NAS02_Vol4 = '\\169.254.113.81\NAS02_Vol4\';
saveFolder =  [NAS02_Vol4,'Derek\SE_analyses_supplementary_files\'];
saveName = [NAS02_Vol4,'Derek\SE_analyses_supplementary_files\chipTracksBinned.mat'];
load([saveFolder,'cisMaps_all.mat'],"fish_maps","sprite_maps","gam_maps","hic_maps");

%%
% we plot the main diagonal as solid (max) for all data, but we exclude it
% from the Correlation analysis next
figure(2); clf;
k=0
for c=1:19
   k=k+1; subplot(19,4,k); a=fish_maps{c};  a(logical(eye(size(a))))=1; skip = isnan(a(2,:)); a(skip,:) = []; a(:,skip)=[];
   a(a<0) =nan;  imagesc(log10(a)); clim([-2.5,-.5]); set(gca,'xtick',[],'ytick',[]); a(logical(eye(size(a))))=nan; fish_maps{c}=a;
   k=k+1; subplot(19,4,k);  a=hic_maps{c};  a(logical(eye(size(a))))=1e3;  a(skip,:) = []; a(:,skip)=[];
   a(a<0) =nan; imagesc(log10(a));  clim([0.1,3]); set(gca,'xtick',[],'ytick',[]);  a(logical(eye(size(a))))=nan; hic_maps{c}=a;
   k=k+1; subplot(19,4,k);  a=sprite_maps{c}; a(logical(eye(size(a))))=1e2;  a(skip,:) = []; a(:,skip)=[];
   a(a<0) =nan; imagesc(log10(a)); clim([-1.1 1.5]); set(gca,'xtick',[],'ytick',[]);  a(logical(eye(size(a))))=nan; sprite_maps{c}=a;
   k=k+1; subplot(19,4,k);  a=gam_maps{c};  a(logical(eye(size(a))))=1; a(skip,:) = []; a(:,skip)=[];
   a(a<0) =nan; imagesc(log10(a));  clim([-1.15 0]); set(gca,'xtick',[],'ytick',[]);  a(logical(eye(size(a))))=nan; gam_maps{c}=a;
end
% 
colormap(flipud(GetColorMap('redToWhiteK')))

%% combine all chromsosomes into a single vector and compute correlation

allSprite = cellfun(@(x) cat(1,x(:)),sprite_maps,'UniformOutput',false);
allSprite = cat(1,allSprite{:});

allGAM = cellfun(@(x) cat(1,x(:)),gam_maps,'UniformOutput',false);
allGAM = cat(1,allGAM{:});

allORCA = cellfun(@(x) cat(1,x(:)),fish_maps,'UniformOutput',false);
allORCA = cat(1,allORCA{:});

allHiC = cellfun(@(x) cat(1,x(:)),hic_maps,'UniformOutput',false);
allHiC = cat(1,allHiC{:});
 
pmat = zeros(2,3);
figure(1); clf; 
subplot(2,3,1); p = PlotCorr(allSprite,allORCA,'log',true,'hex',true); xlabel('SPRITE'); ylabel('ORCA'); pmat(1,1) = p.log10rho;
subplot(2,3,2); p = PlotCorr(allGAM,allORCA,'log',true,'hex',true); xlabel('GAM'); ylabel('ORCA'); pmat(1,2) = p.log10rho;
subplot(2,3,3); p = PlotCorr(allHiC,allORCA,'log',true,'hex',true); xlabel('Hi-C'); ylabel('ORCA'); pmat(1,3) = p.log10rho;
subplot(2,3,4); p = PlotCorr(allHiC,allSprite,'log',true,'hex',true); xlabel('Hi-C');  ylabel('SPRITE');  pmat(2,1) = p.log10rho;
subplot(2,3,5); p = PlotCorr(allHiC,allGAM,'log',true,'hex',true); xlabel('Hi-C');  ylabel('GAM');  pmat(2,2) = p.log10rho;
subplot(2,3,6); p = PlotCorr(allSprite,allGAM,'log',true,'hex',true); ylabel('GAM');  xlabel('SPRITE');  pmat(2,3) = p.log10rho;

figure(3); clf; imagesc(pmat); colorbar; clim([0,1]); hold on;
colormap(flipud(GetColorMap('redToWhiteK')));
text(1,1,'SPRITE v ORCA','HorizontalAlignment','center');
text(1,1,'SPRITE v ORCA','HorizontalAlignment','center');



pmat = zeros(2,3);
figure(1); clf; 
subplot(2,3,1); p = PlotCorr(allSprite,allORCA,'log',true,'hex',true); xlabel('SPRITE'); ylabel('ORCA'); pmat(1,1) = p.log10rho; colorbar;
subplot(2,3,2); p = PlotCorr(allGAM,allORCA,'log',true,'hex',true); xlabel('GAM'); ylabel('ORCA'); pmat(1,2) = p.log10rho;  colorbar;
subplot(2,3,3); p = PlotCorr(allHiC,allORCA,'log',true,'hex',true); xlabel('Hi-C'); ylabel('ORCA'); pmat(1,3) = p.log10rho;  colorbar;
subplot(2,3,4); p = PlotCorr(allHiC,allSprite,'log',true,'hex',true); xlabel('Hi-C');  ylabel('SPRITE');  pmat(2,1) = p.log10rho;  colorbar;
subplot(2,3,5); p = PlotCorr(allHiC,allGAM,'log',true,'hex',true); xlabel('Hi-C');  ylabel('GAM');  pmat(2,2) = p.log10rho;  colorbar;
subplot(2,3,6); p = PlotCorr(allSprite,allGAM,'log',true,'hex',true); ylabel('GAM');  xlabel('SPRITE');  pmat(2,3) = p.log10rho;  colorbar;

%%
if false
f = gcf;
name = ['U:\Manuscripts\SE Clustering Paper\Images\All_chr_orca_hic_sprite_gam.pdf'];
exportgraphics(f,name,'ContentType','vector');

f = gcf;
name = ['U:\Manuscripts\SE Clustering Paper\Images\CorrALL_orca_hic_sprite_gam.pdf'];
exportgraphics(f,name,'ContentType','vector');
end


%%
figure(1); clf;
pmat = zeros(4,4);
k=0;
datas = {allORCA,allGAM,allSprite,allHiC};
names = {'ORCA','GAM','SPRITE','Hi-C'};
for d=1:4
    for e=1:4
        k=k+1;
        figure(2); subplot(4,4,k); 
        p = PlotCorr(datas{e},datas{d},'log',true,'hex',true); xlabel(names{e}); ylabel(names{d});
        pmat(d,e) = p.log10rho;
    end
end

%%
figure(3); clf;
imagesc(pmat);  colorbar; clim([0,1]); hold on;
for d=1:4
    for e=1:4
        % text(d,e,{names{e},names{d},num2str(pmat(d,e),2)},'HorizontalAlignment','center','color','w');
        text(d,e,[num2str(pmat(d,e),2)],'HorizontalAlignment','center','color','w');
    end
end
set(gca,'XTickLabels',names,'YTickLabels',names)
colormap(flipud(GetColorMap('redToWhiteK')));
