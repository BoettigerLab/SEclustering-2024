%% Fig 4 NN analysis


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


%% ===== NN approach
%%
rng default  % For reproducibility
    dTable = seqScores;
    dTable = cat(1,dTable,dTable);
    dTable.out = cpc6;
    dTable(:,3) = []; % start and stop are redundant given the fixed 15kb size
        dTable = chipTable2(:,[1:2,11:end]);
        dTable.speckle = speckleSE;
        dTable.K27ac = chipTable{:,4};
        dTable.lad = lad_score;
        dTable.out = cpc6;
        tbl = dTable;

saveName = [saveFolder,'net6Models_10node_LAD_good.mat'];
if ~exist(saveName,'file')
    
    %% NN (single layer, 10 nodes)  - predict 600 nm contact
    % groups
    %  1 = everything
    %  2 = genome position alone
    %  3 = H3K27ac
    %  4 = OSN
    %  5 = Pol2 Med
    %  6 = speckles + LAD
    
    N= 100;
    netData = cell(N,1);
    netR = nan(N,3);
    netE = nan(N,3);
    for r=1:N
        n = height(tbl);
        tbl.noise = rand(height(tbl),1);   % (a particular random assignment may by chance leave behind a trainable signal (low values on small clusters and high on others) that can buy some correlation, we want to see the typical effect
        % rng default  % For reproducibility
        hpartition = cvpartition(n,'Holdout',0.25); % Nonstratified partition  35 
        idxTrain = training(hpartition);
        trainTable = tbl(idxTrain,:);
        idxTest = test(hpartition);
        testTable = tbl(idxTest,:);
        
      %   figure(2); clf;
        names = {'All','Genome coord.','K27ac','O/S/N','Pol2/Med','Speckle/LAD'};
        selCols ={[1:13], [1:2,13], [10,13], [3:5,13], [6:9,13],  [11:12,13]}; % select columns
        S = length(selCols);
        for s=1:S
            sel = selCols{s};
            netModel = fitcnet(trainTable(:,sel),'out','LayerSizes',10,'Standardize',true);  
            netData{r}{s} = {netModel,trainTable,testTable,idxTest,idxTrain};
            xTest = find(idxTest);
            xTrain = find(idxTrain);
            yTrain = trainTable.out;
            yTrainPred = predict(netModel,trainTable(:,sel));
            yTest = testTable.out;
            yTestPred = predict(netModel,testTable(:,sel));
            rmse = nanmean(abs(yTest-yTestPred));
            rT= PlotCorr(yTest,yTestPred,'showPlot',false);
            netR(r,s) = rT.rho; % pvalue;
            netE(r,s) = rmse;
            % figure(2);  
            % subplot(S,1,s); 
            % plot(xTrain,yTrain,'+','color',.8*ones(1,3)); hold on; plot(xTrain,yTrainPred,'o','color',.8*ones(1,3));
            % plot(xTest,yTest,'b+'); hold on; plot(xTest,yTestPred,'o','color',[.8 .4 0]); 
            % title([names{s},' r=',num2str(rT.rho,2),' p=',num2str(rT.pvalue,2),' RMSE=',num2str(rmse,2)])
            
        end
        disp(['round ',num2str(r)]); pause(.01);
    end
    %%
    save(saveName,'netData','netR','netE');
else
    load(saveName,'netData','netR','netE');
end
%%
figure(4); clf;
violin(netR,'bandwidth',.02,'variableNames',names); ylabel('Pearsons R')

figure(5); clf;
violin(netE,'bandwidth',.02,'variableNames',names); ylabel('RMSE')


figure(4); clf;
boxplot(netR); ylabel('Pearsons R');
set(gca,'xTickLabels',names); box off;

figure(5); clf;
boxplot(netE); ylabel('RMSE');
set(gca,'xTickLabels',names); box off;


%%


figure(4); clf;
violin(netR,'bandwidth',.02,'variableNames',names,'plotMean',false,'faceColor',[.9 .9 .9],'lineColor','none'); ylabel('Pearsons R'); hold on;
boxplot(netR); ylabel('Pearsons R'); hold on;
set(gca,'xTickLabels',names); box off; 

figure(5); clf;

violin(netE,'bandwidth',.02,'variableNames',names,'plotMean',false,'faceColor',[.9 .9 .9],'lineColor','none'); ylabel('RMSE'); hold on;
boxplot(netE); ylabel('RMSE'); hold on;
set(gca,'xTickLabels',names); box off;

med_R = median(netR)
med_RMSE = median(netE)