%% Spontaneous Wake analysis, wakefulness, 06.01.2026
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%% all wake trials 
wakeTrials = contains(WA.recTable.recNames,'Wake') & ...
             cellfun(@isempty,WA.recTable.Remarks);
wakeTable = WA.recTable(wakeTrials,:);

%% subset to analyse for the effect:
rng('shuffle')   % optional: for randomness

animals = unique(wakeTable.Animal);
rowsToKeep = [];

for i = 1:numel(animals)
    idx = find(strcmp(wakeTable.Animal,animals{i}));   % rows of this animal
    nPick = min(3, numel(idx));   % in case some animals have < 3 rows
    rowsToKeep = [rowsToKeep; idx(randperm(numel(idx), nPick))];
end

wakeSubset = wakeTable(rowsToKeep, :);
%%

% wakeSubset = readtable([analysisFolderWake filesep 'wakeSubset.xlsx']);
wakeSubsetAll = readtable([analysisFolderWake filesep 'wakeSubsetAll.xlsx']);
wakeData = wakeSubsetAll(:,["spikes","videoSync","Animal","recNames","Headstage"]);
%% checking correlogram and freq bands:

cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;
maxDendroClusters=3;
normPxxs=cell(height(wakeData),maxDendroClusters);
for i = 1:height(wakeData)
    recName = ['Animal=' wakeData.Animal{i} ',recNames=' wakeData.recNames{i}];
    WA.setCurrentRecording(recName);
    curfreqData = WA.getFreqBandDetection(fMax=80,maxDendroClusters=maxDendroClusters,overwrite=1);
    % WA.plotFreqBandDetection
    curClusterNum=unique(curfreqData.clusters);
    for j = 1:curClusterNum(end)
        normPxxs(i,j)={mean(curfreqData.normsPxx(:,curfreqData.clusters==j),2)};
    end
end
wakeData.normPxxs = normPxxs;
save([analysisFolderWake filesep 'wakeData.mat'],'wakeData')
%% plt unique
for i = 1:height(wakeData)
    figure; plot(x,cell2mat(normPxxs(i,:)))
end
%% plot all
figure; hold on
colors = lines(3);
x = curfreqData.freqHz;
shortinds = x<80;
alphaVal = 0.3;
for c = 1:3
    
    % Convert column to matrix 328 � 27
    mat = cell2mat(normPxxs(:,c)');
    
    plot(x(shortinds), mat(shortinds,:),'Color', [colors(c,:),alphaVal],'LineWidth', 1);
    meanC = mean(mat,2);
    plot(x(shortinds), meanC(shortinds),'Color', colors(c,:),'LineWidth', 2)
end
ylim([0 3])
xlabel('Freq(Hz)'); ylabel('normPxx')
%save
set(gcf,'PaperPositionMode','auto');
fileName=[analysisFolderWake filesep 'clusters_all'];
print(fileName,'-dpdf',['-r300']);

%% ACC to D/B
    hs = WA.recTable.Headstage{WA.currentPRec};
    curSens = cali_result.(hs).sensetivity';
    curZeroG = cali_result.(hs).zeroGbais';
    curLM = WA.getLizardMovements(sensitivity=curSens,zeroGBias=curZeroG);
    DB = WA.getDelta2BetaRatio;


    %%plot ACC to DB
    convWin = 60; %sec
    th= 100;
    dbInd = ~isnan(DB.bufferedDelta2BetaRatio);
    dbtimes= DB.t_ms(dbInd);
    dbratio = DB.bufferedDelta2BetaRatio(dbInd);
    binSize = DB.parDBRatio.movWin -DB.parDBRatio.movOLWin;
    edges = [dbtimes-(binSize/2), dbtimes(end)+(binSize/2)];
    [counts,e] =  histcounts(curLM.t_mov_ms,edges);
    
    h_thresh= counts>th;
    y = ones(1,convWin);
    conv = convn(h_thresh,y,'same');
    f=figure; 
    subplot(2,1,1);plot(conv);hold on;plot(dbratio); legend("mov","DB")
    subplot(2,1,2); scatter(conv,dbratio,'.');
    xlabel('Movement');ylabel("D/B");
    set(f,'PaperPositionMode','auto');
    fileName=[WA.currentPlotFolder filesep 'dbtoACC'];
    print(fileName,'-djpeg',['-r' num2str(WA.figResJPG)]);
    


    % k1 = 2*i-1; k2=2*1;
    % ax1 = nexttile(t,k1);   % dendrogram tile
    % % ax2 = nexttile(t,k2);   % spectra tile
    % WA.plotFreqBandDetection('plotDendrogram',true,'plotSpectralBands',true, ...
    %                      'hDendro',ax1,'hSpectra',ax2, ...
    %                      'savePlots',false);

% end
%% figure out the states parameters:
WA.setCurrentRecording('Animal=PV126,recNames=Wake18');
data = WA.getFreqBandDetection(fMax=80,maxDendroClusters=3);
WA.plotFreqBandDetection;
% take the 50 most type 1 or type 2 from each and plot their pxx
type1i = data.order(1:50);
type2i = data.order(end-49:end);
%for timeseriesviewer:
type1t = data.times(type1i);
type2t = data.times(type2i);
%%
type1pxx = mean(data.sPxx(:,type1i),2);
type2pxx = mean(data.sPxx(:,type2i),2);
figure; plot(data.freqHz,type1pxx);hold on; plot(data.freqHz,type2pxx);legend(["quiet","active"])
title('pxx');xlabel('Freq[Hz]'),ylabel('Power')

type1npxx = mean(data.normsPxx(:,type1i),2);
type2npxx = mean(data.normsPxx(:,type2i),2);
figure; plot(data.freqHz,type1npxx);hold on; plot(data.freqHz,type2npxx);legend(["quiet","active"])
title('normelized pxx');xlabel('Freq[Hz]'),ylabel('nPower')
%% correlation with mov. - acc
% channels = 1:3; startT = 0;win = 1000*60*60*4;
% [V_uV,t_ms] = WA.currentDataObj.getAnalogData(channels,startT,win);
% V_uV = squeeze(V_uV)';


%
cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;
hs = WA.recTable.Headstage{WA.currentPRec};
curSens = cali_result.(hs).sensetivity';
curZeroG = cali_result.(hs).zeroGbais';
LM = WA.getLizardMovements(overwrite=1,sensitivity=curSens,zeroGBias=curZeroG);

%%
bin = 1000;
convWin = 60;
th= 100;
[counts, edges] =  histcounts(LM.t_mov_ms,BinWidth=bin);
h_thresh= counts>th;
% h_times = h.BinEdges;
y = ones(1,convWin);
conv = convn(h_thresh,y,'same');
figure; plot(conv);hold on;plot(DB.bufferedDelta2BetaRatio)
figure; scatter(conv(1:length(DB.bufferedDelta2BetaRatio)),DB.bufferedDelta2BetaRatio,'.');
xlabel('Movement');ylabel("D/B")
% move_conv = conv>0;
% db_stativ = DB.bufferedDelta2BetaRatio;

%% PCA
% running on a single trial, power spec
PxxDataAll=[];
nPxxDataAll=[];
clustersAll=[];
recNamesAll=[];
wakeData.pxxDataAll = cell(height(wakeData),3);
for i = 1:height(wakeData)
    recName = ['Animal=' wakeData.Animal{i} ',recNames=' wakeData.recNames{i}];
    WA.setCurrentRecording(recName);
    curfreqData = WA.getFreqBandDetection;

    % collect the data from the experiment
    pxxData = curfreqData.sPxx';
    PxxDataAll = [PxxDataAll;pxxData]; %save in mat
    % collect the data from the experiment
    npxxData = curfreqData.normsPxx';
    nPxxDataAll = [nPxxDataAll;npxxData]; %save in mat

    % make sure clusters are right
    curNormPxx = curfreqData.normsPxx';     % [nSegments x nFreq] 
    curF = curfreqData.freqHz;   % [nFreq x 1]
    curClust = curfreqData.clusters; % [nSegments x 1]
    % --- delta band mask (0-4 Hz) ---
    deltaMask = (curF >= 0.5) & (curF <= 4);
    % Mean delta power per ORIGINAL cluster label
    origLabels = unique(curClust(~isnan(curClust)));
    clusterDeltaMean = nan(size(origLabels));
    for j = 1:numel(origLabels)
        clusterDeltaMean(j) = mean(curNormPxx(curClust == j,deltaMask),"all", 'omitnan');
    end
    % Rank clusters by delta mean (descending): highest -> lowest
    [~, sortIdx] = sort(clusterDeltaMean, 'descend');
    rankedLabels = origLabels(sortIdx);
    % Build mapping: rankedLabels(1)->1, rankedLabels(2)->2, rankedLabels(3)->3
    newC = curClust;
    for newLab = 1:min(3, numel(rankedLabels))
        newC(curClust == rankedLabels(newLab)) = newLab;
    end 
    
    % add to all clusters
    clustersAll = [clustersAll;newC];

    % add to wakeData/parameter of cell array that can be added to wakeData 
    wakeData.pxxDataAll(i,1) = {pxxData};
    wakeData.pxxDataAll(i,2) = {npxxData};
    wakeData.pxxDataAll(i,3) = {newC};

    % write the exp name. 
    curRecName = repmat({recName},[length(newC),1]);
    recNamesAll = [recNamesAll;curRecName];

    
end
save([analysisFolderWake filesep 'PCAallPxx.mat'],'PxxDataAll','clustersAll','nPxxDataAll')
save([analysisFolderWake filesep 'wakeData.mat'],'wakeData')

%% calc and plot PCA for all same PC
%rmove outliers: 

% Remove extreme total power segments
tp = sum(PxxDataAll,2);
idx = ~isoutlier(tp,'median');
PxxDataAllc = PxxDataAll(idx,:);
clustersAllc = clustersAll(idx);
% % 
% PxxDataAlllog = log10(PxxDataAll);
% PxxDataAlllognorm = PxxDataAlllog ./ sum(PxxDataAlllog,2);
% 
% tp = sum(PxxDataAlllognorm,2);
% idx = ~isoutlier(tp,'median');
% PxxDataAllc = PxxDataAlllognorm(idx,:);

% zscore the data:
PxxDataAllz = zscore(PxxDataAllc);
[coeff, score, latent, explained] = pca(PxxDataAllz);

%%  run the PCA on the normPxx 
tp = sum(nPxxDataAll,2);
idx = ~isoutlier(tp,'median');
nPxxDataAllc = nPxxDataAll(idx,:);
clustersAllc = clustersAll(idx);
[coeff, score, latent, explained,mu] = pca(nPxxDataAllc);


%% add sleep trials?
%get the data from 1 trial 
WA.setCurrentRecording('Animal=PV106,recNames=Night1');
curfreqData = WA.getFreqBandDetection(fMax=80,maxDendroClusters=2,overwrite=1);
nPxx = curfreqData.normsPxx';
clusters=curfreqData.clusters;
% cleaning:
tp = sum(nPxx,2);
idx = ~isoutlier(tp,'median');
nPxxc = nPxx(idx,:);
clustersc = clusters(idx);
% calculate the score acording to the PCA
sleepscore = (nPxxc - mu') * coeff;
sleepclusters=1;
%% scatter plot
labels = clustersAllc; %labels according to the freq clustering
figure;
hold on
colors = lines(length(unique(labels))+2); %adding the sleep colors. 

    for curClust = 1:length(colors)
        idx = labels == curClust;
        scatter(score(idx,1), score(idx,2),20,colors(curClust,:), 'filled','MarkerFaceAlpha',0.5);
        % hist2(score(idx,1), score(idx,2))
    end


hold on
if sleepclusters ==1
    numsleepcluster = unique(clustersc);
    for i = 1:length(numsleepcluster)
        idx = i == clustersc;
        scatter(sleepscore(idx,1), sleepscore(idx,2),20,colors(i+3,:), 'filled','MarkerFaceAlpha',0.5);
        % hist2(score(idx,1), score(idx,2))
    end
end

    xlabel('PC1')
    ylabel('PC2')
    legend('HighD','MiddleD','LowD', 'sleep1','sleep2')
    title('All experiments, all segemnst, normPxx,clean(tp) before PCA')


%% trying to plot in density plots
labels = clustersAllc; %labels according to the freq clustering
colors = lines(length(unique(labels)));
% Inputs:
%   score  : nSegments x nPCs  (use score(:,1) and score(:,2))
%   labels : nSegments x 1     (values 1..3)
pc1 = score(:,1);
pc2 = score(:,2);

K = 3;
cols = lines(K);
alphaMax = 0.7;
nGrid = 150;
qlo = 0.95;      % try 0.90�0.95 for tighter blobs
gamma = 0.25;

% --- common grid ---
pad = 0.05;
xlimAll = [min(pc1) max(pc1)];
ylimAll = [min(pc2) max(pc2)];
xr = diff(xlimAll); yr = diff(ylimAll);

xg = linspace(xlimAll(1)-pad*xr, xlimAll(2)+pad*xr, nGrid);
yg = linspace(ylimAll(1)-pad*yr, ylimAll(2)+pad*yr, nGrid);
[Xg,Yg] = meshgrid(xg,yg);
gridPts = [Xg(:) Yg(:)];

figure; hold on
xlabel('PC1'); ylabel('PC2'); box off

% Optional: faint scatter underlay
% s = scatter(pc1, pc2, 4, 'k', '.'); s.MarkerEdgeAlpha = 0.02;

legH = gobjects(K,1);
legLbl = strings(K,1);

for c = 1:K
    idx = labels == c;
    data = [pc1(idx) pc2(idx)];
    nC = sum(idx);
    legLbl(c) = sprintf('Cluster %d (n=%d)', c, nC);

    % KDE on shared grid
    f = ksdensity(data, gridPts);
    F = reshape(f, size(Xg));

    % threshold background
    thr = quantile(F(:), qlo);
    A = max(F - thr, 0);

    % normalize + shape
    if max(A(:)) > 0
        A = A ./ max(A(:));
    end
    A = (A .^ gamma) * alphaMax;

    % constant-color RGB layer
    RGB = zeros([size(A) 3]);
    RGB(:,:,1) = cols(c,1);
    RGB(:,:,2) = cols(c,2);
    RGB(:,:,3) = cols(c,3);

    % draw blob layer
    h = image(xg, yg, RGB);
    set(h, 'AlphaData', A);

    % legend proxy (dummy patch)
    legH(c) = patch(nan, nan, cols(c,:), ...
        'FaceAlpha', alphaMax, 'EdgeColor', 'none');
end

axis tight
legend(legH, legLbl, 'Location','northeast');
   %% another density plot
pc1 = score(:,1);
pc2 = score(:,2);

K = 3;
cols = lines(K);

nGrid = 220;
pad = 0.05;

nLevels = 9;                 % more = smoother bands
qLevels = linspace(0.70,0.995,nLevels);  % outer -> inner quantiles
alphaOuter = 0.05;
alphaInner = 0.25;           % inner core opacity

% ---- common grid ----
pad = 0.05;
xlimAll = [min(pc1) max(pc1)];
ylimAll = [min(pc2) max(pc2)];
xr = diff(xlimAll); yr = diff(ylimAll);

xg = linspace(xlimAll(1)-pad*xr, xlimAll(2)+pad*xr, nGrid);
yg = linspace(ylimAll(1)-pad*yr, ylimAll(2)+pad*yr, nGrid);
[Xg,Yg] = meshgrid(xg,yg);
gridPts = [Xg(:) Yg(:)];

figure; hold on
xlabel('PC1'); ylabel('PC2'); box off

legH = gobjects(K,1);
legLbl = strings(K,1);

for c = [1,3]
    idx = labels == c;
    data = [pc1(idx) pc2(idx)];
    legLbl(c) = sprintf('Cluster %d (n=%d)', c, sum(idx));

    % KDE on shared grid
    f = ksdensity(data, gridPts);
    F = reshape(f, size(Xg));

    % compute density thresholds from quantiles
    thr = quantile(F(:), qLevels);

    % draw bands from low->high density
    for j = 1:nLevels
        % mask for this level: F >= thr(j)
        M = F >= thr(j);

        % create an RGBA layer using image + AlphaData
        RGB = zeros([size(M) 3]);
        RGB(:,:,1) = cols(c,1);
        RGB(:,:,2) = cols(c,2);
        RGB(:,:,3) = cols(c,3);

        % alpha increases toward inner levels
        a = alphaOuter + (alphaInner-alphaOuter) * (j-1)/(nLevels-1);
        A = a * double(M);

        h = image(xg, yg, RGB);
        set(h, 'AlphaData', A);
    end

    % legend proxy
    legH(c) = patch(nan,nan, cols(c,:), 'FaceAlpha', 0.35, 'EdgeColor','none');
end

axis tight
legend(legH, legLbl, 'Location','northeast');
%% plot all in one grey blob
ind= labels==1 |labels==3;
pc1 = score(ind,1);
pc2 = score(ind,2);

% --- Parameters ---
nGrid = 220;
pad = 0.05;

nLevels = 9;                         % number of density bands
qLevels = linspace(0.75,0.999,nLevels); % outer -> inner quantiles
alphaOuter = 0.03;
alphaInner = 0.45;                    % strong core contrast

grayColor = [0.5 0.5 0.5];            % medium gray

% ---- Common grid ----
xlimAll = [min(pc1) max(pc1)];
ylimAll = [min(pc2) max(pc2)];
xr = diff(xlimAll); 
yr = diff(ylimAll);

xg = linspace(xlimAll(1)-pad*xr, xlimAll(2)+pad*xr, nGrid);
yg = linspace(ylimAll(1)-pad*yr, ylimAll(2)+pad*yr, nGrid);
[Xg,Yg] = meshgrid(xg,yg);
gridPts = [Xg(:) Yg(:)];

% ---- KDE over ALL data ----
data = [pc1 pc2];
f = ksdensity(data, gridPts);
F = reshape(f, size(Xg));

figure; hold on
xlabel('PC1');
ylabel('PC2');
box off

% ---- Draw density bands ----
thr = quantile(F(:), qLevels);

for j = 1:nLevels
    
    mask = F >= thr(j);
    
    RGB = zeros([size(mask) 3]);
    RGB(:,:,1) = grayColor(1);
    RGB(:,:,2) = grayColor(2);
    RGB(:,:,3) = grayColor(3);
    
    a = alphaOuter + (alphaInner-alphaOuter) * (j-1)/(nLevels-1);
    A = a * double(mask);
    
    h = image(xg, yg, RGB);
    set(h, 'AlphaData', A);
end

% axis tight
title(sprintf('All segments (n=%d)', size(score,1)));


%% Plot transient point for the states. (high/low)
crossFreqsAll=[];
S1avg = [];
S3avg = [];
freqdete = WA.getFreqBandDetection;
freqHz = freqdete.freqHz;
for i= 1:height(wakeData)
    curClusters=wakeData.pxxDataAll{i,3};
    curNormsPxx= wakeData.pxxDataAll{i,2}';
    S1= mean(curNormsPxx(:,curClusters==1),2);
    S3= mean(curNormsPxx(:,curClusters==3),2);
    % figure; plot(freqHz,S1); hold on; plot(freqHz,S3)
    if mean(S1(1:3))>mean(S3(1:3))
        crossFreq=freqHz(find(S3-S1>=0,1,'first'));
    else
        crossFreq=freqHz(find(S1-S3>=0,1,'first'));
    end
    crossFreqsAll = [crossFreqsAll;crossFreq];

    % get the avg of values above 1
    curS1Avg = mean(freqHz(S1>1));
    curS3Avg = mean(freqHz(S3>1));
    S1avg = [S1avg;curS1Avg];
    S3avg = [S3avg;curS3Avg];
end
%% plot crossFreq
figure('Color','w');
colors = lines(4);
ax1=subplot(1,3,1:2);
x = ones(length(S1avg),1);
swarmchart(x,S1avg,'filled',color= colors(1,:));
hold on;
swarmchart(x*2,S3avg,'filled',color=repmat(colors(3,:),[length(S1avg),1]));
ylabel("Frequency (Hz)")
xlim([0.5 2.5])
xticks([1,2]); 
xticklabels({'','Cluster 1','','Cluster 2',''})
ax1.XTickLabelRotation=0;
ylims=[0,60];
box(ax1,'off');
ylim(ylims);box(ax1,'off');
title('Clusters Avgs above 1')

ax2=subplot(1,3,3);
swarmchart(ones(size(crossFreqsAll)),crossFreqsAll,'filled','k');
ylim(ylims)
xticks(1);
xticklabels({'',"Cross Freq",''})
xlim([0.5 1.5]);
title('Cross Freq')
box(ax2,'off');
ax2.XTickLabelRotation=0;

% linkaxes([ax1,ax2],'y')


%% plot crossFreq
figure('Color','w');
colors = lines(4);
ax1=subplot(1,3,1:2);
plot([1:2],[S1avg,S3avg],'-','Marker','.','MarkerSize',10)
ylabel("Frequency (Hz)")
xlim([0.5 2.5])
xticks([1,2]); 
xticklabels({'','Cluster 1','','Cluster 2',''})
ax1.XTickLabelRotation=0;
ylims=[0,60];
box(ax1,'off');
ylim(ylims);box(ax1,'off');
title('Clusters Avgs above 1')

ax2=subplot(1,3,3);
swarmchart(ones(size(crossFreqsAll)),crossFreqsAll,'filled','k');
ylim(ylims)
xticks(1);
xticklabels({'',"Cross Freq",''})
xlim([0.5 1.5]);
title('Cross Freq')
box(ax2,'off');
ax2.XTickLabelRotation=0;

% linkaxes([ax1,ax2],'y')
