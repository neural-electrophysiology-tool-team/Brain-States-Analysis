%% spontaneouse and hunter analysis for PCA:

WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
wakeSubsetAll = readtable([analysisFolderWake filesep 'wakeSubsetAll.xlsx']);

%% get hunter subset
hunTrials = contains(WA.recTable.recNames,'Hunter') & ~contains(WA.recTable.Animal,'PV8');
hunterTable = WA.recTable(hunTrials,:);
animals= unique(hunterTable.Animal);
trialsNames = ["Hunter1","Hunter2","Hunter3"];
hunterSubsetAll = hunterTable(ismember(hunterTable.recNames,trialsNames),:);

%% get the data for wake
load([analysisFolderWake filesep 'wakeData.mat'])

%% get the data for the pCA - hunter:
% running on a single trial, power spec
PxxDataAllHun = [];
nPxxDataAllHun = [];
timesAllHun   = [];
freqHzAllHun  = [];
recLabelHun  = [];  % track which recording each row came from
animalLabelHun =[];
clustersHun = [];
%
for i = 1:height(hunterSubsetAll)
    recName = ['Animal=' hunterSubsetAll.Animal{i} ',recNames=' hunterSubsetAll.recNames{i}];
    animal = hunterSubsetAll.Animal(i);
    WA.setCurrentRecording(recName);
    
    curfreqData = WA.getFreqBandDetection(...
        'fMax',          80, ...
        'segmentLength', 2000, ...
        'binDuration',   4000, ...
        'WelchOL',       0.75, ...
        'dftPoints',     2^12 ...
    );
    
    pxxData = curfreqData.sPxx'; % [nBins x nFreqs]
    nPxx =curfreqData.normsPxx';
    nBins   = size(pxxData, 1);

    % Verify frequency axis is consistent across recordings
    if isempty(freqHzAllHun)
        freqHzAllHun = curfreqData.freqHz; % save once from first recording
    else
        assert(numel(curfreqData.freqHz) == numel(freqHzAllHun), ...
            'Frequency axis mismatch in rec %s!', recName);
    end

    PxxDataAllHun = [PxxDataAllHun; pxxData];                         % [totalBins x nFreqs]
    nPxxDataAllHun = [nPxxDataAllHun; nPxx];
    clustersHun = [clustersHun; curfreqData.clusters];
    timesAllHun   = [timesAllHun,   curfreqData.times];               % [1 x totalBins]
    recLabelHun  = [recLabelHun;  repmat(i, nBins, 1)];             % which rec each row is from
    animalLabelHun = [animalLabelHun; repmat(animal, nBins, 1)];
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'allSpectraHunter.mat'], 'PxxDataAllHun','nPxxDataAllHun', 'clustersHun','freqHzAllHun', 'timesAllHun', 'recLabelsHun', 'hunterSubsetAll');

%% get the data for the pCA - Wake:
% running on a single trial, power spec
PxxDataAllWake = [];
nPxxDataAllWake = [];
timesAllWake   = [];
freqHzAllWake  = [];
clustersWake =[];
recLabelsWake  = [];  % track which recording each row came from
animalLabelWake = [];

for i = 1:height(wakeSubsetAll)
    recName = ['Animal=' wakeSubsetAll.Animal{i} ',recNames=' wakeSubsetAll.recNames{i}];
    WA.setCurrentRecording(recName);
    animal = wakeSubsetAll.Animal(i);

    
    curfreqData = WA.getFreqBandDetection(...
        'fMax',          80, ...
        'segmentLength', 2000, ...
        'binDuration',   4000, ...
        'WelchOL',       0.75, ...
        'dftPoints',     2^12, ...
        'overwrite',      1 ...
    );
    
    pxxData = curfreqData.sPxx';        % [nBins x nFreqs]
    npxxData = curfreqData.normsPxx';
    nBins   = size(pxxData, 1);

    % Verify frequency axis is consistent across recordings
    if isempty(freqHzAllWake)
        freqHzAllWake = curfreqData.freqHz; % save once from first recording
    else
        assert(numel(curfreqData.freqHz) == numel(freqHzAllWake), ...
            'Frequency axis mismatch in rec %s!', recName);
    end

    PxxDataAllWake = [PxxDataAllWake; pxxData];
    nPxxDataAllWake = [nPxxDataAllWake; npxxData]; % [totalBins x nFreqs]
    clustersWake = [clustersWake; curfreqData.clusters];
    timesAllWake   = [timesAllWake,   curfreqData.times]; % [1 x totalBins]
    recLabelsWake  = [recLabelsWake;  repmat(i, nBins, 1)];             % which rec each row is from
    animalLabelWake = [animalLabelWake; repmat(animal, nBins, 1)];
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'allSpectraWake.mat'], 'PxxDataAllWake','nPxxDataAllWake', 'freqHzAllWake','clustersWake', 'timesAllWake', ...
    'recLabelsWake', 'animalLabelWake');


%% run the PCA:
PxxDataAll = [PxxDataAllHun ;PxxDataAllWake];
[coeff, scores, ~, ~, explained] = pca(PxxDataAll);
% coeff:   [nFreqs x nComponents]  � frequency loadings
% scores:  [totalBins x nComponents] � each bin in PC space
% recLabels lets you color/separate recordings in the scores plot

%% plot the pca:
% 1. data explaind:
figure;
plot(explained(1:20), '.');
xlabel('PC'); ylabel('% Variance Explained');
title('Scree Plot');

% 2. mean spectrum per rec: Hunter
figure; hold on;
for i = 1:max(recLabelsHun)
    idx = recLabels == i;
    plot(freqHzAll, mean(PxxDataAllHun(idx,:)));
end
xlabel('Frequency (Hz)'); ylabel('Power');
legend(hunterSubsetAll.recNames);
title('Mean spectrum per recording');

% 3. mean spectrum per rec: Wake
figure; hold on;
for i = 1:max(recLabelsWake)
    idx = recLabels == i;
    plot(freqHzAll, mean(PxxDataAllWake(idx,:)));
end
xlabel('Frequency (Hz)'); ylabel('Power');
legend(wakeSubsetAll.recNames);
title('Mean spectrum per recording');

% 4.  PC Loadings � what frequencies drive each PC
figure;
for pc = 1:3
    subplot(3,1,pc);
    plot(freqHzAll(1:size(coeff,2)), coeff(pc,:));
    xlabel('Frequency (Hz)'); ylabel('Loading');
    title(['PC' num2str(pc)]);
    xline(0,'k--');
end

%% 5. plot the pca 1 to 2, color code according to hun/awake
% scatter:
figure;
idx = [ones(size(PxxDataAllHun,1),1);2*ones(size(PxxDataAllWake,1),1)];
s1=scatter(scores(idx==2,1),scores(idx==2,2));
hold on
s2=scatter(scores(idx==1,1),scores(idx==1,2));
hold off


alpha(.3)
% alpha(s2,.3)

xlabel('PC1'); ylabel('PC2'); legend("Hunter","Wake")

% Hist2 in 2 subplots:
figure;
idx = [ones(size(PxxDataAllHun,1),1);2*ones(size(PxxDataAllWake,1),1)];
dx1 = 1; dx2 = 1;
xlim = [-20 160]; ylim = [-70 50];
subplot(1,2,1)
[d1,h1]=hist2(scores(idx==1,1),scores(idx==1,2),'dX1',dx1,'dX2',dx2,'h',gca);
xlabel('PC1'); ylabel('PC2'); title('Hunter')
h1.XLim=xlim; h1.YLim=ylim;
% figure;
subplot(1,2,2)
[d2,h2]=hist2(scores(idx==2,1),scores(idx==2,2),'dX1',dx1,'dX2',dx2,'h',gca);
xlabel('PC1'); ylabel('PC2'); title('Wake');
h2.XLim=xlim; h2.YLim=ylim;


% Hist2 in the same plot
figure;
ax = axes;

% get common edges
edges1 = xlim(1):dx1:xlim(2);
edges2 = ylim(1):dx2:ylim(2);

al=1;
% plot d1 in blue
rgb1 = cat(3, zeros(size(d1)), zeros(size(d1)), ones(size(d1))); % blue
h_d1 = imagesc(ax, edges1, edges2, rgb1);
h_d1.AlphaData = al* mat2gray(log10(1+d1));
hold on;

% plot d2 in red
rgb2 = cat(3, ones(size(d2)), zeros(size(d2)), zeros(size(d2))); % red
h_d2 = imagesc(ax, edges1, edges2, rgb2);
h_d2.AlphaData = al * mat2gray(log10(1+d2));

set(ax, 'YDir', 'normal');
ax.XLim = xlim; ax.YLim = ylim;
xlabel('PC1'); ylabel('PC2');
title('Hunter (blue) vs Wake (red)');


%% Run t-SNE on first 10 PCs rather than raw data
rng(42); % set seed for reproducibility
Y = tsne(scores(:, 1:10));
%% plot T-SNE
%scatter:
idx = [ones(size(PxxDataAllHun,1),1);2*ones(size(PxxDataAllWake,1),1)];

figure;
scatter(Y(idx==1,1), Y(idx==1,2), 20, 'filled');
hold on 
scatter(Y(idx==2,1), Y(idx==2,2), 20, 'filled');
hold off
% colormap(lines(max(idx)));
alpha(.3)
title('t-SNE of LFP spectra');
legend("Hunter","Wake")

% Hist2 in 2 subplots:
figure;
dx1 = 1; dx2 = 1;
xlim = [-80 80]; ylim = [-80 80];

subplot(1,2,1)
[d1,h1]=hist2(Y(idx==1,1),Y(idx==1,2),'dX1',dx1,'dX2',dx2,'h',gca);
title('Hunter')
h1.XLim=xlim; h1.YLim=ylim;
% figure;
subplot(1,2,2)
[d2,h2]=hist2(Y(idx==2,1),Y(idx==2,2),'dX1',dx1,'dX2',dx2,'h',gca);
title('Wake');
h2.XLim=xlim; h2.YLim=ylim;


% Hist2 in the same plot
figure;
ax = axes;
% get common edges
edges1 = xlim(1):dx1:xlim(2);
edges2 = ylim(1):dx2:ylim(2);

al=1;
% plot d1 in blue
rgb1 = cat(3, zeros(size(d1)), zeros(size(d1)), ones(size(d1))); % blue
h_d1 = imagesc(ax, edges1, edges2, rgb1);
h_d1.AlphaData = al* mat2gray(log10(1+d1));
hold on;

% plot d2 in red
rgb2 = cat(3, ones(size(d2)), zeros(size(d2)), zeros(size(d2))); % red
h_d2 = imagesc(ax, edges1, edges2, rgb2);
h_d2.AlphaData = al * mat2gray(log10(1+d2));

set(ax, 'YDir', 'normal');
ax.XLim = xlim; ax.YLim = ylim;
xlabel('PC1'); ylabel('PC2');
title('Hunter (blue) vs Wake (red)');

%% for each animal:

%% Global plots (once)

%% Per-animal loop
dx1 = 500; dx2 = 500;
pcaXlim = [-20 160]; pcaYlim = [-70 50];
animals = unique(animalLabelWake); % assuming same animals in both
% animalsw = unique(animalLabelWake); 
for a = 1:numel(animals)
    thisAnimal = animals{a};

    maskHun  = strcmp(animalLabelHun,  thisAnimal);
    maskWake = strcmp(animalLabelWake, thisAnimal);

    % Data for this animal
    pxxHun  = PxxDataAllHun(maskHun,:);
    pxxWake = PxxDataAllWake(maskWake,:);
    pxxAnimal = [pxxHun; pxxWake];
    nPxxAnimal= [nPxxDataAllHun(maskHun,:);nPxxDataAllWake(maskWake,:)];
    idxAnimal = [ones(size(pxxHun,1),1); 2*ones(size(pxxWake,1),1)];

    % PCA for this animal
    [coeffA, scoresA, ~, ~, explainedA] = pca(pxxAnimal);
    [coeffAnt, scoresAnt, ~, ~, explainedAnt] = pca(nPxxAnimal);
    [coeffAns, scoresAns, ~, ~, explainedAns] = pca(normalize(pxxAnimal,2));

    scoresHun  = scoresA(idxAnimal==1, :);
    scoresWake = scoresA(idxAnimal==2, :);

    figure('Name', ['Animal: ' thisAnimal]);
    sgtitle(['Animal: ' thisAnimal]);

    %-- Panel 1: Scree plot
    subplot(2,4,1);
    plot(explainedA(1:20), '.');
    xlabel('PC'); ylabel('% Variance');
    title('Scree Plot');

    %-- Panel 2: PC Loadings
    subplot(2,4,2);
    hold on;
    for pc = 1:3
        plot(freqHzAll(1:size(coeffA,1)), coeffA(:,pc), 'DisplayName', ['PC' num2str(pc)]);
    end
    xlabel('Frequency (Hz)'); ylabel('Loading');
    title('PC Loadings'); legend;

        %-- Panel 3: Scree plot - normelized
    subplot(2,4,3);
    plot(explainedAns(1:20), '.');
    xlabel('PC'); ylabel('% Variance');
    title('Scree Plot');

    %-- Panel 4: PC Loadings - normelized
    subplot(2,4,4);
    hold on;
    for pc = 1:3
        plot(freqHzAll(1:size(coeffAns,1)), coeffAns(:,pc), 'DisplayName', ['PC' num2str(pc)]);
    end
    xlabel('Frequency (Hz)'); ylabel('Loading');
    title('PC Loadings - normelized'); legend;

    %-- Panel 5: Scatter PC1 vs PC2
    subplot(2,4,5);
    scatter(scoresWake(:,1), scoresWake(:,2), 10, 'r', 'filled'); hold on;
    scatter(scoresHun(:,1),  scoresHun(:,2),  10, 'b', 'filled');
    alpha(0.3);
    xlabel('PC1'); ylabel('PC2');
    legend('Wake','Hunter');
    title('PC1 vs PC2');
    % xlim(pcaXlim); ylim(pcaYlim);

    %-- Panel 6: hist2 Hunter
    subplot(2,4,6);
    [d1,h1] = hist2(scoresHun(:,1), scoresHun(:,2), 'dX1',dx1,'dX2',dx2,'h',gca);
    xlabel('PC1'); ylabel('PC2'); title('Hunter');
    % gca.XLim = pcaXlim; gca.YLim = pcaYlim;

    %-- Panel 7: hist2 Wake
    subplot(2,4,7);
    [d2,h2] = hist2(scoresWake(:,1), scoresWake(:,2), 'dX1',dx1,'dX2',dx2,'h',gca);
    xlabel('PC1'); ylabel('PC2'); title('Wake');
    % gca.XLim = pcaXlim; gca.YLim = pcaYlim;

    %-- Panel 8: Overlaid hist2
    subplot(2,4,8);
    ax = gca;
    edges1 = pcaXlim(1):dx1:pcaXlim(2);
    edges2 = pcaYlim(1):dx2:pcaYlim(2);
    rgb1 = cat(3, zeros(size(d1)), zeros(size(d1)), ones(size(d1)));
    h_d1 = imagesc(ax, edges1, edges2, rgb1);
    h_d1.AlphaData = mat2gray(log10(1+d1));
    hold on;
    rgb2 = cat(3, ones(size(d2)), zeros(size(d2)), zeros(size(d2)));
    h_d2 = imagesc(ax, edges1, edges2, rgb2);
    h_d2.AlphaData = mat2gray(log10(1+d2));
    set(ax, 'YDir', 'normal');
    ax.XLim = pcaXlim; ax.YLim = pcaYlim;
    xlabel('PC1'); ylabel('PC2');
    title('Hunter (blue) vs Wake (red)');
end


%% Per-animal loop - 3x5 figure: rows=PCA type, cols=plot type

% Parameters for each PCA type: [dx1, dx2, xlim, ylim]
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;  dx2_nt  = 1;
dx1_ns  = 1;  dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-20 160];    pcaYlim_nt  = [-70 50];
pcaXlim_ns  = [-20 160];    pcaYlim_ns  = [-70 50];

animals = unique(animalLabelWake);

for a = 1:numel(animals)
    thisAnimal = animals{a};
    maskHun  = strcmp(animalLabelHun,  thisAnimal);
    maskWake = strcmp(animalLabelWake, thisAnimal);

    % Data for this animal
    pxxHun  = PxxDataAllHun(maskHun,:);
    pxxWake = PxxDataAllWake(maskWake,:);
    pxxAnimal = [pxxHun; pxxWake];
    nPxxAnimal = [nPxxDataAllHun(maskHun,:); nPxxDataAllWake(maskWake,:)];
    idxAnimal  = [ones(size(pxxHun,1),1); 2*ones(size(pxxWake,1),1)];

    % --- Run 3 PCAs ---
    [coeffA,   scoresA,   ~, ~, explainedA]   = pca(pxxAnimal);
    [coeffAnt, scoresAnt, ~, ~, explainedAnt] = pca(nPxxAnimal);
    [coeffAns, scoresAns, ~, ~, explainedAns] = pca(normalize(pxxAnimal, 2));

    pcaData = {
        coeffA,   scoresA,   explainedA,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt, scoresAnt, explainedAnt, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns, scoresAns, explainedAns, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['Animal: ' thisAnimal], 'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal]);

    for r = 1:3
        coeff_r    = pcaData{r,1};
        scores_r   = pcaData{r,2};
        explained_r= pcaData{r,3};
        rowLabel   = pcaData{r,4};
        dx1_r      = pcaData{r,5};
        dx2_r      = pcaData{r,6};
        xlim_r     = pcaData{r,7};
        ylim_r     = pcaData{r,8};

        scoresHun_r  = scores_r(idxAnimal==1, :);
        scoresWake_r = scores_r(idxAnimal==2, :);

        %-- Col 1: Scree plot
        subplot(3, 5, (r-1)*5 + 1);
        plot(explained_r(1:20), '.');
        xlabel('PC'); ylabel('% Variance');
        title([rowLabel ' - Scree']);

        %-- Col 2: PC Loadings
        subplot(3, 5, (r-1)*5 + 2);
        hold on;
        for pc = 1:3
            plot(freqHzAll(1:size(coeff_r,1)), coeff_r(:,pc), 'DisplayName', ['PC' num2str(pc)]);
        end
        xlabel('Frequency (Hz)'); ylabel('Loading');
        title([rowLabel ' - Loadings']); legend;


        %-- Col 3: hist2 Hunter
        subplot(3, 5, (r-1)*5 + 3);
        [d1_r, h1_r] = hist2(scoresHun_r(:,1), scoresHun_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Hunter']);
        h1_r.XLim = xlim_r; h1_r.YLim = ylim_r;

        %-- Col 4: hist2 Wake
        subplot(3, 5, (r-1)*5 + 4);
        [d2_r, h2_r] = hist2(scoresWake_r(:,1), scoresWake_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Wake']);
        h2_r.XLim = xlim_r; h2_r.YLim = ylim_r;

        %-- Col 5: Overlaid hist2
        subplot(3, 5, (r-1)*5 + 5);
        ax = gca;

        % Build edges the same way hist2 does internally, so sizes match d1_r and d2_r
        edges1_r = xlim_r(1):dx1_r:(xlim_r(2)+dx1_r);  edges1_r = edges1_r + dx1_r/2;
        edges2_r = ylim_r(1):dx2_r:(ylim_r(2)+dx2_r);  edges2_r = edges2_r + dx2_r/2;

        rgb1 = cat(3, zeros(size(d1_r)), zeros(size(d1_r)), ones(size(d1_r))); % blue = Hunter
        h_d1 = imagesc(ax, edges1_r, edges2_r, rgb1);
        h_d1.AlphaData = mat2gray(log10(1+d1_r));
        set(ax, 'YDir', 'normal');
        hold on;

        rgb2 = cat(3, ones(size(d2_r)), zeros(size(d2_r)), zeros(size(d2_r))); % red = Wake
        h_d2 = imagesc(ax, edges1_r, edges2_r, rgb2);
        h_d2.AlphaData = mat2gray(log10(1+d2_r));

        ax.XLim = xlim_r; ax.YLim = ylim_r;
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Overlap']);
    end
end

%% Per-animal loop - 3x6 figure: rows=PCA type, cols=plot type

% Parameters for each PCA type: [dx1, dx2, xlim, ylim]
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;  dx2_nt  = 1;
dx1_ns  = 1;  dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-40 120];    pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];    pcaYlim_ns  = [-30 30];

animals = unique(animalLabelWake);

for a = 1:numel(animals)
    thisAnimal = animals{a};
    maskHun  = strcmp(animalLabelHun,  thisAnimal);
    maskWake = strcmp(animalLabelWake, thisAnimal);

    % Data for this animal
    pxxHun  = PxxDataAllHun(maskHun,:);
    pxxWake = PxxDataAllWake(maskWake,:);
    pxxAnimal = [pxxHun; pxxWake];
    nPxxAnimal = [nPxxDataAllHun(maskHun,:); nPxxDataAllWake(maskWake,:)];
    idxAnimal  = [ones(size(pxxHun,1),1); 2*ones(size(pxxWake,1),1)];

    % Clusters for this animal
    clustHun  = clustersHun(maskHun);
    clustWake = clustersWake(maskWake);
    nClusters = max([clustHun; clustWake]);
    clusterColors = lines(nClusters);

    % --- Run 3 PCAs ---
    [coeffA,   scoresA,   ~, ~, explainedA]   = pca(pxxAnimal);
    [coeffAnt, scoresAnt, ~, ~, explainedAnt] = pca(nPxxAnimal);
    [coeffAns, scoresAns, ~, ~, explainedAns] = pca(normalize(pxxAnimal, 2));

    pcaData = {
        coeffA,   scoresA,   explainedA,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt, scoresAnt, explainedAnt, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns, scoresAns, explainedAns, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['Animal: ' thisAnimal], 'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal]);

    for r = 1:3
        coeff_r     = pcaData{r,1};
        scores_r    = pcaData{r,2};
        explained_r = pcaData{r,3};
        rowLabel    = pcaData{r,4};
        dx1_r       = pcaData{r,5};
        dx2_r       = pcaData{r,6};
        xlim_r      = pcaData{r,7};
        ylim_r      = pcaData{r,8};

        scoresHun_r  = scores_r(idxAnimal==1, :);
        scoresWake_r = scores_r(idxAnimal==2, :);

        %-- Col 1: PC Loadings
        subplot(3, 5, (r-1)*5 + 1);
        hold on;
        for pc = 1:3
            plot(freqHzAll(1:size(coeff_r,1)), coeff_r(:,pc), 'DisplayName', ['PC' num2str(pc)]);
        end
        xlabel('Frequency (Hz)'); ylabel('Loading');
        title([rowLabel ' - Loadings']); legend;

        %-- Col 2: hist2 Hunter
        subplot(3, 5, (r-1)*5 + 2);
        [d1_r, h1_r] = hist2(scoresHun_r(:,1), scoresHun_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Hunter']);
        h1_r.XLim = xlim_r; h1_r.YLim = ylim_r;

        %-- Col 3: hist2 Wake
        subplot(3, 5, (r-1)*5 + 3);
        [d2_r, h2_r] = hist2(scoresWake_r(:,1), scoresWake_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Wake']);
        h2_r.XLim = xlim_r; h2_r.YLim = ylim_r;

        %-- Col 4: Overlaid hist2
        subplot(3, 5, (r-1)*5 + 4);
        ax = gca;
        edges1_r = xlim_r(1):dx1_r:(xlim_r(2)+dx1_r);  edges1_r = edges1_r + dx1_r/2;
        edges2_r = ylim_r(1):dx2_r:(ylim_r(2)+dx2_r);  edges2_r = edges2_r + dx2_r/2;

        rgb1 = cat(3, zeros(size(d1_r)), zeros(size(d1_r)), ones(size(d1_r)));
        h_d1 = imagesc(ax, edges1_r, edges2_r, rgb1);
        h_d1.AlphaData = mat2gray(log10(1+d1_r));
        set(ax, 'YDir', 'normal');
        hold on;
        rgb2 = cat(3, ones(size(d2_r)), zeros(size(d2_r)), zeros(size(d2_r)));
        h_d2 = imagesc(ax, edges1_r, edges2_r, rgb2);
        h_d2.AlphaData = mat2gray(log10(1+d2_r));
        ax.XLim = xlim_r; ax.YLim = ylim_r;
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Overlap']);

        %-- Col 5: Scatter colored by cluster
        subplot(3, 5, (r-1)*5 + 5);
        hold on;
        for c = 1:nClusters
            maskC_wake = clustWake == c;
            scatter(scoresWake_r(maskC_wake,1), scoresWake_r(maskC_wake,2), ...
                20, clusterColors(c,:), '*', ...
                'DisplayName', ['Wake C' num2str(c)]);
            maskC_hun = clustHun == c;
            scatter(scoresHun_r(maskC_hun,1), scoresHun_r(maskC_hun,2), ...
                20, clusterColors(c,:),'.', ...
                'LineWidth', 1, ...
                'DisplayName', ['Hun C' num2str(c)]);
        end
        alpha(0.4);
        xlabel('PC1'); ylabel('PC2');
        xlim(xlim_r); ylim(ylim_r);
        legend('Location','best','FontSize',6);
        title([rowLabel ' - Clusters']);
        hold off;
    end
end



%%
saveDir = [analysisFolderWake filesep 'PCAeachAnimal'];
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% Get all open figures
figHandles = findall(0, 'Type', 'figure');

% Save each one
for f = 1:numel(figHandles)
    fig = figHandles(f);
    figName = fig.Name;
    
    % Replace spaces/colons that are invalid in filenames
    figName = strrep(figName, ' ', '_');
    figName = strrep(figName, ':', '');
    
    % Save as both .fig and .png
    % savefig(fig, fullfile(saveDir, [figName '.fig']));
    exportgraphics(fig, fullfile(saveDir, [figName '.png']), 'Resolution', 300);
end

%% run on PCA and plot hunter on top of it.
animals={'PV162'};
% animals = unique(animalLabelWake);


for a = 1:numel(animals)
    thisAnimal = animals{a};
    maskHun  = strcmp(animalLabelHun,  thisAnimal);
    maskWake = strcmp(animalLabelWake, thisAnimal);

    % Data for this animal
    pxxHun  = PxxDataAllHun(maskHun,:);
    pxxWake = PxxDataAllWake(maskWake,:);
    % pxxAnimal = [pxxHun; pxxWake];
    nPxxWake = nPxxDataAllWake(maskWake,:);
    nPxxHun = nPxxDataAllHun(maskHun,:);
    % nPxxAnimal = [nPxxWake; nPxxHun];
    idxAnimal  = [ones(size(pxxHun,1),1); 2*ones(size(pxxWake,1),1)];

    % Clusters for this animal
    clustHun  = clustersHun(maskHun);
    clustWake = clustersWake(maskWake);
    nClusters = max([clustHun; clustWake]);
    clusterColors = lines(nClusters);

    % --- Run 3 PCAs on Wake only ---
    [coeffA,   scoresA_wake,   ~, ~, explainedA]   = pca(pxxWake);
    [coeffAnt, scoresAnt_wake, ~, ~, explainedAnt] = pca(nPxxWake);
    [coeffAns, scoresAns_wake, ~, ~, explainedAns] = pca(normalize(pxxWake, 2));

    % --- Project Hunter data onto Wake PCA space ---
    % For raw:
    scoresA_hun   = (pxxHun - mean(pxxWake)) * coeffA;

    % For norm to trial:
    scoresAnt_hun = (nPxxHun - mean(nPxxWake)) * coeffAnt;

    % For norm to segment:
    pxxWake_norm  = normalize(pxxWake, 2);
    pxxHun_norm   = normalize(pxxHun,  2);
    scoresAns_hun = (pxxHun_norm - mean(pxxWake_norm)) * coeffAns;

    pcaData = {
        coeffA,   scoresA_wake,   scoresA_hun,   explainedA,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt, scoresAnt_wake, scoresAnt_hun, explainedAnt, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns, scoresAns_wake, scoresAns_hun, explainedAns, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
        };

    figure('Name', ['Animal: ' thisAnimal], 'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal]);


        for r = 1:3
            coeff_r       = pcaData{r,1};
            scoresWake_r  = pcaData{r,2};   % from pca()
            scoresHun_r   = pcaData{r,3};   % projected manually
            explained_r   = pcaData{r,4};
            rowLabel      = pcaData{r,5};
            dx1_r         = pcaData{r,6};
            dx2_r         = pcaData{r,7};
            xlim_r        = pcaData{r,8};
            ylim_r        = pcaData{r,9};

            % scoresHun_r  = scores_r(idxAnimal==1, :);
            % scoresWake_r = scores_r(idxAnimal==2, :);

            %-- Col 1: PC Loadings
            subplot(3, 5, (r-1)*5 + 1);
            hold on;
            for pc = 1:3
                plot(freqHzAll(1:size(coeff_r,1)), coeff_r(:,pc), 'DisplayName', ['PC' num2str(pc)]);
            end
            xlabel('Frequency (Hz)'); ylabel('Loading');
            title([rowLabel ' - Loadings']); legend;

            %-- Col 2: hist2 Hunter
            subplot(3, 5, (r-1)*5 + 2);
            [d1_r, h1_r] = hist2(scoresHun_r(:,1), scoresHun_r(:,2), ...
                'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
                'limits1',xlim_r,'limits2',ylim_r);
            xlabel('PC1'); ylabel('PC2');
            title([rowLabel ' - Hunter']);
            h1_r.XLim = xlim_r; h1_r.YLim = ylim_r;

            %-- Col 3: hist2 Wake
            subplot(3, 5, (r-1)*5 + 3);
            [d2_r, h2_r] = hist2(scoresWake_r(:,1), scoresWake_r(:,2), ...
                'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
                'limits1',xlim_r,'limits2',ylim_r);
            xlabel('PC1'); ylabel('PC2');
            title([rowLabel ' - Wake']);
    
            %-- Col 4: Overlaid hist2
            subplot(3, 5, (r-1)*5 + 4);
            ax = gca;
            edges1_r = xlim_r(1):dx1_r:(xlim_r(2)+dx1_r);  edges1_r = edges1_r + dx1_r/2;
            edges2_r = ylim_r(1):dx2_r:(ylim_r(2)+dx2_r);  edges2_r = edges2_r + dx2_r/2;

            % Normalize each to [0 1]
            n1 = mat2gray(log10(1+d1_r));  % Hunter -> blue
            n2 = mat2gray(log10(1+d2_r));  % Wake   -> red

            % Build single RGB image: R=Wake, G=0, B=Hunter
            rgbCombined = cat(3, n2, zeros(size(n1)), n1);

            % Alpha = max of both � empty bins fully transparent, occupied bins opaque
            alphaMap = max(n1, n2);

            h_rgb = imagesc(ax, edges1_r, edges2_r, rgbCombined);
            h_rgb.AlphaData = alphaMap;
            set(ax, 'YDir', 'normal', 'Color', 'white');
            ax.XLim = xlim_r; ax.YLim = ylim_r;
            xlabel('PC1'); ylabel('PC2');
            title([rowLabel ' - Overlap (red=Wake, blue=Hunter)']);
            

            %-- Col 5: Scatter colored by cluster
            subplot(3, 5, (r-1)*5 + 5);
            hold on;
            for c = 1:nClusters
                maskC_hun = clustHun == c;
                scatter(scoresHun_r(maskC_hun,1), scoresHun_r(maskC_hun,2), ...
                    20, clusterColors(c,:), 'diamond','MarkerEdgeAlpha',0.2, ...
                    'LineWidth', 1, ...
                    'DisplayName', ['Hun C' num2str(c)]);
                
                maskC_wake = clustWake == c;
                scatter(scoresWake_r(maskC_wake,1), scoresWake_r(maskC_wake,2), ...
                    20, clusterColors(c,:), '*','MarkerEdgeAlpha',0.2, ...
                    'DisplayName', ['Wake C' num2str(c)]);


            end
            xlabel('PC1'); ylabel('PC2');
            xlim(xlim_r); ylim(ylim_r);
            legend('Location','best','FontSize',6);
            title([rowLabel ' - Clusters']);
            hold off;
        end

end