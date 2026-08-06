%% spontaneouse and hunter analysis for PCA:

WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
% wakeSubsetAll = readtable([analysisFolderWake filesep 'wakeSubsetAll.xlsx']);

%% get wake subset - all first trial   
wakeTrials = contains(WA.recTable.recNames,'Wake') & ...
            ~contains(WA.recTable.Animal,'PV8') & ...
             isnan(WA.recTable.Exclude) & ...
             cellfun(@isempty,WA.recTable.Remarks);
wakeTable = WA.recTable(wakeTrials,:);
[G, AnimalsWake] = findgroups(wakeTable.Animal);
idx = splitapply(@(x){x(1:min(4,numel(x)))}, ...
                 (1:height(wakeTable))', G);

idx = vertcat(idx{:});
% trialsNames = ["Wake1","Wake2","Wake3","Wake4"];
% trialsid  = ismember(wakeTable.recNames,trialsNames);
wakeSubsetAll = wakeTable(idx ,:);

%% get hunter subset
hunTrials = contains(WA.recTable.recNames,'Hunter') & ...
            ~contains(WA.recTable.Animal,'PV8') &...
            isnan(WA.recTable.Exclude);
hunterTable = WA.recTable(hunTrials,:);

[G, AnimalsHun] = findgroups(hunterTable.Animal);
idx = splitapply(@(x){x(1:min(7,numel(x)))}, ...
                 (1:height(hunterTable))', G);
idx = vertcat(idx{:});
hunterSubsetAll = hunterTable(idx ,:);
% hunterSubsetAll = hunterTable(ismember(hunterTable.recNames,trialsNames),:);

%% get sleep subset
sleepTrials = contains(WA.recTable.recNames,'Night') & ...
            ~contains(WA.recTable.Animal,'PV8') &...
            isnan(WA.recTable.Exclude)& ...
             cellfun(@isempty,WA.recTable.Remarks);
sleepTable = WA.recTable(sleepTrials,:);

[G, AnimalsSleep] = findgroups(sleepTable.Animal);
idx = splitapply(@(x){x(1:min(3,numel(x)))}, ...
                 (1:height(sleepTable))', G);
idx = vertcat(idx{:});
sleepSubsetAll = sleepTable(idx ,:);
%% get the data for the pCA - hunter:
% running on a single trial, power spec
PxxDataAllHun = [];
nPxxDataAllHun = [];
timesAllHun   = [];
freqHzAllHun  = [];
recLabelsHun  = [];  % track which recording each row came from
animalLabelHun =[];
clustersHun = [];

for i = 1:height(hunterSubsetAll)
    recName = ['Animal=' hunterSubsetAll.Animal{i} ',recNames=' hunterSubsetAll.recNames{i}];
    animal = hunterSubsetAll.Animal(i);
    WA.setCurrentRecording(recName);
    
    curfreqData = WA.getFreqBandDetection(...
        'fMax',          80, ...
        'segmentLength', 1000, ...
        'binDuration',   5000, ...
        'WelchOL',       0.75, ...
        'overwrite',     1,...
        'dftPoints',     2^10 ...
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
    recLabelsHun  = [recLabelsHun;  repmat(hunterSubsetAll.recNames(i), nBins, 1)];             % which rec each row is from
    animalLabelHun = [animalLabelHun; repmat(animal, nBins, 1)];
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'allSpectraHunter.mat'], 'PxxDataAllHun','nPxxDataAllHun','animalLabelHun', 'clustersHun' ...
    ,'freqHzAllHun', 'timesAllHun', 'recLabelsHun', 'hunterSubsetAll');

%% get the data for the PCA - Wake:
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
        'segmentLength', 1000, ...
        'binDuration',   5000, ...
        'WelchOL',       0.75, ...
        'dftPoints',     2^10, ...
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
    recLabelsWake  = [recLabelsWake;  repmat(wakeSubsetAll.recNames(i), nBins, 1)];             % which rec each row is from
    animalLabelWake = [animalLabelWake; repmat(animal, nBins, 1)];
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'allSpectraWake.mat'], 'PxxDataAllWake','nPxxDataAllWake', 'freqHzAllWake','clustersWake', 'timesAllWake', ...
    'recLabelsWake', 'animalLabelWake');

%% get the pre and post data from the hunters

PxxDataAllBefore = [];
PxxDataAllAfter = [];
% nPxxDataAllHun = [];
% timesAlldots   = [];
freqHzAlldots  = [];
recLabelBefore  = [];  % track which recording each row came from
recLabelAfter = [];
animalLabelBefore =[];
animalLabelAfter =[];
 
PPanimals=AnimalsHun(~strcmp(AnimalsHun,"PV149"));
%
for i = 1:height(hunterSubsetAll)
    animal = hunterSubsetAll.Animal(i);
    if ismember(animal,PPanimals)
        recName = ['Animal=' hunterSubsetAll.Animal{i} ',recNames=' hunterSubsetAll.recNames{i}];
        WA.setCurrentRecording(recName);
        arenaData = WA.getArenaCSVs;
    
        % get the data point for before and after each trial 
        
        segmentLength = 1000;
        binDuration = 5000;
        ch = WA.recTable.defaulLFPCh(WA.currentPRec);
        tStart = 0;
        fMax = 80;
        WelchOL=  0.75;
        dftPoints=2^10;

        WA.getFilters;
        win=WA.currentDataObj.recordingDuration_ms;
        win=floor(win/binDuration)*binDuration; %making win an integer number of segment length
        MLong=WA.currentDataObj.getData(ch,tStart,win);
        % downsample:
        FMLong=WA.filt.F.getFilteredData(MLong);

        timestamps = arenaData.startTrigSh;
  
        segmentSamples = round(segmentLength/1000 * WA.filt.FFs);
        samplesOL      = round(segmentSamples * WelchOL);
        samplesBin     = round(binDuration/1000 * WA.filt.FFs);

        % nBins=numel(FMLong)/samplesBin;

        % FMLongB=reshape(FMLong,[samplesBin,nBins]);

            nTimestamps = numel(timestamps);
            nFreqs = floor(dftPoints/2) + 1;
            pxx_before = zeros(nFreqs, nTimestamps);
            pxx_after  = zeros(nFreqs, nTimestamps);

            for j = 1:nTimestamps
                centerSample = round(timestamps(j)/1000 * WA.filt.FFs);

                idxBefore = (centerSample - samplesBin):(centerSample - 1);
                idxAfter  =  centerSample:(centerSample + samplesBin - 1);

                if idxBefore(1) < 1 || idxAfter(end) > numel(FMLong)
                    warning('Timestamp %d (%.1f ms) too close to edge, skipping.', j, timestamps(j));
                    pxx_before(:,j) = NaN;
                    pxx_after(:,j)  = NaN;
                    continue;
                end

                [pxx_before(:,j), f] = pwelch(squeeze(FMLong(idxBefore)), segmentSamples, samplesOL, dftPoints, WA.filt.FFs);
                [pxx_after(:,j),  ~] = pwelch(squeeze(FMLong(idxAfter)),  segmentSamples, samplesOL, dftPoints, WA.filt.FFs);

            end

            p=find(f<fMax);
            ppBef=find(sum(pxx_before(p,:))<0.4e6); %reject signals with very high amplitudes (probably noise)
            sPxxBef=pxx_before(p,ppBef);
            ppaft=find(sum(pxx_after(p,:))<0.4e6); %reject signals with very high amplitudes (probably noise)
            sPxxAft=pxx_after(p,ppaft);

            freqHz=f(p);
            
            nDotsbef = size(sPxxBef,2);
            nDotaft = size(sPxxAft,2);

        PxxDataAllBefore = [PxxDataAllBefore; sPxxBef'];                         % [totalBins x nFreqs]
        PxxDataAllAfter = [PxxDataAllAfter; sPxxAft'];
        % nPxxDataAllHun = [nPxxDataAllHun; nPxx];
        % clustersHun = [clustersHun; curfreqData.clusters];
        % timesAlldots   = [timesAlldots,   curfreqData.times];               % [1 x totalBins]
        recLabelBefore  = [recLabelBefore;  repmat(hunterSubsetAll.recNames(i), nDotsbef, 1)];             % which rec each row is from
        recLabelAfter = [recLabelAfter ;  repmat(hunterSubsetAll.recNames(i), nDotaft, 1)];             % which rec each row is from
        animalLabelBefore = [animalLabelBefore; repmat(animal, nDotsbef, 1)];
        animalLabelAfter = [animalLabelAfter; repmat(animal, nDotaft, 1)];

    end
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'allSpectraHunterBugs.mat'], 'PxxDataAllBefore','PxxDataAllAfter','animalLabelBefore' ...
    ,'freqHz', 'animalLabelAfter','recLabelBefore','recLabelAfter');

%% get the data for the PCA - SLEEP:
% running on a single trial, power spec
PxxDataAllSleep = [];
nPxxDataAllSleep = [];
timesAllSleep   = [];
freqHzAllSleep  = [];
clustersSleep =[];
clustersSleepAno=[];
recLabelsSleep  = [];  % track which recording each row came from
animalLabelSleep = [];

%
for i = 1:height(sleepSubsetAll)
    recName = ['Animal=' sleepSubsetAll.Animal{i} ',recNames=' sleepSubsetAll.recNames{i}];
    WA.setCurrentRecording(recName);
    animal = sleepSubsetAll.Animal(i);
    % get the sleep timings: 
    WA.getDelta2BetaRatio;
    % WA.plotDelta2BetaRatio;
    curAC = WA.getDelta2BetaAC;
    tStartSleep = curAC.tStartSleep;
    tEndSleep = curAC.tEndSleep;
    win= tEndSleep-tStartSleep;
    
    curfreqData = WA.getFreqBandDetection(...
        'fMax',          80, ...
        'segmentLength', 1000, ...
        'binDuration',   5000, ...
        'WelchOL',       0.75, ...
        'dftPoints',     2^10, ...
        'tStart',       tStartSleep,...
        'win',          win ...
            );
    % 'overwrite', 1 ...
    
    
    pxxData = curfreqData.sPxx';        % [nBins x nFreqs]
    npxxData = curfreqData.normsPxx';
    nBins   = size(pxxData, 1);
    
    % find the SWS and REM according to the freqs.
    meanclus1 = mean(curfreqData.normsPxx(:,curfreqData.clusters==1),2);
    meanclus2 = mean(curfreqData.normsPxx(:,curfreqData.clusters==2),2);
    mainfreq1= mean(curfreqData.freqHz(meanclus1>1));
    mainfreq2= mean(curfreqData.freqHz(meanclus2>1));
    curclustAno = cell(length(curfreqData.clusters),1);
    if mainfreq2>curfreqData.crossFreq && mainfreq1<curfreqData.crossFreq
        %put in clusters that clust1 is SWS and clust2 is REM
        curclustAno(curfreqData.clusters==1) ={"SWS"}; %SWS is 3
        curclustAno(curfreqData.clusters==2) ={"REM"}; %REM is 4
    elseif mainfreq1>curfreqData.crossFreq && mainfreq2<curfreqData.crossFreq
        %put in clusters that clust2 is SWS and clust1 is REM
        curclustAno(curfreqData.clusters==2) ={"SWS"}; %SWS is 3
        curclustAno(curfreqData.clusters==1) ={"REM"}; %REM is 4

    end


    % Verify frequency axis is consistent across recordings
    if isempty(freqHzAllSleep)
        freqHzAllSleep = curfreqData.freqHz; % save once from first recording
    else
        assert(numel(curfreqData.freqHz) == numel(freqHzAllSleep), ...
            'Frequency axis mismatch in rec %s!', recName);
    end

    PxxDataAllSleep = [PxxDataAllSleep; pxxData];
    nPxxDataAllSleep = [nPxxDataAllSleep; npxxData]; % [totalBins x nFreqs]
    clustersSleep = [clustersSleep; curfreqData.clusters];
    clustersSleepAno = [clustersSleepAno;curclustAno];
    timesAllSleep  = [timesAllSleep,   curfreqData.times]; % [1 x totalBins]
    recLabelsSleep  = [recLabelsSleep;  repmat(sleepSubsetAll.recNames(i), nBins, 1)];             % which rec each row is from
    animalLabelSleep = [animalLabelSleep; repmat(animal, nBins, 1)];
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'allSpectraSleep.mat'], 'PxxDataAllSleep','nPxxDataAllSleep', 'freqHzAllSleep','clustersSleep', ...
    'clustersSleepAno','timesAllSleep', 'recLabelsSleep', 'animalLabelSleep');

%% load data
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

load([analysisFolderWake filesep 'allSpectraWake.mat'])
load([analysisFolderWake filesep 'allSpectraHunter.mat'])
load([analysisFolderWake filesep 'allSpectraHunterBugs.mat'])
load([analysisFolderWake filesep 'allSpectraSleep.mat'])
freqHzAll = freqHzAllHun;
% remove 0 freq
f_ind=freqHzAll>0.5;
freqHzAll = freqHzAll(f_ind);
freqHzAllHun = freqHzAllHun(f_ind);
freqHzAllWake= freqHzAllWake(f_ind);
nPxxDataAllWake=nPxxDataAllWake(:,f_ind);
nPxxDataAllHun=nPxxDataAllHun(:,f_ind);
nPxxDataAllSleep = nPxxDataAllSleep(:,f_ind);

PxxDataAllAfter=PxxDataAllAfter(:,f_ind);
PxxDataAllBefore=PxxDataAllBefore(:,f_ind);
PxxDataAllSleep=PxxDataAllSleep(:,f_ind);
PxxDataAllWake=PxxDataAllWake(:,f_ind);
PxxDataAllHun=PxxDataAllHun(:,f_ind);

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

            % Normalize each to [0 1]
            n1 = mat2gray(log10(1+d1_r));  % Hunter -> blue
            n2 = mat2gray(log10(1+d2_r));  % Wake   -> red

            % White background: where both are zero, set to white (1,1,1)
            % Blend: start from white and subtract toward color
            R = 1 - n1;          % blue pulls red channel down
            G = 1 - n1 - n2;     % both pull green down
            B = 1 - n2;          % red pulls blue channel down
            G = max(G, 0);       % clamp to avoid negative values

            rgbCombined = cat(3, R, G, B);

            h_rgb = imagesc(ax, edges1_r, edges2_r, rgbCombined);
            set(ax, 'YDir', 'normal', 'Color', 'white');
            ax.XLim = xlim_r; ax.YLim = ylim_r;
            xlabel('PC1'); ylabel('PC2');
            title([rowLabel ' - Overlap (red=Wake, blue=Hunter)']);
            

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
% animals={'PV162'};
animals = unique(animalLabelWake);
% Parameters for each PCA type: [dx1, dx2, xlim, ylim]
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;  dx2_nt  = 1;
dx1_ns  = 1;  dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-40 120];    pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];    pcaYlim_ns  = [-30 30];

if all(freqHzAllWake==freqHzAllHun)
    freqHzAll = freqHzAllWake;
end

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

            % White background: where both are zero, set to white (1,1,1)
            % Blend: start from white and subtract toward color
            R = 1 - n1;          % blue pulls red channel down
            G = 1 - n1 - n2;     % both pull green down
            B = 1 - n2;          % red pulls blue channel down
            G = max(G, 0);       % clamp to avoid negative values

            rgbCombined = cat(3, R, G, B);

            h_rgb = imagesc(ax, edges1_r, edges2_r, rgbCombined);
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

%% 3D PCA =

% animals={'PV162'};
animals = unique(animalLabelWake);
% Parameters for each PCA type: [dx1, dx2, xlim, ylim]
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;  dx2_nt  = 1;
dx1_ns  = 1;  dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-40 120];    pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];    pcaYlim_ns  = [-30 30];

if all(freqHzAllWake==freqHzAllHun)
    freqHzAll = freqHzAllWake;
end

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

    figure('Name', ['3D PCA - ' thisAnimal]);
    for r = 1:3
        scoresWake_r = pcaData{r,2};
        scoresHun_r  = pcaData{r,3};
        rowLabel     = pcaData{r,5};

        subplot(1,3,r);
        hold on;
        % according to origin:
        scatter3(scoresWake_r(:,1), scoresWake_r(:,2), scoresWake_r(:,3), ...
            20, 'filled');
        scatter3(scoresHun_r(:,1), scoresHun_r(:,2), scoresHun_r(:,3), ...
            20, clusterColors(c,:), 'o');
        % 
        % for c = 1:nClusters
        %     maskC_wake = clustWake == c;
        %     scatter3(scoresWake_r(maskC_wake,1), scoresWake_r(maskC_wake,2), scoresWake_r(maskC_wake,3), ...
        %         20, clusterColors(c,:), 'filled', 'DisplayName', ['Wake C' num2str(c)]);
        % 
        %     maskC_hun = clustHun == c;
        %     scatter3(scoresHun_r(maskC_hun,1), scoresHun_r(maskC_hun,2), scoresHun_r(maskC_hun,3), ...
        %         20, clusterColors(c,:), 'o', 'DisplayName', ['Hun C' num2str(c)]);
        % end
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
        title(rowLabel);
        legend('Location','best','FontSize',6);
        grid on;
        view(45, 30);
        hold off;
    end
    sgtitle(['3D PCA - ' thisAnimal]);
end


%% Project all animals onto one reference animal's Wake PCA

refAnimal = 'PV162'; % the animal whose Wake PCA is the reference

% Parameters for each PCA type: [dx1, dx2, xlim, ylim]
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;    dx2_nt  = 1;
dx1_ns  = 1;    dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-40 120];        pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];         pcaYlim_ns  = [-30 30];

if all(freqHzAllWake == freqHzAllHun)
    freqHzAll = freqHzAllWake;
end

%% Step 1: Compute reference PCA from refAnimal Wake only
maskRef = strcmp(animalLabelWake, refAnimal);

pxxWake_ref   = PxxDataAllWake(maskRef,:);
nPxxWake_ref  = nPxxDataAllWake(maskRef,:);
pxxWake_ref_norm = normalize(pxxWake_ref, 2);

[coeffA_ref,   ~, ~, ~, explainedA_ref]   = pca(pxxWake_ref);
[coeffAnt_ref, ~, ~, ~, explainedAnt_ref] = pca(nPxxWake_ref);
[coeffAns_ref, ~, ~, ~, explainedAns_ref] = pca(pxxWake_ref_norm);

meanA_ref   = mean(pxxWake_ref);
meanAnt_ref = mean(nPxxWake_ref);
meanAns_ref = mean(pxxWake_ref_norm);

%% Step 2: Loop over all animals and project onto reference PCA
animals = unique(animalLabelWake);

for a = 1:numel(animals)
    thisAnimal = animals{a};
    maskHun  = strcmp(animalLabelHun,  thisAnimal);
    maskWake = strcmp(animalLabelWake, thisAnimal);

    % Data for this animal
    pxxHun  = PxxDataAllHun(maskHun,:);
    pxxWake = PxxDataAllWake(maskWake,:);
    nPxxWake = nPxxDataAllWake(maskWake,:);
    nPxxHun  = nPxxDataAllHun(maskHun,:);
    pxxWake_norm = normalize(pxxWake, 2);
    pxxHun_norm  = normalize(pxxHun,  2);

    % Clusters for this animal
    clustHun  = clustersHun(maskHun);
    clustWake = clustersWake(maskWake);
    nClusters = max([clustHun; clustWake]);
    clusterColors = lines(nClusters);

    % --- Project both Wake and Hunter onto REFERENCE PCA space ---
    scoresA_wake   = (pxxWake      - meanA_ref)   * coeffA_ref;
    scoresA_hun    = (pxxHun       - meanA_ref)   * coeffA_ref;

    scoresAnt_wake = (nPxxWake     - meanAnt_ref) * coeffAnt_ref;
    scoresAnt_hun  = (nPxxHun      - meanAnt_ref) * coeffAnt_ref;

    scoresAns_wake = (pxxWake_norm - meanAns_ref) * coeffAns_ref;
    scoresAns_hun  = (pxxHun_norm  - meanAns_ref) * coeffAns_ref;

    pcaData = {
        coeffA_ref,   scoresA_wake,   scoresA_hun,   explainedA_ref,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt_ref, scoresAnt_wake, scoresAnt_hun, explainedAnt_ref, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns_ref, scoresAns_wake, scoresAns_hun, explainedAns_ref, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['Animal: ' thisAnimal ' (ref: ' refAnimal ')'], ...
           'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal ' � projected onto ' refAnimal ' Wake PCA']);

    for r = 1:3
        coeff_r      = pcaData{r,1};
        scoresWake_r = pcaData{r,2};
        scoresHun_r  = pcaData{r,3};
        explained_r  = pcaData{r,4};
        rowLabel     = pcaData{r,5};
        dx1_r        = pcaData{r,6};
        dx2_r        = pcaData{r,7};
        xlim_r       = pcaData{r,8};
        ylim_r       = pcaData{r,9};

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

        n1 = mat2gray(log10(1+d1_r));  % Hunter -> blue
        n2 = mat2gray(log10(1+d2_r));  % Wake   -> red

        R = 1 - n1;
        G = max(1 - n1 - n2, 0);
        B = 1 - n2;
        rgbCombined = cat(3, R, G, B);

        imagesc(ax, edges1_r, edges2_r, rgbCombined);
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
                20, clusterColors(c,:), 'diamond', 'MarkerEdgeAlpha', 0.2, ...
                'LineWidth', 1, 'DisplayName', ['Hun C' num2str(c)]);

            maskC_wake = clustWake == c;
            scatter(scoresWake_r(maskC_wake,1), scoresWake_r(maskC_wake,2), ...
                20, clusterColors(c,:), '*', 'MarkerEdgeAlpha', 0.2, ...
                'DisplayName', ['Wake C' num2str(c)]);
        end
        xlabel('PC1'); ylabel('PC2');
        xlim(xlim_r); ylim(ylim_r);
        legend('Location','best','FontSize',6);
        title([rowLabel ' - Clusters']);
        hold off;
    end
end



%% save all figures
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

%%  Project all animals onto one reference animal's Wake PCA - With Pre+Post!!
%% Project all animals onto one reference animal's Wake PCA

refAnimal = {'PV143','PV162'}; % the animal whose Wake PCA is the reference

% Parameters for each PCA type: [dx1, dx2, xlim, ylim]
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;    dx2_nt  = 1;
dx1_ns  = 1;    dx2_ns  = 1;

pcaXlim_raw = [-2 2]*10000; pcaYlim_raw = [-2 3]*10000;
pcaXlim_nt  = [-40 120];        pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];         pcaYlim_ns  = [-30 30];

% Colors for the extra Before / After datasets
colorBefore = [0   0.6 0  ];   % green
colorAfter  = [1   0.4 0  ];   % orange

if all(freqHzAllWake == freqHzAllHun)
    freqHzAll = freqHzAllWake;
end

%% Step 1: Compute reference PCA from refAnimal Wake only
maskRef = ismember(animalLabelWake, refAnimal);

pxxWake_ref   = PxxDataAllWake(maskRef,:);
nPxxWake_ref  = nPxxDataAllWake(maskRef,:);
pxxWake_ref_norm = normalize(pxxWake_ref, 2);

[coeffA_ref,   ~, ~, ~, explainedA_ref]   = pca(pxxWake_ref);
[coeffAnt_ref, ~, ~, ~, explainedAnt_ref] = pca(nPxxWake_ref);
[coeffAns_ref, ~, ~, ~, explainedAns_ref] = pca(pxxWake_ref_norm);

meanA_ref   = mean(pxxWake_ref);
meanAnt_ref = mean(nPxxWake_ref);
meanAns_ref = mean(pxxWake_ref_norm);

%% Step 2: Loop over all animals and project onto reference PCA
animals = unique(animalLabelWake);
% animals={'PV149'};
for a = 1:numel(animals)
    thisAnimal = animals{a};
    maskHun    = strcmp(animalLabelHun,    thisAnimal);
    maskWake   = strcmp(animalLabelWake,   thisAnimal);
    maskBefore = strcmp(animalLabelBefore, thisAnimal);
    maskAfter  = strcmp(animalLabelAfter,  thisAnimal);

    % Data for this animal
    pxxHun   = PxxDataAllHun(maskHun,:);
    pxxWake  = PxxDataAllWake(maskWake,:);
    nPxxWake = nPxxDataAllWake(maskWake,:);
    nPxxHun  = nPxxDataAllHun(maskHun,:);
    pxxWake_norm = normalize(pxxWake, 2);
    pxxHun_norm  = normalize(pxxHun,  2);

    % New datasets (raw only -> normalize-to-segment computed here)
    pxxBefore = PxxDataAllBefore(maskBefore,:);
    pxxAfter  = PxxDataAllAfter (maskAfter, :);
    pxxBefore_norm = normalize(pxxBefore, 2);
    pxxAfter_norm  = normalize(pxxAfter , 2);

    % Clusters for this animal
    clustHun  = clustersHun(maskHun);
    clustWake = clustersWake(maskWake);
    nClusters = max([clustHun; clustWake]);
    clusterColors = lines(nClusters);

    % --- Project Wake / Hunter onto REFERENCE PCA space ---
    scoresA_wake   = (pxxWake      - meanA_ref)   * coeffA_ref;
    scoresA_hun    = (pxxHun       - meanA_ref)   * coeffA_ref;

    scoresAnt_wake = (nPxxWake     - meanAnt_ref) * coeffAnt_ref;
    scoresAnt_hun  = (nPxxHun      - meanAnt_ref) * coeffAnt_ref;

    scoresAns_wake = (pxxWake_norm - meanAns_ref) * coeffAns_ref;
    scoresAns_hun  = (pxxHun_norm  - meanAns_ref) * coeffAns_ref;

    % --- Project Before / After onto REFERENCE PCA space ---
    % Raw projection (uses raw reference PCA)
    if ~isempty(pxxBefore)
        scoresA_before   = (pxxBefore      - meanA_ref)   * coeffA_ref;
        scoresAns_before = (pxxBefore_norm - meanAns_ref) * coeffAns_ref;
    else
        scoresA_before = []; scoresAns_before = [];
    end
    if ~isempty(pxxAfter)
        scoresA_after    = (pxxAfter       - meanA_ref)   * coeffA_ref;
        scoresAns_after  = (pxxAfter_norm  - meanAns_ref) * coeffAns_ref;
    else
        scoresA_after = []; scoresAns_after = [];
    end

    % Columns 10 & 11 hold Before/After scores; empty for the trial-norm row
    pcaData = {
        coeffA_ref,   scoresA_wake,   scoresA_hun,   explainedA_ref,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw, scoresA_before,   scoresA_after;
        % coeffAnt_ref, scoresAnt_wake, scoresAnt_hun, explainedAnt_ref, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt,  [],               [];
        coeffAns_ref, scoresAns_wake, scoresAns_hun, explainedAns_ref, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns,  scoresAns_before, scoresAns_after;
    };
    nRows = size(pcaData,1);

    figure('Name', ['Animal: ' thisAnimal ' (ref: ' strjoin(refAnimal) ') PP'], ...
           'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal ' -- projected onto ' strjoin(refAnimal) ' Wake PCA, with Pre-Post']);

    for r = 1:nRows
        coeff_r        = pcaData{r,1};
        scoresWake_r   = pcaData{r,2};
        scoresHun_r    = pcaData{r,3};
        explained_r    = pcaData{r,4};
        rowLabel       = pcaData{r,5};
        dx1_r          = pcaData{r,6};
        dx2_r          = pcaData{r,7};
        xlim_r         = pcaData{r,8};
        ylim_r         = pcaData{r,9};
        scoresBefore_r = pcaData{r,10};
        scoresAfter_r  = pcaData{r,11};

        %-- Col 1: PC Loadings
        subplot(nRows, 5, (r-1)*5 + 1);
        hold on;
        for pc = 1:3
            plot(freqHzAll(1:size(coeff_r,1)), coeff_r(:,pc), 'DisplayName', ['PC' num2str(pc)]);
        end
        xlabel('Frequency (Hz)'); ylabel('Loading');
        title([rowLabel ' - Loadings']); legend;

        %-- Col 2: hist2 Hunter (+ Before/After dots)
        subplot(nRows, 5, (r-1)*5 + 2);
        [d1_r, h1_r] = hist2(scoresHun_r(:,1), scoresHun_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Hunter']);
        h1_r.XLim = xlim_r; h1_r.YLim = ylim_r;
        hold(h1_r,'on');
        if ~isempty(scoresBefore_r)
            scatter(h1_r, scoresBefore_r(:,1), scoresBefore_r(:,2), 15, ...
                colorBefore, 'filled',  ...
                'MarkerFaceAlpha',0.7, 'DisplayName','Before');
        end
        if ~isempty(scoresAfter_r)
            scatter(h1_r, scoresAfter_r(:,1), scoresAfter_r(:,2), 15, ...
                colorAfter, 'filled', ...
                'MarkerFaceAlpha',0.7, 'DisplayName','After');
        end
        hold(h1_r,'off');

        %-- Col 3: hist2 Wake
        subplot(nRows, 5, (r-1)*5 + 3);
        [d2_r, h2_r] = hist2(scoresWake_r(:,1), scoresWake_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Wake']);
        h2_r.XLim = xlim_r; h2_r.YLim = ylim_r;

        %-- Col 4: Overlaid hist2 (+ Before/After dots)
        subplot(nRows, 5, (r-1)*5 + 4);
        ax = gca;
        edges1_r = xlim_r(1):dx1_r:(xlim_r(2)+dx1_r);  edges1_r = edges1_r + dx1_r/2;
        edges2_r = ylim_r(1):dx2_r:(ylim_r(2)+dx2_r);  edges2_r = edges2_r + dx2_r/2;

        n1 = mat2gray(log10(1+d1_r));  % Hunter -> blue
        n2 = mat2gray(log10(1+d2_r));  % Wake   -> red

        R = 1 - n1;
        G = max(1 - n1 - n2, 0);
        B = 1 - n2;
        rgbCombined = cat(3, R, G, B);

        imagesc(ax, edges1_r, edges2_r, rgbCombined);
        set(ax, 'YDir', 'normal', 'Color', 'white');
        ax.XLim = xlim_r; ax.YLim = ylim_r;
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Overlap (red=Wake, blue=Hunter)']);
        hold(ax,'on');
        if ~isempty(scoresBefore_r)
            scatter(ax, scoresBefore_r(:,1), scoresBefore_r(:,2), 15, ...
                colorBefore, 'filled',  ...
                'MarkerFaceAlpha',0.8, 'DisplayName','Before');
        end
        if ~isempty(scoresAfter_r)
            scatter(ax, scoresAfter_r(:,1), scoresAfter_r(:,2), 15, ...
                colorAfter, 'filled',  ...
                'MarkerFaceAlpha',0.8, 'DisplayName','After');
        end
        hold(ax,'off');

        %-- Col 5: Scatter colored by cluster (+ Before/After dots)
        subplot(nRows, 5, (r-1)*5 + 5);
        hold on;
        for c = 1:nClusters
            maskC_hun = clustHun == c;
            scatter(scoresHun_r(maskC_hun,1), scoresHun_r(maskC_hun,2), ...
                20, clusterColors(c,:), 'diamond', 'MarkerEdgeAlpha', 0.2, ...
                'LineWidth', 1, 'DisplayName', ['Hun C' num2str(c)]);

            maskC_wake = clustWake == c;
            scatter(scoresWake_r(maskC_wake,1), scoresWake_r(maskC_wake,2), ...
                20, clusterColors(c,:), '*', 'MarkerEdgeAlpha', 0.2, ...
                'DisplayName', ['Wake C' num2str(c)]);
        end
        if ~isempty(scoresBefore_r)
            scatter(scoresBefore_r(:,1), scoresBefore_r(:,2), 15, ...
                colorBefore, 'filled', 'o', ...
                'DisplayName','Before');
        end
        if ~isempty(scoresAfter_r)
            scatter(scoresAfter_r(:,1), scoresAfter_r(:,2), 15, ...
                colorAfter, 'filled', 'square',  ...
                'DisplayName','After');
        end
        xlabel('PC1'); ylabel('PC2');
        xlim(xlim_r); ylim(ylim_r);
        legend('Location','best','FontSize',6);
        title([rowLabel ' - Clusters']);
        hold off;
    end
end


%% save figs
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


%% plot PC3 and 4 
%% Plot Before/After only -- PC3 vs PC4, all animals in one figure

useNormSeg = true;   % true = Norm. to segment PCA, false = Raw PCA

if useNormSeg
    coeff_use = coeffAns_ref;
    mean_use  = meanAns_ref;
    spaceLbl  = 'Norm. to segment';
else
    coeff_use = coeffA_ref;
    mean_use  = meanA_ref;
    spaceLbl  = 'Raw';
end

animals = unique(animalLabelWake);
nA = numel(animals);
nCols = ceil(sqrt(nA));
nRows = ceil(nA / nCols);

figure('Name', ['Before/After PC2 vs PC3 (' spaceLbl ', ref: ' refAnimal ')'], ...
       'Units','normalized', 'Position',[0 0 1 1]);
sgtitle(['Before/After projected onto ' refAnimal ' Wake PCA (' spaceLbl ') -- PC2 vs PC3']);

for a = 1:nA
    thisAnimal = animals{a};
    maskBefore = strcmp(animalLabelBefore, thisAnimal);
    maskAfter  = strcmp(animalLabelAfter,  thisAnimal);

    pxxBefore = PxxDataAllBefore(maskBefore,:);
    pxxAfter  = PxxDataAllAfter (maskAfter, :);

    if useNormSeg
        if ~isempty(pxxBefore), pxxBefore = normalize(pxxBefore, 2); end
        if ~isempty(pxxAfter),  pxxAfter  = normalize(pxxAfter,  2); end
    end

    if ~isempty(pxxBefore)
        scoresBefore = (pxxBefore - mean_use) * coeff_use;
    else
        scoresBefore = [];
    end
    if ~isempty(pxxAfter)
        scoresAfter = (pxxAfter - mean_use) * coeff_use;
    else
        scoresAfter = [];
    end

    subplot(nRows, nCols, a);
    hold on;
    if ~isempty(scoresBefore)
        scatter(scoresBefore(:,2), scoresBefore(:,3), 20, ...
            colorBefore, 'filled', 'MarkerFaceAlpha',0.7, 'DisplayName','Before');
    end
    if ~isempty(scoresAfter)
        scatter(scoresAfter(:,2), scoresAfter(:,3), 20, ...
            colorAfter, 'filled', 'MarkerFaceAlpha',0.7, 'DisplayName','After');
    end
    xlabel('PC2'); ylabel('PC3');
    title(thisAnimal);
    legend('Location','best','FontSize',6);
    grid on;
    hold off;
end



%% PCA with sleep; no pre-post:
%% Project all animals onto one reference animal's Wake PCA

refAnimal = 'PV162';

% Parameters for each PCA type
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;    dx2_nt  = 1;
dx1_ns  = 1;    dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-40 120];        pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];         pcaYlim_ns  = [-30 30];

if all(freqHzAllWake == freqHzAllHun)
    freqHzAll = freqHzAllWake;
end

%% Step 1: Compute reference PCA from refAnimal Wake only
maskRef = strcmp(animalLabelWake, refAnimal);

pxxWake_ref      = PxxDataAllWake(maskRef,:);
nPxxWake_ref     = nPxxDataAllWake(maskRef,:);
pxxWake_ref_norm = normalize(pxxWake_ref, 2);

[coeffA_ref,   ~, ~, ~, explainedA_ref]   = pca(pxxWake_ref);
[coeffAnt_ref, ~, ~, ~, explainedAnt_ref] = pca(nPxxWake_ref);
[coeffAns_ref, ~, ~, ~, explainedAns_ref] = pca(pxxWake_ref_norm);

meanA_ref   = mean(pxxWake_ref);
meanAnt_ref = mean(nPxxWake_ref);
meanAns_ref = mean(pxxWake_ref_norm);

%% Step 2: Loop over all animals and project onto reference PCA
animals = unique(animalLabelWake);
% animals={'PV162'};
for a = 1:numel(animals)
    thisAnimal = animals{a};

    maskHun   = strcmp(animalLabelHun,   thisAnimal);
    maskWake  = strcmp(animalLabelWake,  thisAnimal);
    maskSleep = strcmp(animalLabelSleep, thisAnimal);

    % Data for this animal
    pxxHun   = PxxDataAllHun(maskHun,:);
    pxxWake  = PxxDataAllWake(maskWake,:);
    pxxSleep = PxxDataAllSleep(maskSleep,:);

    nPxxHun   = nPxxDataAllHun(maskHun,:);
    nPxxWake  = nPxxDataAllWake(maskWake,:);
    nPxxSleep = nPxxDataAllSleep(maskSleep,:);

    pxxHun_norm   = normalize(pxxHun,   2);
    pxxWake_norm  = normalize(pxxWake,  2);
    pxxSleep_norm = normalize(pxxSleep, 2);

    % Clusters
    clustHun   = clustersHun(maskHun);
    clustWake  = clustersWake(maskWake);
    clustSleep = clustersSleepAno(maskSleep);

    nClustersWakeHun = max([clustHun; clustWake]);
    % nClustersSleep   = max(clustSleep);
    nClustersSleep =2;
    clustersNames = {"REM","SWS"};

    clusterColorsWakeHun = lines(nClustersWakeHun);
    clusterColorsSleep   = autumn(nClustersSleep);

    % Project Wake, Hunter and Sleep onto reference PCA space
    scoresA_wake   = (pxxWake  - meanA_ref) * coeffA_ref;
    scoresA_hun    = (pxxHun   - meanA_ref) * coeffA_ref;
    scoresA_sleep  = (pxxSleep - meanA_ref) * coeffA_ref;

    scoresAnt_wake  = (nPxxWake  - meanAnt_ref) * coeffAnt_ref;
    scoresAnt_hun   = (nPxxHun   - meanAnt_ref) * coeffAnt_ref;
    scoresAnt_sleep = (nPxxSleep - meanAnt_ref) * coeffAnt_ref;

    scoresAns_wake  = (pxxWake_norm  - meanAns_ref) * coeffAns_ref;
    scoresAns_hun   = (pxxHun_norm   - meanAns_ref) * coeffAns_ref;
    scoresAns_sleep = (pxxSleep_norm - meanAns_ref) * coeffAns_ref;

    pcaData = {
        coeffA_ref,   scoresA_wake,   scoresA_hun,   scoresA_sleep,   explainedA_ref,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt_ref, scoresAnt_wake, scoresAnt_hun, scoresAnt_sleep, explainedAnt_ref, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns_ref, scoresAns_wake, scoresAns_hun, scoresAns_sleep, explainedAns_ref, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['Animal: ' thisAnimal ' (ref: ' refAnimal ')'], ...
           'Units','normalized', 'Position',[0 0 1 1]);

    sgtitle(['Animal: ' thisAnimal ' projected onto ' refAnimal ' Wake PCA']);

    for r = 1:3
        coeff_r       = pcaData{r,1};
        scoresWake_r  = pcaData{r,2};
        scoresHun_r   = pcaData{r,3};
        scoresSleep_r = pcaData{r,4};
        rowLabel      = pcaData{r,6};
        dx1_r         = pcaData{r,7};
        dx2_r         = pcaData{r,8};
        xlim_r        = pcaData{r,9};
        ylim_r        = pcaData{r,10};

        %-- Col 1: PC Loadings
        subplot(3, 6, (r-1)*6 + 1);
        hold on;
        for pc = 1:3
            plot(freqHzAll(1:size(coeff_r,1)), coeff_r(:,pc), ...
                'DisplayName', ['PC' num2str(pc)]);
        end
        xlabel('Frequency (Hz)');
        ylabel('Loading');
        title([rowLabel ' - Loadings']);
        legend;

        %-- Col 2: hist2 Hunter
        subplot(3, 6, (r-1)*6 + 2);
        [dHun_r, hHun_r] = hist2(scoresHun_r(:,1), scoresHun_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r, ...
            'plotProb',1);%,...
            % 'logColorScale',0);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Hunter']);
        hHun_r.XLim = xlim_r; hHun_r.YLim = ylim_r;

        %-- Col 3: hist2 Wake
        subplot(3, 6, (r-1)*6 + 3);
        [dWake_r, hWake_r] = hist2(scoresWake_r(:,1), scoresWake_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r, ...
            'plotProb',1);%,...
            % 'logColorScale',0);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Wake']);
        hWake_r.XLim = xlim_r; hWake_r.YLim = ylim_r;

        %-- Col 4: hist2 Sleep
        subplot(3, 6, (r-1)*6 + 4);
        [dSleep_r, hSleep_r] = hist2(scoresSleep_r(:,1), scoresSleep_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r, ...
            'plotProb',1);%,...
            % 'logColorScale',0);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Sleep']);
        hSleep_r.XLim = xlim_r; hSleep_r.YLim = ylim_r;

        %-- Col 5: Overlaid hist2
        subplot(3, 6, (r-1)*6 + 5);
        ax = gca;

        edges1_r = xlim_r(1):dx1_r:(xlim_r(2)+dx1_r);
        edges2_r = ylim_r(1):dx2_r:(ylim_r(2)+dx2_r);

        edges1_r = edges1_r + dx1_r/2;
        edges2_r = edges2_r + dx2_r/2;

        nHun   = mat2gray(log10(1 + dHun_r));
        nWake  = mat2gray(log10(1 + dWake_r));
        nSleep = mat2gray(log10(1 + dSleep_r));
        % nHun   = mat2gray(dHun_r);
        % nWake  = mat2gray(dWake_r);
        % nSleep = mat2gray(dSleep_r);

        % Wake = red, Hunter = blue, Sleep = green
        R = max(nWake,  0);
        G = max(nSleep, 0);
        B = max(nHun,   0);

        rgbCombined = cat(3, R, G, B);

        imagesc(ax, edges1_r, edges2_r, rgbCombined);
        set(ax, 'YDir', 'normal', 'Color', 'white');
        ax.XLim = xlim_r;
        ax.YLim = ylim_r;

        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Overlap: Wake=red, Hunter=blue, Sleep=green']);

        %-- Col 6: Scatter colored by cluster
        subplot(3, 6, (r-1)*6 + 6);
        hold on;

        for c = 1:nClustersWakeHun
            maskC_hun = clustHun == c;
            scatter(scoresHun_r(maskC_hun,1), scoresHun_r(maskC_hun,2), ...
                10, clusterColorsWakeHun(c,:), 'diamond', ...
                'MarkerEdgeAlpha', 0.2, 'LineWidth', 1, ...
                'DisplayName', ['Hun C' num2str(c)]);

            maskC_wake = clustWake == c;
            scatter(scoresWake_r(maskC_wake,1), scoresWake_r(maskC_wake,2), ...
                10, clusterColorsWakeHun(c,:), '*', ...
                'MarkerEdgeAlpha', 0.2, ...
                'DisplayName', ['Wake C' num2str(c)]);
        end

        for c = 1:nClustersSleep
            % maskC_sleep = clustSleep == c;
            curCluster=clustersNames{c};
            maskC_sleep = string(clustSleep(~cellfun(@isempty,clustSleep)))==curCluster;
            scatter(scoresSleep_r(maskC_sleep,1), scoresSleep_r(maskC_sleep,2), ...
                5, clusterColorsSleep(c,:), 'o', ...
                'MarkerEdgeAlpha', 0.25, ...
                'DisplayName', curCluster);
        end

        xlabel('PC1'); ylabel('PC2');
        xlim(xlim_r); ylim(ylim_r);
        legend('Location','best','FontSize',6);
        title([rowLabel ' - Clusters']);
        hold off;
    end
end


%% PCA with sleep; no pre-post:
%% 1) Each animal's own Wake PCA  2) All-animals Wake PCA

% Parameters for each PCA type
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;    dx2_nt  = 1;
dx1_ns  = 1;    dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-40 120];        pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];         pcaYlim_ns  = [-30 30];

if all(freqHzAllWake == freqHzAllHun)
    freqHzAll = freqHzAllWake;
end

animals = unique(animalLabelWake);

% ============================================================
% APPROACH 2: Compute ONE global PCA from ALL animals' Wake data
% ============================================================
pxxWake_all      = PxxDataAllWake;
nPxxWake_all     = nPxxDataAllWake;
pxxWake_all_norm = normalize(pxxWake_all, 2);

[coeffA_global,   ~, ~, ~, explainedA_global]   = pca(pxxWake_all);
[coeffAnt_global, ~, ~, ~, explainedAnt_global] = pca(nPxxWake_all);
[coeffAns_global, ~, ~, ~, explainedAns_global] = pca(pxxWake_all_norm);

meanA_global   = mean(pxxWake_all);
meanAnt_global = mean(nPxxWake_all);
meanAns_global = mean(pxxWake_all_norm);

% ============================================================
% LOOP over animals � generate two figures per animal
% ============================================================
for a = 1:numel(animals)
    thisAnimal = animals{a};

    maskHun   = strcmp(animalLabelHun,   thisAnimal);
    maskWake  = strcmp(animalLabelWake,  thisAnimal);
    maskSleep = strcmp(animalLabelSleep, thisAnimal);

    % Raw data for this animal
    pxxHun   = PxxDataAllHun(maskHun,:);
    pxxWake  = PxxDataAllWake(maskWake,:);
    pxxSleep = PxxDataAllSleep(maskSleep,:);

    nPxxHun   = nPxxDataAllHun(maskHun,:);
    nPxxWake  = nPxxDataAllWake(maskWake,:);
    nPxxSleep = nPxxDataAllSleep(maskSleep,:);

    pxxHun_norm   = normalize(pxxHun,   2);
    pxxWake_norm  = normalize(pxxWake,  2);
    pxxSleep_norm = normalize(pxxSleep, 2);

    % Clusters
    clustHun   = clustersHun(maskHun);
    clustWake  = clustersWake(maskWake);
    clustSleep = clustersSleepAno(maskSleep);

    nClustersWakeHun = max([clustHun; clustWake]);
    nClustersSleep   = 2;
    clustersNames    = {"REM","SWS"};

    clusterColorsWakeHun = lines(nClustersWakeHun);
    clusterColorsSleep   = autumn(nClustersSleep);

    % ----------------------------------------------------------
    % APPROACH 1: PCA from THIS animal's own Wake data
    % ----------------------------------------------------------
    [coeffA_own,   ~, ~, ~, explainedA_own]   = pca(pxxWake);
    [coeffAnt_own, ~, ~, ~, explainedAnt_own] = pca(nPxxWake);
    [coeffAns_own, ~, ~, ~, explainedAns_own] = pca(pxxWake_norm);

    meanA_own   = mean(pxxWake);
    meanAnt_own = mean(nPxxWake);
    meanAns_own = mean(pxxWake_norm);

    % Project onto own PCA
    scoresA_wake_own   = (pxxWake  - meanA_own) * coeffA_own;
    scoresA_hun_own    = (pxxHun   - meanA_own) * coeffA_own;
    scoresA_sleep_own  = (pxxSleep - meanA_own) * coeffA_own;

    scoresAnt_wake_own  = (nPxxWake  - meanAnt_own) * coeffAnt_own;
    scoresAnt_hun_own   = (nPxxHun   - meanAnt_own) * coeffAnt_own;
    scoresAnt_sleep_own = (nPxxSleep - meanAnt_own) * coeffAnt_own;

    scoresAns_wake_own  = (pxxWake_norm  - meanAns_own) * coeffAns_own;
    scoresAns_hun_own   = (pxxHun_norm   - meanAns_own) * coeffAns_own;
    scoresAns_sleep_own = (pxxSleep_norm - meanAns_own) * coeffAns_own;

    pcaData_own = {
        coeffA_own,   scoresA_wake_own,   scoresA_hun_own,   scoresA_sleep_own,   explainedA_own,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt_own, scoresAnt_wake_own, scoresAnt_hun_own, scoresAnt_sleep_own, explainedAnt_own, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns_own, scoresAns_wake_own, scoresAns_hun_own, scoresAns_sleep_own, explainedAns_own, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['[Own PCA] Animal: ' thisAnimal], ...
           'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal ' � PCA from own Wake data']);
    plotPCAFigure(pcaData_own, freqHzAll, clustHun, clustWake, clustSleep, ...
                  nClustersWakeHun, nClustersSleep, clustersNames, ...
                  clusterColorsWakeHun, clusterColorsSleep);

    % ----------------------------------------------------------
    % APPROACH 2: Project onto global (all-animals) Wake PCA
    % ----------------------------------------------------------
    scoresA_wake_g   = (pxxWake  - meanA_global) * coeffA_global;
    scoresA_hun_g    = (pxxHun   - meanA_global) * coeffA_global;
    scoresA_sleep_g  = (pxxSleep - meanA_global) * coeffA_global;

    scoresAnt_wake_g  = (nPxxWake  - meanAnt_global) * coeffAnt_global;
    scoresAnt_hun_g   = (nPxxHun   - meanAnt_global) * coeffAnt_global;
    scoresAnt_sleep_g = (nPxxSleep - meanAnt_global) * coeffAnt_global;

    scoresAns_wake_g  = (pxxWake_norm  - meanAns_global) * coeffAns_global;
    scoresAns_hun_g   = (pxxHun_norm   - meanAns_global) * coeffAns_global;
    scoresAns_sleep_g = (pxxSleep_norm - meanAns_global) * coeffAns_global;

    pcaData_global = {
        coeffA_global,   scoresA_wake_g,   scoresA_hun_g,   scoresA_sleep_g,   explainedA_global,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt_global, scoresAnt_wake_g, scoresAnt_hun_g, scoresAnt_sleep_g, explainedAnt_global, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns_global, scoresAns_wake_g, scoresAns_hun_g, scoresAns_sleep_g, explainedAns_global, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['[Global PCA] Animal: ' thisAnimal], ...
           'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal ' � projected onto all-animals Wake PCA']);
    plotPCAFigure(pcaData_global, freqHzAll, clustHun, clustWake, clustSleep, ...
                  nClustersWakeHun, nClustersSleep, clustersNames, ...
                  clusterColorsWakeHun, clusterColorsSleep);
end


%% ============================================================
%% LOCAL FUNCTION � shared plotting logic (paste at end of script
%% or save as plotPCAFigure.m)
%% ============================================================
function plotPCAFigure(pcaData, freqHzAll, clustHun, clustWake, clustSleep, ...
                       nClustersWakeHun, nClustersSleep, clustersNames, ...
                       clusterColorsWakeHun, clusterColorsSleep)

    for r = 1:3
        coeff_r       = pcaData{r,1};
        scoresWake_r  = pcaData{r,2};
        scoresHun_r   = pcaData{r,3};
        scoresSleep_r = pcaData{r,4};
        explained_r   = pcaData{r,5};
        rowLabel      = pcaData{r,6};
        dx1_r         = pcaData{r,7};
        dx2_r         = pcaData{r,8};
        xlim_r        = pcaData{r,9};
        ylim_r        = pcaData{r,10};

        %-- Col 1: PC Loadings
        subplot(3, 6, (r-1)*6 + 1);
        hold on;
        for pc = 1:3
            plot(freqHzAll(1:size(coeff_r,1)), coeff_r(:,pc), ...
                'DisplayName', sprintf('PC%d (%.1f%%)', pc, explained_r(pc)));
        end
        xlabel('Frequency (Hz)');
        ylabel('Loading');
        title([rowLabel ' - Loadings']);
        legend;

        %-- Col 2: hist2 Hunter
        subplot(3, 6, (r-1)*6 + 2);
        [dHun_r, hHun_r] = hist2(scoresHun_r(:,1), scoresHun_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r,'plotProb',1);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Hunter']);
        hHun_r.XLim = xlim_r; hHun_r.YLim = ylim_r;

        %-- Col 3: hist2 Wake
        subplot(3, 6, (r-1)*6 + 3);
        [dWake_r, hWake_r] = hist2(scoresWake_r(:,1), scoresWake_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r,'plotProb',1);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Wake']);
        hWake_r.XLim = xlim_r; hWake_r.YLim = ylim_r;

        %-- Col 4: hist2 Sleep
        subplot(3, 6, (r-1)*6 + 4);
        [dSleep_r, hSleep_r] = hist2(scoresSleep_r(:,1), scoresSleep_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r,'plotProb',1);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Sleep']);
        hSleep_r.XLim = xlim_r; hSleep_r.YLim = ylim_r;

        %-- Col 5: Overlaid RGB hist2
        subplot(3, 6, (r-1)*6 + 5);
        ax = gca;

        edges1_r = (xlim_r(1) : dx1_r : (xlim_r(2)+dx1_r)) + dx1_r/2;
        edges2_r = (ylim_r(1) : dx2_r : (ylim_r(2)+dx2_r)) + dx2_r/2;

        nHun   = mat2gray(log10(1 + dHun_r));
        nWake  = mat2gray(log10(1 + dWake_r));
        nSleep = mat2gray(log10(1 + dSleep_r));

        rgbCombined = cat(3, max(nWake,0), max(nSleep,0), max(nHun,0));

        imagesc(ax, edges1_r, edges2_r, rgbCombined);
        set(ax, 'YDir','normal', 'Color','white');
        ax.XLim = xlim_r; ax.YLim = ylim_r;
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Overlap: Wake=red, Hunter=blue, Sleep=green']);

        %-- Col 6: Scatter colored by cluster
        subplot(3, 6, (r-1)*6 + 6);
        hold on;

        for c = 1:nClustersWakeHun
            scatter(scoresHun_r(clustHun==c, 1), scoresHun_r(clustHun==c, 2), ...
                10, clusterColorsWakeHun(c,:), 'diamond', ...
                'MarkerEdgeAlpha',0.2,'LineWidth',1, ...
                'DisplayName', ['Hun C' num2str(c)]);
            scatter(scoresWake_r(clustWake==c, 1), scoresWake_r(clustWake==c, 2), ...
                10, clusterColorsWakeHun(c,:), '*', ...
                'MarkerEdgeAlpha',0.2, ...
                'DisplayName', ['Wake C' num2str(c)]);
        end

        validSleep = ~cellfun(@isempty, clustSleep);
        for c = 1:nClustersSleep
            maskC_sleep = strcmp(string(clustSleep(validSleep)), clustersNames{c});
            scatter(scoresSleep_r(maskC_sleep,1), scoresSleep_r(maskC_sleep,2), ...
                5, clusterColorsSleep(c,:), 'o', ...
                'MarkerEdgeAlpha',0.25, ...
                'DisplayName', clustersNames{c});
        end

        xlabel('PC1'); ylabel('PC2');
        xlim(xlim_r); ylim(ylim_r);
        legend('Location','best','FontSize',6);
        title([rowLabel ' - Clusters']);
        hold off;
    end
end