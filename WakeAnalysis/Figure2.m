%% Figure 2 - Wake vs Sleep figrue
%% get the data
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%% Panel A - examples of sleep and wake

%% Panel B - correlogram both sleep and wake
% good example:
% wake - PV157,Wake3,tstart 0.5
% sleep - PV157, Night1, tstart 7 (or 5.5)
% second example:
% Wake - PV159,Wake4, tstart 0,
% sleep  - PV159, Night2, tStart 5.5

% data parameters:
fMax = 80;
segmentLength = 1000;
binDuration = 5000;
WelchOL = 0.75;
dftPoints = 2^10;
win = 0.5*60*60*1000;

% get data from Wake
WA.setCurrentRecording('Animal=PV157,recNames=Wake3');
wakeFreqData = WA.getFreqBandDetection(...
    'fMax',          fMax, ...
    'segmentLength', segmentLength, ...
    'binDuration',   binDuration, ...
    'WelchOL',       WelchOL, ...
    'dftPoints',     dftPoints, ...
    'overwrite',     1,...
    'tStart',        0.5*60*60*1000, ...
    'win', win ...
    );
WA.plotFreqBandDetection;

% get the data from sleep

WA.setCurrentRecording('Animal=PV157,recNames=Night1');
sleepFreqData = WA.getFreqBandDetection(...
    'fMax',          fMax, ...
    'segmentLength', segmentLength, ...
    'binDuration',   binDuration, ...
    'WelchOL',       WelchOL, ...
    'dftPoints',     dftPoints, ...
    'overwrite',     1, ...
    'tStart',        5.5*60*60*1000, ...
    'win', win ...
    );
    WA.plotFreqBandDetection;

%% colors - same palette as the bout-duration violin figures
% row 1 = LowDB, row 2 = HighDB. Which dendrogram cluster is HighDB is decided
% per state from the mean delta/beta ratio of its bins (see drawStateTree below).
clusterColorsWake  = [0.30 0.62 0.85;    % Wake  LowDB  - light blue (was 0.40 0.70 0.90)
                      0.93 0.60 0.40];   % Wake  HighDB - light orange
clusterColorsSleep = [0.00 0.45 0.74;    % Sleep LowDB  - dark blue
                      0.85 0.33 0.10];   % Sleep HighDB - dark orange

% bands for the delta/beta ranking of the clusters. Only the ranking matters,
% but keep them in line with WA.getDelta2BetaRatio.
dbBands.delta = [0.5 4];    % Hz
dbBands.beta  = [10 30];    % Hz

%% plot them together with the dendro tree:
figW = 798;  figH = 650;                 % extra height at the bottom for the horizontal colour bar
f = figure('Position',[100 100 figW figH]);
margin = 0.05;
leftWidth  = 0.18;
cbarSpace   = 0.08;                      % vertical room below the matrix for the colour bar and its tick labels
cbarThick   = 0.02;                      % colour bar thickness (fraction of figure height)
cbarLenFrac = 0.25;                      % colour bar length as a fraction of the correlogram width
usableHeight = 1 - 2*margin - cbarSpace;
axBottom = margin + cbarSpace;           % bottom edge of the trees and the matrix
ax3Width = usableHeight * figH / figW;   % matrix axes exactly square in pixels, so the trees stay aligned with the rows
ax1 = axes('Parent', f, 'Position', [margin, axBottom + usableHeight/2, leftWidth, usableHeight/2]);
ax2 = axes('Parent', f, 'Position', [margin, axBottom,                  leftWidth, usableHeight/2]);
ax3 = axes('Parent', f, 'Position', [margin + leftWidth, axBottom, ax3Width, usableHeight]);
set(ax1, 'PositionConstraint', 'innerposition');
set(ax2, 'PositionConstraint', 'innerposition');
set(ax3, 'PositionConstraint', 'innerposition');

maxDendroClusters = 2;

% wake tree into ax1 - cluster colors: HighDB orange, LowDB blue (light shades)
[orderW, clustersW, cMapW] = drawStateTree(wakeFreqData, ax1, ax3, f, ...
    maxDendroClusters, clusterColorsWake, dbBands, 'Wake');
ylabel(ax1, 'Wake')

% sleep tree into ax2 - same, dark shades
[orderS, clusterS, cMapS] = drawStateTree(sleepFreqData, ax2, ax3, f, ...
    maxDendroClusters, clusterColorsSleep, dbBands, 'Sleep');
ylabel(ax2, 'Sleep')

% no numbers on the tree axes either (delete this line to keep the distance scale)
set([ax1 ax2], 'XTick', [], 'YTick', []);

% Combined ordered matrix into ax3
normsPXXboth = [wakeFreqData.normsPxx, sleepFreqData.normsPxx];
corrMatboth = corrcoef(normsPXXboth);
nWake = length(orderW);
nSleep = length(orderS);
orderBoth = [orderW, orderS + nWake];
reordered_corr = corrMatboth(orderBoth, orderBoth);

delete(ax3)
ax3 = axes('Parent', f, 'Position', [margin + leftWidth, axBottom, ax3Width, usableHeight]);
set(ax3, 'PositionConstraint', 'innerposition');
axes(ax3)
imagesc(reordered_corr);
axis square;
colormap(parula);
clim([-1 1]);
set(ax3, 'XTick', [], 'YTick', []);      % no bin numbers on the left / bottom

% horizontal colour bar below the matrix: right-aligned with it, 1/4 of its width
cb = colorbar(ax3, 'Location', 'southoutside');
cbarLen = cbarLenFrac * ax3Width;
cb.Position = [margin + leftWidth + ax3Width - cbarLen, margin + 0.03, cbarLen, cbarThick];
cb.Ticks    = [-1 0 1];
cb.FontSize = 8;

hold(ax3, 'on');

% --- 1) black divider lines between wake and sleep blocks ---
xline(ax3, nWake+0.5, 'k-', 'LineWidth', 1);
yline(ax3, nWake+0.5, 'k-', 'LineWidth', 1);

% --- 2) "Wake" / "Sleep" titles above ax3, split proportionally ---
nTotal = nWake + nSleep;
wakeFrac  = nWake  / nTotal;
sleepFrac = nSleep / nTotal;
text(ax3, wakeFrac/2, 1.03, 'Wake', 'Units','normalized', ...
    'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',12);
text(ax3, wakeFrac + sleepFrac/2, 1.03, 'Sleep', 'Units','normalized', ...
    'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize',12);

% --- 2b) the same labels on the right side, rotated (wake block is the top rows: YDir is reversed) ---
text(ax3, 1.01, 1 - wakeFrac/2, 'Wake', 'Units','normalized', 'Rotation', 90, ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontWeight','bold', 'FontSize',12);
text(ax3, 1.01, sleepFrac/2, 'Sleep', 'Units','normalized', 'Rotation', 90, ...
    'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontWeight','bold', 'FontSize',12);

% --- 3) redraw cluster-state boxes from wake/sleep trees onto ax3, with the state colors ---
[pStartW, pEndW, cluW] = getClusterBoundaries(clustersW, orderW);
[pStartS, pEndS, cluS] = getClusterBoundaries(clusterS, orderS);
pStartS = pStartS + nWake;   % shift sleep block into combined coordinates
pEndS   = pEndS + nWake;

for i = 1:numel(pStartW)
    rectangle(ax3, 'Position',[pStartW(i) pStartW(i) pEndW(i)-pStartW(i) pEndW(i)-pStartW(i)], ...
        'EdgeColor', cMapW(cluW(i),:), 'LineWidth', 3);
end
for i = 1:numel(pStartS)
    rectangle(ax3, 'Position',[pStartS(i) pStartS(i) pEndS(i)-pStartS(i) pEndS(i)-pStartS(i)], ...
        'EdgeColor', cMapS(cluS(i),:), 'LineWidth', 3);
end

hold(ax3, 'off');

% save figure:
% set(f,'PaperPosition',[1 5 4 4]);
fileName = [analysisFolderWake filesep 'SleepWakeCorrelogram'];
print(fileName,'-dpdf','-r300');

%% Panel C - mean normalized PSD of each cluster, same colors as the correlogram
% needs clustersW / cMapW and clusterS / cMapS from the section above
showSEM = true;                          % shade mean +/- SEM across the bins of each cluster

fPSD  = figure('Position',[100 100 420 320]);
axPSD = axes('Parent', fPSD);
hold(axPSD, 'on');
plotClusterMeanPSD(axPSD, wakeFreqData,  clustersW, cMapW, clusterColorsWake,  'Wake',  showSEM);
plotClusterMeanPSD(axPSD, sleepFreqData, clusterS,  cMapS, clusterColorsSleep, 'Sleep', showSEM);
hold(axPSD, 'off');
xlim(axPSD, [0 fMax]);
xlabel(axPSD, 'Frequency (Hz)');
ylabel(axPSD, 'nPSD');
box(axPSD, 'off');
legend(axPSD, 'Location', 'best');

% save figure:
set(fPSD, 'PaperPosition', [1 5 4 3]);
fileName = [analysisFolderWake filesep 'ClusterMeanNormPSD'];
print(fPSD, fileName, '-dpdf', '-r300');

% --- helper: mean (+/- SEM) normsPxx of every cluster of one state, in the cluster's color ---
function h = plotClusterMeanPSD(ax, freqData, clusters, cMap, stateColors, stateName, showSEM)
    fHz  = freqData.freqHz(:)';
    npxx = freqData.normsPxx;                              % expected [nFreq x nBins]
    if size(npxx,1) ~= numel(fHz) && size(npxx,2) == numel(fHz)
        npxx = npxx';
    end
    nClu = max(clusters);
    h = gobjects(nClu, 1);
    for c = 1:nClu
        X = npxx(:, clusters == c);
        m = mean(X, 2)';
        if isequal(cMap(c,:), stateColors(2,:))            % color map row 2 = HighDB (see drawStateTree)
            cluName = 'HighDB';
        else
            cluName = 'LowDB';
        end
        if showSEM
            s = std(X, 0, 2)' / sqrt(size(X, 2));
            fill(ax, [fHz, fliplr(fHz)], [m + s, fliplr(m - s)], cMap(c,:), ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        h(c) = plot(ax, fHz, m, 'Color', cMap(c,:), 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%s %s', stateName, cluName));
    end
end

% --- helper: draw one state's dendrogram with cluster colors chosen by delta/beta ---
function [order, clusters, cMap] = drawStateTree(freqData, axTree, axMat, f, maxClusters, stateColors, bands, stateName)
% The cluster with the highest mean delta/beta gets stateColors(2,:) (HighDB),
% every other cluster gets stateColors(1,:) (LowDB).
% cMap(i,:) is the color of cluster i, so it can be reused for the boxes on the matrix.

    C = corrcoef(freqData.normsPxx);
    dendroArgs = {'linkMetric','euclidean', 'linkMethod','ward', 'maxClusters',maxClusters, ...
                  'toPlotBinaryTree',1, 'plotOrderLabels',0};

    % pass 1: cluster membership only, drawn into a temporary figure that is closed right away
    fTmp  = figure('Visible','off');
    axTmp = [axes(fTmp, 'Position',[0.05 0.1 0.20 0.8]), axes(fTmp, 'Position',[0.30 0.1 0.60 0.8])];
    [~, ~, clusters] = DendrogramMatrix(C, dendroArgs{:}, 'figureHandle',fTmp, 'hDendro',axTmp);
    close(fTmp);

    % mean delta/beta of each cluster (from the raw spectra)
    fHz = freqData.freqHz(:);
    pxx = freqData.sPxx;                                   % expected [nFreq x nBins]
    if size(pxx,1) ~= numel(fHz) && size(pxx,2) == numel(fHz)
        pxx = pxx';
    end
    pDelta = fHz >= bands.delta(1) & fHz < bands.delta(2);
    pBeta  = fHz >= bands.beta(1)  & fHz < bands.beta(2);
    db     = mean(pxx(pDelta,:), 1) ./ mean(pxx(pBeta,:), 1);   % 1 x nBins
    nClu   = max(clusters);
    meanDB = arrayfun(@(c) mean(db(clusters == c)), 1:nClu);
    [~, rankIdx] = sort(meanDB);                           % ascending: last = highest delta/beta
    cMap = repmat(stateColors(1,:), nClu, 1);              % LowDB
    cMap(rankIdx(end),:) = stateColors(2,:);               % HighDB
    fprintf('%s: mean delta/beta per cluster = %s -> HighDB is cluster %d\n', ...
        stateName, mat2str(meanDB, 3), rankIdx(end));
    cMapFull = [cMap; lines(64)];                          % padded, in case DendrogramMatrix indexes past nClu

    % pass 2: the real tree, with the state colors
    [~, order, clusters] = DendrogramMatrix(C, dendroArgs{:}, 'figureHandle',f, ...
        'hDendro',[axTree, axMat], 'cMap',cMapFull);
end

% --- helper: replicates DendrogramMatrix's internal boundary logic ---
function [pStart, pEnd, blockClusters] = getClusterBoundaries(clusters, order)
    ordClu = clusters(order);
    ordClu = ordClu(:)';
    pStart = [1, find(diff(ordClu) ~= 0) + 1];
    pEnd   = [pStart(2:end) - 1, numel(ordClu)];
    blockClusters = ordClu(pStart);     % cluster id of each block -> row of the color map
end

 

%% panels  C, D, E, F
%% prepare pData parameters
% clearvars
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%
trialsIdx=WA.recTable.ForPaper==1;
recs = WA.recTable(trialsIdx,:);

% parameters states:: 
medianFiltWin= 1000*5; %default = 1000*5
longOrdFiltWin = 1000*200; %defualt = 1000*200
longOrdFiltOrd = 0.5; %default = 0.4
minStateDur= 1000*7; % defualt = 7*1000
tStartDef=1000*60*60*5; % tstart and win for night recording.
winDef = 1000*60*60*3; 


%parameters D/b
movWin = 5000;
movOLWin = 4500;

load([analysisFolderWake filesep 'wakeSleepData.mat'])
pDataOld = pData;

pData = recs(:,["Animal","recNames","Remarks","isAwake","Date","Headstage","spikes"]);
pData.type = categorical(repmat({''}, height(pData), 1));
pData.type(contains(pData.recNames, 'Wake'))  = 'Wake';
pData.type(contains(pData.recNames, 'Night')) = 'Sleep';
pData.P2V = zeros(height(pData),1);
pData.bouts = cell(height(pData),1);
pData.HighDBdur = zeros(height(pData),1);
pData.LowDBdur = zeros(height(pData),1);
pData.totalDur = zeros(height(pData),1);
pData.HighDBAvg = zeros(height(pData),1);
pData.LowDBavg = zeros(height(pData),1);
pData.transPoint = zeros(height(pData),1);
pData.DBTh = zeros(height(pData),1);


% check for recordings allready calculated, assuming the parameters havent
% changed:
newparameters=1;
if ~newparameters
    [isOld, idx] = ismember(pData(:, {'Animal','recNames'}), ...
                          pDataOld(:, {'Animal','recNames'}), 'rows');
    % Copy over already-computed rows
    pData(isOld, :) = pDataOld(idx(isOld), :);

    rowsToProcess = find(~isOld);
else
    rowsToProcess = [1:height(pData)]';
end

%% run on all recs and calculate pData:

for i = rowsToProcess'
    recName=['Animal=' pData.Animal{i} ',recNames=' pData.recNames{i}];
    WA.setCurrentRecording(recName);
    DB=WA.getDelta2BetaRatio;
    if DB.parDBRatio.movOLWin ~=movOLWin | DB.parDBRatio.movWin ~=movWin
        WA.getDelta2BetaRatio('overwrite',1, ...
            'movWin',movWin, ...
            'movOLWin',movOLWin);
    end
    
    if pData.type(i)=='Sleep'
        %add: if the isAwake is more than 1 number, take it as tstart and
        %win in hrs.
        if isscalar(str2num(pData.isAwake{i}))
            tStart = tStartDef;
            win = winDef;
        else 
            waketimes = str2num(pData.isAwake{i});
            tStart = waketimes(1)*60*60*1000;
            win = waketimes(2)*60*60*1000;
        end 
        SS = WA.getStatesAC( ... % get the states timings
            'overwrite',    1,...
            'medianFiltWin',medianFiltWin, ...
            'longOrdFiltWin', longOrdFiltWin, ...
            'longOrdFiltOrd', longOrdFiltOrd,...
            'minStateDur',minStateDur, ...
            'tStart', tStart, ... % three hrs into the rec
            'win',win ... % Three hrs of sleep
            );
        
        % WA.plotStatesAC
        AC = WA.getDelta2BetaAC('overwrite',1); % get the AC
        pData.P2V(i) = AC.peak2VallyDiff;
        
       

    elseif pData.type(i)=='Wake'
        if isscalar(pData.isAwake{i})
            tStart = 0;
            win = WA.currentDataObj.recordingDuration_ms;
        else 
            waketimes = str2num(pData.isAwake{i});
            tStart = waketimes(1)*60*60*1000;
            win = waketimes(2)*60*60*1000;
        end

        SS = WA.getStatesAC( ... % get the states timings
            'overwrite',    1,...
            'medianFiltWin',medianFiltWin, ...
            'longOrdFiltWin', longOrdFiltWin, ...
            'longOrdFiltOrd', longOrdFiltOrd,...
            'minStateDur',minStateDur, ...
            'tStart', tStart, ... % three hrs into the rec
            'win',win ... % Three hrs of sleep
            );
        
        % WA.plotStatesAC
        try
            AC = WA.getDelta2BetaAC('overwrite',1); % get the AC
            pData.P2V(i) = AC.peak2VallyDiff;
        catch
            disp('getDelta2BetaAC was not be completed for this recording. moving forward. ')
            pData.P2V(i) = 0;
        end

      
    end
    % calculate the duration of each states
    % get the th:
    pData.DBTh(i) = SS.Th;

    % -------------------------------------------------------------------
       % get the transition point - ONLY IF ALREADY RAN THE PCA DATAAAA!!!
       curFreqBand = WA.getFreqBandDetection;
       if ~isempty(curFreqBand.crossFreq)
           pData.transPoint(i) = curFreqBand.crossFreq;
       end
    % -------------------------------------------------------------------

    %get bouts
    BS = WA.getStateBouts('overwrite',1);

    pData.bouts{i}= BS.bouts(:,["state","label","duration_s","AvgDB"]);

    % get sum of state and rec duration
    pData.HighDBdur(i) = sum(BS.bouts.duration_s(BS.bouts.state==1));
    pData.LowDBdur(i) = sum(BS.bouts.duration_s(BS.bouts.state==2));
    pData.totalDur(i) = sum(BS.bouts.duration_s);
    
    % get the mean DB for each state:
    pData.HighDBAvg(i) = mean(BS.bouts.AvgDB(BS.bouts.state==1),'omitnan');
    pData.LowDBavg(i) = mean(BS.bouts.AvgDB(BS.bouts.state==2),'omitnan');
      

end

save([analysisFolderWake filesep 'wakeSleepData.mat'], "pData");

%% load pData
clearvars
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
load([analysisFolderWake filesep 'wakeSleepData.mat'])

%% plot AC
f = figure;
violinplot(pData,"type",'P2V');
hold on;
swarmchart(pData.type, pData.P2V, 15, 'k', 'filled', 'XJitterWidth', 0.3);
hold off;
ylabel('Periodicity (P2V)')

% save figure:
set(f,'PaperPosition',[1 5 2 2]);
fileName = [analysisFolderWake filesep 'sleepWakeAC'];
print(fileName,'-dpdf','-r300');

%% plot duration

% pool all inner tables into one long table, tagging each row with its Animal
allBouts = table();
for i = 1:height(pData)
    curBouts = pData.bouts{i};              % nx4 table for this row
    curBouts.Animal = repmat(pData.Animal(i), height(curBouts), 1);
    curBouts.recNames = repmat(pData.recNames(i), height(curBouts), 1);
    curBouts.type = repmat(pData.type(i), height(curBouts), 1);
    allBouts = [allBouts; curBouts];
end

allBouts.label = categorical(allBouts.label); % fix category order if desired

%% previous ver: plot bouts DOTS- sleep vs wake
% allBouts.group = categorical(strcat(string(allBouts.type), "-", string(allBouts.label)));
% 
% groupList = categories(allBouts.group);  % should be: Sleep_HighDB, Sleep_LowDB, Wake_HighDB, Wake_LowDB
% 
% % base colors: blue-family for HighDB (1,3), orange-family for LowDB (2,4)
% colorMat = [
%     0.85 0.33 0.10;   % 1: Sleep_HighDB  - orange (dark)
%     0.00 0.45 0.74;   % 2: Sleep_LowDB   - blue (dark)
%     0.93 0.60 0.40;   % 3: Wake_HighDB   - orange (light)
%     0.40 0.70 0.90    % 4: Wake_LowDB    - blue (light)
% ];
% 
% [~, idx] = ismember(allBouts.group, groupList);
% pointColors = colorMat(idx, :);
% 
% figure;
% swarmchart(allBouts.group, allBouts.duration_s, 5, pointColors, 'filled');
% ylabel('Bout duration (s)');
% xlabel('State');
% title('Bout duration distribution');
% yscale("log")
% 
% % save figure:
% fileName = [analysisFolderWake filesep 'BoutDurationAll'];
% print(fileName,'-dpdf','-r300');
% 
% % per animals
% allBouts.Animal = categorical(allBouts.Animal);
% animalList = categories(allBouts.Animal);
% groupList = categories(allBouts.group);
% cmap = lines(numel(animalList));
% 
% nAnimals = numel(animalList);
% offsets = linspace(-0.3, 0.3, nAnimals);  % spread animals within each group's x-slot
% 
% figure;
% hold on;
% h = gobjects(nAnimals,1);
% for i = 1:nAnimals
%     idx = allBouts.Animal == animalList{i};
% 
%     % numeric x position = group code + this animal's offset
%     xNum = double(allBouts.group(idx)) + offsets(i);
% 
%     h(i) = swarmchart(xNum, allBouts.duration_s(idx), 10, cmap(i,:), 'filled', ...
%         'XJitterWidth', 0.15);
% end
% 
% % set x-ticks to group names at the center of each cluster
% xticks(1:numel(groupList));
% xticklabels(groupList);
% xlim([0.5, numel(groupList)+0.5]);
% 
% legend(h, animalList, 'Location', 'best');
% ylabel('Bout duration (s)');
% xlabel('State');
% title('Bout duration distribution, per animal');allBouts.Animal = categorical(allBouts.Animal);
% animalList = categories(allBouts.Animal);
% groupList = categories(allBouts.group);
% cmap = lines(numel(animalList));
% 
% % save figure:
% fileName = [analysisFolderWake filesep 'BoutDurationPerAnimal'];
% print(fileName,'-dpdf','-r300');

%% Single-animal bout duration distributions: violin outline + swarm points
% Figure 1: vertical (groups on x-axis, duration on y-axis) -- like the original swarmchart
% Figure 2: horizontal, all 4 groups overlaid on the same row

thisAnimal = "PV157";   % <-- set to the animal you want to plot
allBouts.group = categorical(strcat(string(allBouts.type), "-", string(allBouts.label)));
% --- filter to this animal ---
isAnimal = allBouts.Animal == thisAnimal;
subBouts = allBouts(isAnimal, :);

groupList = categories(allBouts.group);   % Sleep_HighDB, Sleep_LowDB, Wake_HighDB, Wake_LowDB
nGroups   = numel(groupList);

colorMat = [
    0.85 0.33 0.10;   % 1: Sleep_HighDB - orange (dark)
    0.00 0.45 0.74;   % 2: Sleep_LowDB  - blue (dark)
    0.93 0.60 0.40;   % 3: Wake_HighDB  - orange (light)
    0.40 0.70 0.90    % 4: Wake_LowDB   - blue (light)
];

%% ---------- Figure 1: vertical, one column per group ----------
% violinHalfWidth = 0.4;    % max half-width of the density outline
% jitterHalfWidth = 0.35;   % max half-width for scattered points
% 
% fi = figure; hold on;
% 
% for i = 1:nGroups
%     grpData = subBouts.duration_s(subBouts.group == groupList{i});
%     if isempty(grpData)
%         continue
%     end
% 
%     % --- violin-style outline: KDE of duration, mirrored around baseline x = i ---
%     [f, xi] = ksdensity(grpData);
%     f = f / max(f) * violinHalfWidth;
%     xLeft  = i - f;
%     xRight = i + f;
%     fill([xRight, fliplr(xLeft)], [xi, fliplr(xi)], colorMat(i,:), ...
%         'FaceAlpha', 0.20, 'EdgeColor', colorMat(i,:), 'LineWidth', 1.3);
% 
%     % --- individual points, jittered in x, width scaled by local density ---
%     % fPoint = interp1(xi, f, grpData, 'linear', 0);
%     % jitterWidth = min(fPoint, jitterHalfWidth);
%     % xJitter = i + (rand(size(grpData)) - 0.5) .* 2 .* jitterWidth;
%     % scatter(xJitter, grpData, 10, colorMat(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
% end
% 
% xticks(1:nGroups);
% xticklabels(strrep(groupList, '_', ' '));
% xlim([0.5, nGroups + 0.5]);
% ylabel('Bout duration (s)');
% xlabel('State');
% title(sprintf('Bout duration distribution -- %s', thisAnimal));
% hold off;
% % yscale("log")
% % ylim([-100 200])
% % save figure:
%  set(fi,'PaperPosition',[1 5 3 3]);
% fileName = [analysisFolderWake filesep 'BoutDurationoneAnimalExamplecut'];
% print(fileName,'-dpdf','-r300');
% 
% 
% %% ---------- Figure 2: horizontal half-violin (raincloud-style), by TYPE, overlaid ----------
% allBouts.type = categorical(allBouts.type);   % ensure categorical (Sleep / Wake)
% subBouts.type = categorical(subBouts.type);
% typeList = categories(allBouts.type);
% nTypes   = numel(typeList);
% 
% typeColorMat = [
%     0.20 0.20 0.60;   % Sleep - dark indigo
%     0.90 0.55 0.10    % Wake  - orange
% ];   % edit/reorder if categories(allBouts.type) comes out in a different order
% 
% densityHeight = 0.4;    % how tall the density hump rises above baseline
% jitterHeight  = 0.4;   % how far points scatter below baseline
% baseline = 1;            % single shared row
% 
% fi = figure; hold on;
% h = gobjects(nTypes, 1);
% 
% for i = 1:nTypes
%     tData = subBouts.duration_s(subBouts.type == typeList{i});
%     if isempty(tData)
%         continue
%     end
% 
%     % --- upper half only: density hump rising from baseline (histogram-like) ---
%     [f, xi] = ksdensity(tData);
%     f = f / max(f) * densityHeight;
%     yTop = baseline + f;
%     h(i) = fill([xi, fliplr(xi)], [yTop, baseline*ones(size(xi))], typeColorMat(i,:), ...
%         'FaceAlpha', 0.25, 'EdgeColor', typeColorMat(i,:), 'LineWidth', 1.3);
% 
%     % --- points scattered below baseline only (raincloud strip) ---
%     fPoint = interp1(xi, f, tData, 'linear', 0);
%     jitterWidth = min(fPoint, jitterHeight);
%     yJitter =  baseline + abs((rand(size(tData)) - 0.5)) .* 2 .* jitterWidth;
%     scatter(tData, yJitter, 5, typeColorMat(i,:), 'filled', 'MarkerFaceAlpha', 0.5);
% end
% 
% yticks([]);
% ylim([baseline - 0.01, baseline + densityHeight + 0.01]);
% xlabel('Bout duration (s)');
% title(sprintf('Bout duration distribution (overlaid, by type) -- %s', thisAnimal));
% legend(h, strrep(typeList, '_', ' '), 'Location', 'best');
% hold off;
% 
% % save figure:
% set(fi,'PaperPosition',[1 5 4 4]);
% fileName = [analysisFolderWake filesep 'BoutDurationoneAnimalExampleHorizontalSW'];
% print(fileName,'-dpdf','-r300');

%% ---------- Figure 3: vertical split violins, log scale ----------
% Two violins (HighDB, LowDB). Each one is split down the middle:
%   left half = Wake, right half = Sleep.
% No dots, same fill/outline style as the previous violin figures.
% The density is estimated on log10(duration), so the shape is correct on a log axis.
%
% Runs after the sections above (needs subBouts / allBouts, thisAnimal, analysisFolderWake).

plotData = subBouts;  titleTag = thisAnimal;        % single animal
% plotData = allBouts; titleTag = "all animals";    % pooled across animals

violinHalfWidth = 0.4;     % max half-width of each half-violin
commonScale     = false;   % false: every half is scaled to the same max width (as in the previous figures)
                           % true : all halves share one scale, so widths are comparable between halves

labelNames = ["HighDB", "LowDB"];   % one violin per label
typeNames  = ["Wake", "Sleep"];     % left half, right half
sideSign   = [-1, 1];               % -1 = left of the centre line, +1 = right

% same palette as colorMat above; rows = label, columns = type (Wake, Sleep)
halfColors = { [0.93 0.60 0.40], [0.85 0.33 0.10];     % HighDB: Wake light orange, Sleep dark orange
               [0.40 0.70 0.90], [0.00 0.45 0.74] };   % LowDB : Wake light blue,   Sleep dark blue

labelStr = string(plotData.label);
typeStr  = string(plotData.type);
nLabels  = numel(labelNames);
nTypes   = numel(typeNames);

% --- density of every half, estimated in log10 space ---
dens = cell(nLabels, nTypes);   % density values
yy   = cell(nLabels, nTypes);   % matching durations (s)
for i = 1:nLabels
    for j = 1:nTypes
        sel = contains(labelStr, labelNames(i), 'IgnoreCase', true) & ...
              contains(typeStr,  typeNames(j),  'IgnoreCase', true);
        d = plotData.duration_s(sel);
        d = d(d > 0 & isfinite(d));              % log scale: positive values only
        if numel(d) < 2
            warning('Not enough bouts for %s - %s, half-violin skipped.', labelNames(i), typeNames(j));
            continue
        end
        [f, xi]   = ksdensity(log10(d));
        dens{i,j} = f;
        yy{i,j}   = 10.^xi;
    end
end

if commonScale
    globalMax = max(cellfun(@(f) max([f, 0]), dens(:)));
end

% --- plot ---
fi = figure; hold on;
h = gobjects(nLabels, nTypes);

for i = 1:nLabels
    for j = 1:nTypes
        f = dens{i,j};
        if isempty(f)
            continue
        end
        if commonScale
            w = f / globalMax * violinHalfWidth;
        else
            w = f / max(f) * violinHalfWidth;
        end

        xOuter = i + sideSign(j) * w;            % curved outer edge
        xInner = i * ones(size(w));              % straight edge on the centre line
        h(i,j) = fill([xOuter, xInner], [yy{i,j}, fliplr(yy{i,j})], halfColors{i,j}, ...
            'FaceAlpha', 0.20, 'EdgeColor', halfColors{i,j}, 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s - %s', labelNames(i), typeNames(j)));
    end
end

set(gca, 'YScale', 'log');
xticks(1:nLabels);
xticklabels(labelNames);
xlim([0.5, nLabels + 0.5]);
ylabel('Bout duration (s)');
xlabel('State');
title(sprintf('Bout duration distribution -- %s', titleTag));
valid = isgraphics(h);
legend(h(valid), 'Location', 'best');
hold off;

% save figure:
set(fi, 'PaperPosition', [1 5 3 3]);
fileName = [analysisFolderWake filesep 'BoutDurationSplitViolinLog'];
print(fileName, '-dpdf', '-r300');
%% ---------- Setup ----------
% Assumes 'allBouts' (with .Animal, .type, .label, .group, .duration_s) and
% 'thisAnimal' are already defined, as in the earlier scripts.
% thisAnimal = "PV157";
isAnimal = allBouts.Animal == thisAnimal;
subBouts = allBouts(isAnimal, :);

groupList = categories(allBouts.group);   % Sleep_HighDB, Sleep_LowDB, Wake_HighDB, Wake_LowDB

colorMat = [
    0.85 0.33 0.10;   % Sleep_HighDB - orange (dark)
    0.00 0.45 0.74;   % Sleep_LowDB  - blue (dark)
    0.93 0.60 0.40;   % Wake_HighDB  - orange (light)
    0.40 0.70 0.90    % Wake_LowDB   - blue (light)
];

densityHeight = 0.4;
jitterHeight  = 0.4;
baseline      = 1;

highGroups = groupList(contains(groupList, 'High', 'IgnoreCase', true));
lowGroups  = groupList(contains(groupList, 'Low',  'IgnoreCase', true));

%% ---------- Figure 1: raincloud plot, split into High / Low subplots ----------
fi = figure;

ax1 =subplot(2,1,1); hold on;
hHigh = gobjects(numel(highGroups),1);
for i = 1:numel(highGroups)
    g = highGroups{i};
    color = colorMat(strcmp(groupList, g), :);
    data = subBouts.duration_s(subBouts.group == g);
    if isempty(data), continue; end
    hHigh(i) = drawRaincloudHalf(data, color, baseline, densityHeight, jitterHeight);
end
yticks([]); ylim([baseline - 0.01, baseline + densityHeight + 0.01]);
xlabel('Bout duration (s)');
title('High DB');
legend(hHigh, strrep(highGroups, '_', ' '), 'Location', 'best');

ax2 =subplot(2,1,2); hold on;
hLow = gobjects(numel(lowGroups),1);
for i = 1:numel(lowGroups)
    g = lowGroups{i};
    color = colorMat(strcmp(groupList, g), :);
    data = subBouts.duration_s(subBouts.group == g);
    if isempty(data), continue; end
    hLow(i) = drawRaincloudHalf(data, color, baseline, densityHeight, jitterHeight);
end
yticks([]); ylim([baseline - 0.01, baseline + densityHeight + 0.01]);
xlabel('Bout duration (s)');
title('Low DB');
linkaxes([ax1,ax2],'x');
xlim([-20 500])
legend(hLow, strrep(lowGroups, '_', ' '), 'Location', 'best');

sgtitle(sprintf('Bout duration distribution -- %s', thisAnimal));

% save figure
 set(fi,'PaperPosition',[1 5 3 3]);
fileName = [analysisFolderWake filesep 'BoutDurationOneAnimalSW_HighLow'];
print(fileName, '-dpdf', '-r300');

%% ---------- Figure 2: per-animal Sleep-minus-Wake density difference, High / Low subplots ----------
animalList = categories(categorical(allBouts.Animal));
nAnimals   = numel(animalList);
nGrid      = 200;

sleepHighName = groupList{find(contains(groupList,'Sleep','IgnoreCase',true) & contains(groupList,'High','IgnoreCase',true), 1)};
wakeHighName  = groupList{find(contains(groupList,'Wake', 'IgnoreCase',true) & contains(groupList,'High','IgnoreCase',true), 1)};
sleepLowName  = groupList{find(contains(groupList,'Sleep','IgnoreCase',true) & contains(groupList,'Low', 'IgnoreCase',true), 1)};
wakeLowName   = groupList{find(contains(groupList,'Wake', 'IgnoreCase',true) & contains(groupList,'Low', 'IgnoreCase',true), 1)};

fi = figure;

% --- High DB: Sleep_HighDB minus Wake_HighDB, per animal ---
ax1 = subplot(2,1,1); hold on;
highMask = allBouts.group == sleepHighName | allBouts.group == wakeHighName;
xiHigh = linspace(min(allBouts.duration_s(highMask)), max(allBouts.duration_s(highMask)), nGrid);
diffMatHigh = nan(nAnimals, nGrid);
for a = 1:nAnimals
    sleepData = allBouts.duration_s(contains(allBouts.Animal, animalList{a}) & allBouts.group == sleepHighName);
    wakeData  = allBouts.duration_s(contains(allBouts.Animal, animalList{a}) & allBouts.group == wakeHighName);
    if numel(sleepData) < 2 || numel(wakeData) < 2
        continue   % not enough bouts for this animal to estimate a density
    end
    fSleep = ksdensity(sleepData, xiHigh);
    fWake  = ksdensity(wakeData, xiHigh);
    diffMatHigh(a, :) = fSleep - fWake;
    plot(xiHigh, diffMatHigh(a,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end
plot(xiHigh, mean(diffMatHigh, 1, 'omitnan'), 'k', 'LineWidth', 2.5);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Bout duration (s)');
ylabel('Density diff (Sleep - Wake)');
title('High DB');

% --- Low DB: Sleep_LowDB minus Wake_LowDB, per animal ---
ax2=subplot(2,1,2); hold on;
lowMask = allBouts.group == sleepLowName | allBouts.group == wakeLowName;
xiLow = linspace(min(allBouts.duration_s(lowMask)), max(allBouts.duration_s(lowMask)), nGrid);
diffMatLow = nan(nAnimals, nGrid);
for a = 1:nAnimals
    sleepData = allBouts.duration_s(contains(allBouts.Animal, animalList{a}) & allBouts.group == sleepLowName);
    wakeData  = allBouts.duration_s(contains(allBouts.Animal, animalList{a}) & allBouts.group == wakeLowName);
    if numel(sleepData) < 2 || numel(wakeData) < 2
        continue
    end
    fSleep = ksdensity(sleepData, xiLow);
    fWake  = ksdensity(wakeData, xiLow);
    diffMatLow(a, :) = fSleep - fWake;
    plot(xiLow, diffMatLow(a,:), 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
end
plot(xiLow, mean(diffMatLow, 1, 'omitnan'), 'k', 'LineWidth', 2.5);
yline(0, '--', 'Color', [0.5 0.5 0.5]);
xlabel('Bout duration (s)');
ylabel('Density diff (Sleep - Wake)');
title('Low DB');
linkaxes([ax1,ax2],'x')
xlim([0 150])
sgtitle('Per-animal Sleep - Wake bout-duration density difference');

% save figure
% save figure
 set(fi,'PaperPosition',[1 5 3 3]);
fileName2 = [analysisFolderWake filesep 'BoutDurationDiffSleepWake_HighLow_allAnimals'];
print(fileName2, '-dpdf', '-r300');

%% ---------- Local functions ----------
function h = drawRaincloudHalf(data, color, baseline, densityHeight, jitterHeight)
% Draws one upper-half density hump (rising from baseline) plus a
% density-width-scaled jittered point cloud filling the area under it.
    [f, xi] = ksdensity(data);
    f = f / max(f) * densityHeight;
    yTop = baseline + f;
    h = fill([xi, fliplr(xi)], [yTop, baseline*ones(size(xi))], color, ...
        'FaceAlpha', 0.25, 'EdgeColor', color, 'LineWidth', 1.3);

    % fPoint = interp1(xi, f, data, 'linear', 0);
    % jitterWidth = min(fPoint, jitterHeight);
    % yJitter = baseline + abs((rand(size(data)) - 0.5)) .* 2 .* jitterWidth;
    % scatter(data, yJitter, 5, color, 'filled', 'MarkerFaceAlpha', 0.5);
end


%% Fano factor & coefficient of variation of bout duration, by group
% Fano factor = variance / mean
% CV          = std / mean
%
% Produces two tables:
%   statsPerAnimal - one row per Animal x group
%   statsAllAnimals - one row per group, pooling bouts across all animals

allBouts.group  = categorical(allBouts.group);
allBouts.Animal = categorical(allBouts.Animal);

% ---------- per animal, per group ----------
statsPerAnimal = groupsummary(allBouts, {'Animal', 'group'}, ...
    {"mean","var","std"}, 'duration_s');

statsPerAnimal.Properties.VariableNames{'mean_duration_s'} = 'meanDur';
statsPerAnimal.Properties.VariableNames{'var_duration_s'}  = 'varDur';
statsPerAnimal.Properties.VariableNames{'std_duration_s'} = 'stdDur';
statsPerAnimal.fanoFactor = statsPerAnimal.varDur ./ statsPerAnimal.meanDur;
statsPerAnimal.CV         = statsPerAnimal.stdDur ./ statsPerAnimal.meanDur;

disp('--- Fano factor & CV, per animal per group ---');
disp(statsPerAnimal);

% ---------- pooled across all animals, per group ----------
statsAllAnimals = groupsummary(allBouts, {'group'}, ...
    {"mean", "var", "std"}, 'duration_s');

statsAllAnimals.Properties.VariableNames{'mean_duration_s'} = 'meanDur';
statsAllAnimals.Properties.VariableNames{'var_duration_s'}  = 'varDur';
statsAllAnimals.Properties.VariableNames{'std_duration_s'} = 'stdDur';
statsAllAnimals.fanoFactor = statsAllAnimals.varDur ./ statsAllAnimals.meanDur;
statsAllAnimals.CV         = statsAllAnimals.stdDur ./ statsAllAnimals.meanDur;

disp('--- Fano factor & CV, pooled across all animals, per group ---');
disp(statsAllAnimals);

% ---------- optional: save both tables ----------
writetable(statsPerAnimal,  fullfile(analysisFolderWake, 'fanoCV_perAnimal.csv'));
writetable(statsAllAnimals, fullfile(analysisFolderWake, 'fanoCV_allAnimals.csv'));


%% plot percantage
pData.HighDBPer = pData.HighDBdur./pData.totalDur;
pData.LowDBPer = pData.LowDBdur./pData.totalDur;
% sum = pData.HighDBPer+pData.LowDBPer;

fi = figure;
hold on;

sleepType = pData.type == 'Sleep';

swarmchart(ones(nnz(sleepType),1), pData.HighDBPer(sleepType)*100, 20, 'filled','XJitter','density');
swarmchart(2*ones(nnz(~sleepType),1), pData.HighDBPer(~sleepType)*100, 20, 'filled');

xticks([1 2]);
xticklabels(["Sleep", "Wake"]);
xlim([0.5 2.5]);
yline(50, '--','Color',[0.4 0.4 0.4]);

ylabel('High DB percentage of the rec time');
title('Bout duration percentage');

% save figure
set(fi,'PaperPosition',[1 5 2 2]);
fileName2 = [analysisFolderWake filesep 'BoutsPercentagesAllAnimals'];
print(fileName2, '-dpdf', '-r300');
%%  old - plot mean DB for each state
% groupOrder = {'Sleep','Wake'};
% pData.type = categorical(pData.type, groupOrder);
% 
% % Reshape to long format
% longType  = [pData.type; pData.type];
% longValue = [pData.HighDBAvg; pData.LowDBavg];
% longLabel = [repmat({'highDB'}, height(pData), 1); repmat({'lowDB'}, height(pData), 1)];
% longLabel = categorical(longLabel, {'highDB','lowDB'});
% 
% % Base x position from type, with small offset per sub-group
% baseX = double(longType);
% offset = 0.30;
% xPos = baseX - offset .* (longLabel == 'highDB') + offset .* (longLabel == 'lowDB');
% 
% colorMat = [
%     0.85 0.33 0.10;   % Sleep_HighDB  - orange (dark)
%     0.00 0.45 0.74;   % Sleep_LowDB   - blue (dark)
%     0.93 0.60 0.40;   % Wake_HighDB   - orange (light)
%     0.40 0.70 0.90    % Wake_LowDB    - blue (light)
% ];
% 
% % Per-point color matrix
% pointColors = zeros(numel(xPos), 3);
% pointColors(longType=='Sleep' & longLabel=='highDB', :) = repmat(colorMat(1,:), sum(longType=='Sleep' & longLabel=='highDB'), 1);
% pointColors(longType=='Sleep' & longLabel=='lowDB',  :) = repmat(colorMat(2,:), sum(longType=='Sleep' & longLabel=='lowDB'),  1);
% pointColors(longType=='Wake'  & longLabel=='highDB', :) = repmat(colorMat(3,:), sum(longType=='Wake'  & longLabel=='highDB'), 1);
% pointColors(longType=='Wake'  & longLabel=='lowDB',  :) = repmat(colorMat(4,:), sum(longType=='Wake'  & longLabel=='lowDB'),  1);
% 
% figure; hold on;
% swarmchart(xPos, longValue, 25, pointColors, 'filled', 'o', 'XJitterWidth', 0.25);
% % Xticks at each sub-group position, labeled directly
% xticks([1-offset, 1+offset, 2-offset, 2+offset]);
% xticklabels({'Sleep-HighDB','Sleep-LowDB','Wake-HighDB','Wake-LowDB'});
% xtickangle(0);
% 
% xlim([0.5, numel(groupOrder)+0.5]);
% % ylabel('Avg DB Ratio');
% % save figure
% fileName2 = [analysisFolderWake filesep 'bbAvgperState'];
% print(fileName2, '-dpdf', '-r300');
%%  plot mean DB for each state
pData.avgDB = mean([pData.HighDBAvg,pData.LowDBavg],2);

fi = figure;
swarmchart(categorical(pData.type), pData.avgDB, 25, 'filled')
xlabel('Type')
ylabel('avg dB')
title('avgDB by State')
box off

% % save figure
set(fi,'PaperPosition',[1 5 2 2]);
fileName2 = [analysisFolderWake filesep 'bbAvgbytype'];
print(fileName2, '-dpdf', '-r300');

%%  plot Threshold for each type - sleep or wake
% pData.avgDB = mean([pData.HighDBAvg,pData.LowDBavg],2);

fi = figure;
swarmchart(categorical(pData.type), pData.DBTh, 25, 'filled')
xlabel('Type')
ylabel('D/B')
title('DB threshold betweh high and low DB by State')
box off

% % save figure
set(fi,'PaperPosition',[1 5 2 2]);
fileName2 = [analysisFolderWake filesep 'dbThperType'];
print(fileName2, '-dpdf', '-r300');
%%  plot Cross Frequency for each type - sleep or wake
% pData.avgDB = mean([pData.HighDBAvg,pData.LowDBavg],2);

fi = figure;
swarmchart(categorical(pData.type), pData.transPoint, 25, 'filled')
xlabel('Type')
ylabel('Frequency')
title('Cross Frequency between normsPXX')
box off

% % save figure
set(fi,'PaperPosition',[1 5 2 2]);
fileName2 = [analysisFolderWake filesep 'CrossFreq'];
print(fileName2, '-dpdf', '-r300');

%% Panels ??? - PCA 
clearvars
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

% parameters for freqbands
overwrite = 1;  % set to 1 to force recomputation of freq-band data / bouts

fMax = 80;
segmentLength = 1000;
binDuration = 5000;
WelchOL = 0.75;
dftPoints = 2^10;
winDuration = 1*60*60*1000;  % target window length: 1 hour, in ms
tStartDef=1000*60*60*5; % tstart for night recording.

%% get wake subset - all trials
wakeTrials = WA.recTable.ForPaper==1 & cellfun(@(x) contains(x, 'Wake'), WA.recTable.recNames);
wakeTable = WA.recTable(wakeTrials,:);

% running on a single trial, power spec - WAKE
wakeRows = cell(height(wakeTable), 1);
freqHzAllWake = [];

for i = 1:height(wakeTable)
    recName = ['Animal=' wakeTable.Animal{i} ',recNames=' wakeTable.recNames{i}];
    WA.setCurrentRecording(recName);
    if isscalar(wakeTable.isAwake{i})
        tStart = 0;
        win = WA.currentDataObj.recordingDuration_ms;
    else
        waketimes = str2num(wakeTable.isAwake{i});
        tStart = waketimes(1)*60*60*1000;
        win = waketimes(2)*60*60*1000;
        if win < winDuration
            win = winDuration;
        end
    end

    animal = wakeTable.Animal(i);

    curfreqData = WA.getFreqBandDetection(...
        'fMax',          fMax, ...
        'segmentLength', segmentLength, ...
        'binDuration',   binDuration, ...
        'WelchOL',       WelchOL, ...
        'dftPoints',     dftPoints, ...
        'tStart',        tStart, ...
        'win',           win, ...
        'overwrite',     overwrite ...
    );

    pxxData  = curfreqData.sPxx';      % [nBins x nFreqs]
    npxxData = curfreqData.normsPxx';
    nBins    = size(pxxData, 1);

    if isempty(freqHzAllWake)
        freqHzAllWake = curfreqData.freqHz;
    else
        assert(numel(curfreqData.freqHz) == numel(freqHzAllWake), ...
            'Frequency axis mismatch in rec %s!', recName);
    end

    % --- state labels from bouts (external source), matched by segment START time ---
    statesStruct = WA.getStateBouts;
    segLabels = assignStateLabels(curfreqData.times, statesStruct.bouts);

    wakeRows{i} = table(repmat(animal, nBins,1), repmat(wakeTable.recNames(i), nBins,1), ...
        pxxData, npxxData, curfreqData.clusters, curfreqData.times', segLabels, ...
        'VariableNames', {'animal','recs','Pxx','nPxx','clusters','times','label'});
end

WpcaT = vertcat(wakeRows{:});
save([analysisFolderWake filesep 'WakePCAdata.mat'], 'WpcaT', 'freqHzAllWake');

%% get sleep subset
sleepTrials = WA.recTable.ForPaper==1 & contains(WA.recTable.recNames,'Night');
sleepTable = WA.recTable(sleepTrials,:);

% get the data for the PCA - SLEEP
sleepRows = cell(height(sleepTable), 1);
freqHzAllSleep = [];

for i = 1:height(sleepTable)
    recName = ['Animal=' sleepTable.Animal{i} ',recNames=' sleepTable.recNames{i}];
    WA.setCurrentRecording(recName);
    animal = sleepTable.Animal(i);

    % curAC = WA.getDelta2BetaAC;
    % tStartSleep = curAC.tStartSleep + 0.5*60*60*1000;
    if isscalar(str2num(sleepTable.isAwake{i}))
        tStartSleep = tStartDef;
    else
        waketimes = str2num(sleepTable.isAwake{i});
        tStartSleep = waketimes(1)*60*60*1000;
    end
    win = winDuration;

    curfreqData = WA.getFreqBandDetection(...
        'fMax',          fMax, ...
        'segmentLength', segmentLength, ...
        'binDuration',   binDuration, ...
        'WelchOL',       WelchOL, ...
        'dftPoints',     dftPoints, ...
        'tStart',        tStartSleep, ...
        'win',           win, ...
        'overwrite',     overwrite ...
    );

    pxxData  = curfreqData.sPxx';
    npxxData = curfreqData.normsPxx';
    nBins    = size(pxxData, 1);

    if isempty(freqHzAllSleep)
        freqHzAllSleep = curfreqData.freqHz;
    else
        assert(numel(curfreqData.freqHz) == numel(freqHzAllSleep), ...
            'Frequency axis mismatch in rec %s!', recName);
    end

    % --- state labels from bouts (external source), matched by segment START time ---
    statesStruct = WA.getStateBouts;
    segLabels = assignStateLabels(curfreqData.times, statesStruct.bouts);

    sleepRows{i} = table(repmat(animal, nBins,1), repmat(sleepTable.recNames(i), nBins,1), ...
        pxxData, npxxData, curfreqData.clusters, curfreqData.times', segLabels, ...
        'VariableNames', {'animal','recs','Pxx','nPxx','clusters','times','label'});
end

SpcaT = vertcat(sleepRows{:});
save([analysisFolderWake filesep 'SleepPCAdata.mat'], 'SpcaT', 'freqHzAllSleep');

%% load data, remove freq 0, add normalized-segment Pxx

WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

load([analysisFolderWake filesep 'SleepPCAdata.mat'])
load([analysisFolderWake filesep 'WakePCAdata.mat'])

freqHzAll = freqHzAllSleep;
f_ind = freqHzAll > 0.5;
freqHzAll     = freqHzAll(f_ind);
freqHzAllWake = freqHzAllWake(f_ind);

SpcaT.Pxx  = SpcaT.Pxx(:, f_ind);
SpcaT.nPxx = SpcaT.nPxx(:, f_ind);
WpcaT.Pxx  = WpcaT.Pxx(:, f_ind);
WpcaT.nPxx = WpcaT.nPxx(:, f_ind);

SpcaT.pxxNormSegm = normalize(SpcaT.Pxx, 2);
WpcaT.pxxNormSegm = normalize(WpcaT.Pxx, 2);

SpcaT = SpcaT(:, {'animal','recs','Pxx','pxxNormSegm','nPxx','clusters','times','label'});
WpcaT = WpcaT(:, {'animal','recs','Pxx','pxxNormSegm','nPxx','clusters','times','label'});

%% 
function labels = assignStateLabels(segStartTimes, bouts)
% Assigns a behavioral-state label to each analysis segment based on
% which bout (from WA.getStatesBouts) contains the segment's START time.
%
%   segStartTimes : 1xN (or Nx1) segment start times [ms]
%   bouts         : table with startTime_ms, endTime_ms, label
%
% Segments whose start time falls outside every bout get 'unlabeled'.

    N = numel(segStartTimes);
    labels = repmat({'unlabeled'}, N, 1);
    segStart = segStartTimes(:);

    boutLabels = string(bouts.label); % normalizes char/cellstr/categorical/string

    for b = 1:height(bouts)
        inBout = segStart >= bouts.startTime_ms(b) & segStart < bouts.endTime_ms(b);
        labels(inBout) = {char(boutLabels(b))};
    end
end


%% plot spectrograms by clusters:
Wsum = groupsummary(WpcaT, {'animal', 'recs', 'clusters'}, "mean", {'Pxx','pxxNormSegm'});
Ssum = groupsummary(SpcaT, {'animal', 'recs', 'clusters'}, "mean", {'Pxx','pxxNormSegm'});

colorMat = [
    0.85 0.33 0.10;   % Sleep_HighDB - orange (dark)
    0.00 0.45 0.74;   % Sleep_LowDB  - blue (dark)
    0.93 0.60 0.40;   % Wake_HighDB  - orange (light)
    0.40 0.70 0.90    % Wake_LowDB   - blue (light)
];

labels = {'Sleep\_HighDB','Sleep\_LowDB','Wake\_HighDB','Wake\_LowDB'};

% --- Figure 1: pxxNormSegm ---
plotGroupedSpectra(freqHzAll, Ssum, Wsum, 'mean_pxxNormSegm', colorMat, labels, ...
    'Normalized Power', 'Wake/Sleep spectra (normalized)');

fileName2 = [analysisFolderWake filesep 'spectrogramByClustersNorm'];
print(fileName2, '-dpdf', '-r300');

% --- Figure 2: raw Pxx ---
plotGroupedSpectra(freqHzAll, Ssum, Wsum, 'mean_Pxx', colorMat, labels, ...
    'Power', 'Wake/Sleep spectra (raw Pxx)');
fileName2 = [analysisFolderWake filesep 'spectrogramByClustersRaw'];
print(fileName2, '-dpdf', '-r300');

function plotGroupedSpectra(freqHzAll, Ssum, Wsum, varName, colorMat, labels, yLbl, ttl)
    figure; hold on;
    h = gobjects(4,1);

    % Sleep clusters 1-2, Wake clusters 1-2 -> colorMat rows 1-4
    tables  = {Ssum, Ssum, Wsum, Wsum};
    clusters = [1, 2, 1, 2];

    for i = 1:4
        T = tables{i};
        c = clusters(i);
        data = T.(varName)(T.clusters==c, :)'; % freq x nLines

        % individual lines at 50% opacity
        for k = 1:size(data,2)
            p = plot(freqHzAll, data(:,k), 'LineWidth', 0.75);
            p.Color = [colorMat(i,:), 0.2]; % 4th value = alpha
        end

        % mean line, full opacity, on top
        h(i) = plot(freqHzAll, mean(data,2), ...
            'Color', colorMat(i,:), 'LineWidth', 2.5);
    end

    legend(h, labels, 'Location', 'best');
    xlabel('Frequency (Hz)');
    ylabel(yLbl);
    title(ttl);
    hold off;
end

%% plot spectra grouped by state LABEL (HighDB / LowDB) instead of clusters

useUnlabeled = false;   % set true to also plot the 'unlabeled' segments

if useUnlabeled
    Wlab = WpcaT;
    Slab = SpcaT;
else
    Wlab = WpcaT(string(WpcaT.label) ~= "unlabeled", :);
    Slab = SpcaT(string(SpcaT.label) ~= "unlabeled", :);
end

Wsum = groupsummary(Wlab, {'animal','recs','label'}, "mean", {'Pxx','pxxNormSegm'});
Ssum = groupsummary(Slab, {'animal','recs','label'}, "mean", {'Pxx','pxxNormSegm'});

% optional: drop animal/rec/label groups built from very few segments
minSegs = 1;                       % raise this (e.g. 10) if you want
Wsum = Wsum(Wsum.GroupCount >= minSegs, :);
Ssum = Ssum(Ssum.GroupCount >= minSegs, :);

% each row: {summary table, label to select, color, legend entry}
groupSpec = { ...
    Ssum, "HighDB", [0.85 0.33 0.10], 'Sleep\_HighDB'; ...
    Ssum, "LowDB",  [0.00 0.45 0.74], 'Sleep\_LowDB' ; ...
    Wsum, "HighDB", [0.93 0.60 0.40], 'Wake\_HighDB' ; ...
    Wsum, "LowDB",  [0.40 0.70 0.90], 'Wake\_LowDB'  };

if useUnlabeled
    groupSpec = [groupSpec; ...
        {Ssum, "unlabeled", [0.35 0.35 0.35], 'Sleep\_unlabeled'; ...
         Wsum, "unlabeled", [0.70 0.70 0.70], 'Wake\_unlabeled' }];
end

plotLabeledSpectra(freqHzAll, groupSpec, 'mean_pxxNormSegm', ...
    'Normalized Power', 'Wake/Sleep spectra by state label (normalized)');
fileName2 = [analysisFolderWake filesep 'spectrogramByLabelNorm'];
print(fileName2, '-dpdf', '-r300');

plotLabeledSpectra(freqHzAll, groupSpec, 'mean_Pxx', ...
    'Power', 'Wake/Sleep spectra by state label (raw Pxx)');
fileName2 = [analysisFolderWake filesep 'spectrogramByLabelRaw'];
print(fileName2, '-dpdf', '-r300');

function plotLabeledSpectra(freqHzAll, groupSpec, varName, yLbl, ttl)
    nG = size(groupSpec, 1);
    figure; hold on;
    h = gobjects(nG, 1);

    for i = 1:nG
        T   = groupSpec{i,1};
        lbl = groupSpec{i,2};
        col = groupSpec{i,3};

        sel = string(T.label) == lbl;
        if ~any(sel), continue; end

        data = T.(varName)(sel, :)';        % freq x nLines

        % individual animal/rec lines, faded
        for k = 1:size(data,2)
            p = plot(freqHzAll, data(:,k), 'LineWidth', 0.65);
            p.Color = [col, 0.2];
        end

        % group mean on top
        h(i) = plot(freqHzAll, mean(data, 2, 'omitnan'), ...
            'Color', col, 'LineWidth', 2.5);
    end

    valid = isgraphics(h);
    legend(h(valid), groupSpec(valid,4), 'Location', 'best');
    xlabel('Frequency (Hz)');
    ylabel(yLbl);
    yscale log
    title(ttl);
    hold off;
end


%% PCA Part - 

out = plotPCAwakeBased(WpcaT, SpcaT, freqHzAll, 'panelSize',[300 220], ...
    dataField='pxxNormSegm');   % bigger maps
