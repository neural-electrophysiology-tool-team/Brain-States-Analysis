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

% plot them together with the dendro tree:
f = figure('Position',[100 100 798 601]);
margin = 0.05;
leftWidth  = 0.18;
cbarWidth  = 0.025;
cbarGap    = 0.015;
usableWidth  = 1 - 2*margin;
usableHeight = 1 - 2*margin;
ax3Width = usableWidth - leftWidth - cbarWidth - cbarGap;
ax1 = axes('Parent', f, 'Position', [margin, margin + usableHeight/2, leftWidth, usableHeight/2]);
ax2 = axes('Parent', f, 'Position', [margin, margin,                 leftWidth, usableHeight/2]);
ax3 = axes('Parent', f, 'Position', [margin + leftWidth, margin, ax3Width, usableHeight]);
set(ax1, 'PositionConstraint', 'innerposition');
set(ax2, 'PositionConstraint', 'innerposition');
set(ax3, 'PositionConstraint', 'innerposition');

maxDendroClusters = 2;
cMap = lines(64);   % shared colormap so tree colors == box colors, wake and sleep alike

% wake tree into ax1
wakeC = corrcoef(wakeFreqData.normsPxx);
[~, orderW, clustersW] = DendrogramMatrix(wakeC,'linkMetric','euclidean','linkMethod','ward', ...
                        'maxClusters',maxDendroClusters, 'figureHandle',f, ...
                        'toPlotBinaryTree',1,'hDendro',[ax1,ax3],'plotOrderLabels',0,'cMap',cMap);
ylabel('Wake')

% sleep tree into ax2
sleepC = corrcoef(sleepFreqData.normsPxx);
[~, orderS, clusterS] = DendrogramMatrix(sleepC,'linkMetric','euclidean','linkMethod','ward', ...
                        'maxClusters',maxDendroClusters, 'figureHandle',f, ...
                        'toPlotBinaryTree',1,'hDendro',[ax2,ax3],'plotOrderLabels',0,'cMap',cMap);
ylabel('Sleep')

% Combined ordered matrix into ax3
normsPXXboth = [wakeFreqData.normsPxx, sleepFreqData.normsPxx];
corrMatboth = corrcoef(normsPXXboth);
nWake = length(orderW);
nSleep = length(orderS);
orderBoth = [orderW, orderS + nWake];
reordered_corr = corrMatboth(orderBoth, orderBoth);

delete(ax3)
ax3 = axes('Parent', f, 'Position', [margin + leftWidth, margin, ax3Width, usableHeight]);
set(ax3, 'PositionConstraint', 'innerposition');
axes(ax3)
imagesc(reordered_corr);
axis square;
colormap(parula);
clim([-1 1]);

cb = colorbar(ax3);
cb.Position = [margin + leftWidth + ax3Width + cbarGap*0.3, margin + 0.15*usableHeight, cbarWidth, 0.7*usableHeight];

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

% --- 3) redraw cluster-state boxes from wake/sleep trees onto ax3 ---
[pStartW, pEndW] = getClusterBoundaries(clustersW, orderW);
[pStartS, pEndS] = getClusterBoundaries(clusterS, orderS);
pStartS = pStartS + nWake;   % shift sleep block into combined coordinates
pEndS   = pEndS + nWake;

for i = 1:numel(pStartW)
    rectangle(ax3, 'Position',[pStartW(i) pStartW(i) pEndW(i)-pStartW(i) pEndW(i)-pStartW(i)], ...
        'EdgeColor', cMap(i,:), 'LineWidth', 3);
end
for i = 1:numel(pStartS)
    rectangle(ax3, 'Position',[pStartS(i) pStartS(i) pEndS(i)-pStartS(i) pEndS(i)-pStartS(i)], ...
        'EdgeColor', cMap(i,:), 'LineWidth', 3);
end

hold(ax3, 'off');

% save figure:
fileName = [analysisFolderWake filesep 'SleepWakeCorrelogram'];
print(fileName,'-dpdf','-r300');

% --- helper: replicates DendrogramMatrix's internal boundary logic ---
function [pStart, pEnd] = getClusterBoundaries(clusters, order)
    ordClu = clusters(order);
    pStart = [1, find(diff(ordClu)==1)'+1];
    pEnd   = [pStart(2:end)-1, numel(clusters)];
end


%% panels  C, D, E, F
%% Sleep and wake comparison - only on recs where animals were awake
clearvars
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

%% get the AC, and states for each

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

WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
load([analysisFolderWake filesep 'wakeSleepData.mat'])

%% plot AC
figure;
violinplot(pData,"type",'P2V');
hold on;
swarmchart(pData.type, pData.P2V, 15, 'k', 'filled', 'XJitterWidth', 0.3);
hold off;

% save figure:
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

% %% plot bouts - sleep vs wake
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
violinHalfWidth = 0.4;    % max half-width of the density outline
jitterHalfWidth = 0.35;   % max half-width for scattered points

figure; hold on;

for i = 1:nGroups
    grpData = subBouts.duration_s(subBouts.group == groupList{i});
    if isempty(grpData)
        continue
    end

    % --- violin-style outline: KDE of duration, mirrored around baseline x = i ---
    [f, xi] = ksdensity(grpData);
    f = f / max(f) * violinHalfWidth;
    xLeft  = i - f;
    xRight = i + f;
    fill([xRight, fliplr(xLeft)], [xi, fliplr(xi)], colorMat(i,:), ...
        'FaceAlpha', 0.20, 'EdgeColor', colorMat(i,:), 'LineWidth', 1.3);

    % --- individual points, jittered in x, width scaled by local density ---
    fPoint = interp1(xi, f, grpData, 'linear', 0);
    jitterWidth = min(fPoint, jitterHalfWidth);
    xJitter = i + (rand(size(grpData)) - 0.5) .* 2 .* jitterWidth;
    scatter(xJitter, grpData, 10, colorMat(i,:), 'filled', 'MarkerFaceAlpha', 0.6);
end

xticks(1:nGroups);
xticklabels(strrep(groupList, '_', ' '));
xlim([0.5, nGroups + 0.5]);
ylabel('Bout duration (s)');
xlabel('State');
title(sprintf('Bout duration distribution -- %s', thisAnimal));
hold off;

% save figure:
fileName = [analysisFolderWake filesep 'BoutDurationoneAnimalExample'];
print(fileName,'-dpdf','-r300');


%% ---------- Figure 2: horizontal half-violin (raincloud-style), by TYPE, overlaid ----------
allBouts.type = categorical(allBouts.type);   % ensure categorical (Sleep / Wake)
subBouts.type = categorical(subBouts.type);
typeList = categories(allBouts.type);
nTypes   = numel(typeList);

typeColorMat = [
    0.20 0.20 0.60;   % Sleep - dark indigo
    0.90 0.55 0.10    % Wake  - orange
];   % edit/reorder if categories(allBouts.type) comes out in a different order

densityHeight = 0.4;    % how tall the density hump rises above baseline
jitterHeight  = 0.4;   % how far points scatter below baseline
baseline = 1;            % single shared row

figure; hold on;
h = gobjects(nTypes, 1);

for i = 1:nTypes
    tData = subBouts.duration_s(subBouts.type == typeList{i});
    if isempty(tData)
        continue
    end

    % --- upper half only: density hump rising from baseline (histogram-like) ---
    [f, xi] = ksdensity(tData);
    f = f / max(f) * densityHeight;
    yTop = baseline + f;
    h(i) = fill([xi, fliplr(xi)], [yTop, baseline*ones(size(xi))], typeColorMat(i,:), ...
        'FaceAlpha', 0.25, 'EdgeColor', typeColorMat(i,:), 'LineWidth', 1.3);

    % --- points scattered below baseline only (raincloud strip) ---
    fPoint = interp1(xi, f, tData, 'linear', 0);
    jitterWidth = min(fPoint, jitterHeight);
    yJitter =  baseline + abs((rand(size(tData)) - 0.5)) .* 2 .* jitterWidth;
    scatter(tData, yJitter, 5, typeColorMat(i,:), 'filled', 'MarkerFaceAlpha', 0.5);
end

yticks([]);
ylim([baseline - 0.01, baseline + densityHeight + 0.01]);
xlabel('Bout duration (s)');
title(sprintf('Bout duration distribution (overlaid, by type) -- %s', thisAnimal));
legend(h, strrep(typeList, '_', ' '), 'Location', 'best');
hold off;

% save figure:
fileName = [analysisFolderWake filesep 'BoutDurationoneAnimalExampleHorizontalSW'];
print(fileName,'-dpdf','-r300');


%% ---------- Setup ----------
% Assumes 'allBouts' (with .Animal, .type, .label, .group, .duration_s) and
% 'thisAnimal' are already defined, as in the earlier scripts.
thisAnimal = "PV157";
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
figure;

subplot(2,1,1); hold on;
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

subplot(2,1,2); hold on;
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
legend(hLow, strrep(lowGroups, '_', ' '), 'Location', 'best');

sgtitle(sprintf('Bout duration distribution -- %s', thisAnimal));

% save figure
fileName = [analysisFolderWake filesep 'BoutDurationoneAnimalExampleHorizontalSW_HighLow'];
print(fileName, '-dpdf', '-r300');

%% ---------- Figure 2: per-animal Sleep-minus-Wake density difference, High / Low subplots ----------
animalList = categories(categorical(allBouts.Animal));
nAnimals   = numel(animalList);
nGrid      = 200;

sleepHighName = groupList{find(contains(groupList,'Sleep','IgnoreCase',true) & contains(groupList,'High','IgnoreCase',true), 1)};
wakeHighName  = groupList{find(contains(groupList,'Wake', 'IgnoreCase',true) & contains(groupList,'High','IgnoreCase',true), 1)};
sleepLowName  = groupList{find(contains(groupList,'Sleep','IgnoreCase',true) & contains(groupList,'Low', 'IgnoreCase',true), 1)};
wakeLowName   = groupList{find(contains(groupList,'Wake', 'IgnoreCase',true) & contains(groupList,'Low', 'IgnoreCase',true), 1)};

figure;

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

    fPoint = interp1(xi, f, data, 'linear', 0);
    jitterWidth = min(fPoint, jitterHeight);
    yJitter = baseline + abs((rand(size(data)) - 0.5)) .* 2 .* jitterWidth;
    scatter(data, yJitter, 5, color, 'filled', 'MarkerFaceAlpha', 0.5);
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


figure;
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
fileName2 = [analysisFolderWake filesep 'BoutsPercentagesAllAnimals'];
print(fileName2, '-dpdf', '-r300');
%%  plot mean DB for each state
groupOrder = {'Sleep','Wake'};
pData.type = categorical(pData.type, groupOrder);

% Reshape to long format
longType  = [pData.type; pData.type];
longValue = [pData.HighDBAvg; pData.LowDBavg];
longLabel = [repmat({'highDB'}, height(pData), 1); repmat({'lowDB'}, height(pData), 1)];
longLabel = categorical(longLabel, {'highDB','lowDB'});

% Base x position from type, with small offset per sub-group
baseX = double(longType);
offset = 0.30;
xPos = baseX - offset .* (longLabel == 'highDB') + offset .* (longLabel == 'lowDB');

colorMat = [
    0.85 0.33 0.10;   % Sleep_HighDB  - orange (dark)
    0.00 0.45 0.74;   % Sleep_LowDB   - blue (dark)
    0.93 0.60 0.40;   % Wake_HighDB   - orange (light)
    0.40 0.70 0.90    % Wake_LowDB    - blue (light)
];

% Per-point color matrix
pointColors = zeros(numel(xPos), 3);
pointColors(longType=='Sleep' & longLabel=='highDB', :) = repmat(colorMat(1,:), sum(longType=='Sleep' & longLabel=='highDB'), 1);
pointColors(longType=='Sleep' & longLabel=='lowDB',  :) = repmat(colorMat(2,:), sum(longType=='Sleep' & longLabel=='lowDB'),  1);
pointColors(longType=='Wake'  & longLabel=='highDB', :) = repmat(colorMat(3,:), sum(longType=='Wake'  & longLabel=='highDB'), 1);
pointColors(longType=='Wake'  & longLabel=='lowDB',  :) = repmat(colorMat(4,:), sum(longType=='Wake'  & longLabel=='lowDB'),  1);

figure; hold on;
swarmchart(xPos, longValue, 25, pointColors, 'filled', 'o', 'XJitterWidth', 0.25);
% Xticks at each sub-group position, labeled directly
xticks([1-offset, 1+offset, 2-offset, 2+offset]);
xticklabels({'Sleep-HighDB','Sleep-LowDB','Wake-HighDB','Wake-LowDB'});
xtickangle(0);

xlim([0.5, numel(groupOrder)+0.5]);
ylabel('Avg DB Ratio');


%% Panels ??? - PCA 
clearvars
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
%% parameters for freqbands
% data parameters:
fMax = 80;
segmentLength = 1000;
binDuration = 5000;
WelchOL = 0.75;
dftPoints = 2^10;
% win = 0.5*60*60*1000;

%% get wake subset - all trials 
wakeTrials = WA.recTable.ForPaper==1 & cellfun(@(x) contains(x, 'Wake'), WA.recTable.recNames);
wakeTable = WA.recTable(wakeTrials,:);
% [G, AnimalsWake] = findgroups(wakeTable.Animal);
% idx = splitapply(@(x){x(1:min(4,numel(x)))}, ...
%                  (1:height(wakeTable))', G);
% 
% idx = vertcat(idx{:});
% % trialsNames = ["Wake1","Wake2","Wake3","Wake4"];
% % trialsid  = ismember(wakeTable.recNames,trialsNames);
% wakeSubsetAll = wakeTable(idx ,:);
% running on a single trial, power spec
PxxDataAllWake = [];
nPxxDataAllWake = [];
timesAllWake   = [];
freqHzAllWake  = [];
clustersWake =[];
recLabelsWake  = [];  % track which recording each row came from
animalLabelWake = [];

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
    end

    animal = wakeTable.Animal(i);
    
    curfreqData = WA.getFreqBandDetection(...
        'fMax',          fMax, ...
        'segmentLength', segmentLength, ...
        'binDuration',   binDuration, ...
        'WelchOL',       WelchOL, ...
        'dftPoints',     dftPoints, ...
        'tStart',   tStart, ...
        'win',          win,...
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
    recLabelsWake  = [recLabelsWake;  repmat(wakeTable.recNames(i), nBins, 1)];             % which rec each row is from
    animalLabelWake = [animalLabelWake; repmat(animal, nBins, 1)];
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'WakePCAdata.mat'], 'PxxDataAllWake','nPxxDataAllWake', 'freqHzAllWake','clustersWake', 'timesAllWake', ...
    'recLabelsWake', 'animalLabelWake');


%% get sleep subset
sleepTrials = WA.recTable.ForPaper==1& contains(WA.recTable.recNames,'Night');
sleepTable = WA.recTable(sleepTrials,:);

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
for i = 1:height(sleepTable)
    recName = ['Animal=' sleepTable.Animal{i} ',recNames=' sleepTable.recNames{i}];
    WA.setCurrentRecording(recName);
    animal = sleepTable.Animal(i);
    % get the sleep timings: 
    WA.getDelta2BetaRatio;
    % WA.plotDelta2BetaRatio;
    curAC = WA.getDelta2BetaAC;
    tStartSleep = curAC.tStartSleep+0.5*60*60*1000;
    % tEndSleep = curAC.tEndSleep;
    win= 1*60*60*1000;
    
    curfreqData = WA.getFreqBandDetection(...
        'fMax',          fMax, ...
        'segmentLength', segmentLength, ...
        'binDuration',   binDuration, ...
        'WelchOL',       WelchOL, ...
        'dftPoints',     dftPoints, ...
        'tStart',       tStartSleep,...
        'win',          win, ...
        'overwrite', 1 ...
            );
        
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
    recLabelsSleep  = [recLabelsSleep;  repmat(sleepTable.recNames(i), nBins, 1)];             % which rec each row is from
    animalLabelSleep = [animalLabelSleep; repmat(animal, nBins, 1)];
end

% Save everything needed for PCA
save([analysisFolderWake filesep 'SleepPCAdata.mat'], 'PxxDataAllSleep','nPxxDataAllSleep', 'freqHzAllSleep','clustersSleep', ...
    'clustersSleepAno','timesAllSleep', 'recLabelsSleep', 'animalLabelSleep');

%% load data and remove freq 0
% WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
% analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

load([analysisFolderWake filesep 'SleepPCAdata.mat'])
% load([analysisFolderWake filesep 'allSpectraHunter.mat'])
% load([analysisFolderWake filesep 'allSpectraHunterBugs.mat'])
load([analysisFolderWake filesep 'WakePCAdata.mat'])

freqHzAll = freqHzAllSleep;
% remove 0 freq
f_ind=freqHzAll>0.5;
freqHzAll = freqHzAll(f_ind);
% freqHzAllHun = freqHzAllHun(f_ind);
freqHzAllWake= freqHzAllWake(f_ind);
nPxxDataAllWake=nPxxDataAllWake(:,f_ind);
% nPxxDataAllHun=nPxxDataAllHun(:,f_ind);
nPxxDataAllSleep = nPxxDataAllSleep(:,f_ind);

% PxxDataAllAfter=PxxDataAllAfter(:,f_ind);
% PxxDataAllBefore=PxxDataAllBefore(:,f_ind);
PxxDataAllSleep=PxxDataAllSleep(:,f_ind);
PxxDataAllWake=PxxDataAllWake(:,f_ind);
% PxxDataAllHun=PxxDataAllHun(:,f_ind);


%% PCA with sleep; no pre-post:
%% 1) Each animal's own Wake PCA  2) All-animals Wake PCA

% Parameters for each PCA type
dx1_raw = 500;  dx2_raw = 500;
dx1_nt  = 1;    dx2_nt  = 1;
dx1_ns  = 1;    dx2_ns  = 1;

pcaXlim_raw = [-0.5 4.5]*10000; pcaYlim_raw = [-1.5 3]*10000;
pcaXlim_nt  = [-40 120];        pcaYlim_nt  = [-50 50];
pcaXlim_ns  = [-40 40];         pcaYlim_ns  = [-30 30];

% if all(freqHzAllWake == freqHzAllHun)
%     freqHzAll = freqHzAllWake;
% end

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

    % maskHun   = strcmp(animalLabelHun,   thisAnimal);
    maskWake  = strcmp(animalLabelWake,  thisAnimal);
    maskSleep = strcmp(animalLabelSleep, thisAnimal);

    % Raw data for this animal
    % pxxHun   = PxxDataAllHun(maskHun,:);
    pxxWake  = PxxDataAllWake(maskWake,:);
    pxxSleep = PxxDataAllSleep(maskSleep,:);

    % nPxxHun   = nPxxDataAllHun(maskHun,:);
    nPxxWake  = nPxxDataAllWake(maskWake,:);
    nPxxSleep = nPxxDataAllSleep(maskSleep,:);

    % pxxHun_norm   = normalize(pxxHun,   2);
    pxxWake_norm  = normalize(pxxWake,  2);
    pxxSleep_norm = normalize(pxxSleep, 2);

    % Clusters
    % clustHun   = clustersHun(maskHun);
    clustWake  = clustersWake(maskWake);
    clustSleep = clustersSleepAno(maskSleep);

    % nClustersWakeHun = max([clustHun; clustWake]);
    nClustersWake = 2;
    nClustersSleep   = 2;
    clustersNames    = {"REM","SWS"};

    clusterColorsWake = lines(nClustersWake);
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
    % scoresA_hun_own    = (pxxHun   - meanA_own) * coeffA_own;
    scoresA_sleep_own  = (pxxSleep - meanA_own) * coeffA_own;

    scoresAnt_wake_own  = (nPxxWake  - meanAnt_own) * coeffAnt_own;
    % scoresAnt_hun_own   = (nPxxHun   - meanAnt_own) * coeffAnt_own;
    scoresAnt_sleep_own = (nPxxSleep - meanAnt_own) * coeffAnt_own;

    scoresAns_wake_own  = (pxxWake_norm  - meanAns_own) * coeffAns_own;
    % scoresAns_hun_own   = (pxxHun_norm   - meanAns_own) * coeffAns_own;
    scoresAns_sleep_own = (pxxSleep_norm - meanAns_own) * coeffAns_own;

    pcaData_own = {
        coeffA_own,   scoresA_wake_own,   scoresA_sleep_own,   explainedA_own,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt_own, scoresAnt_wake_own, scoresAnt_sleep_own, explainedAnt_own, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns_own, scoresAns_wake_own, scoresAns_sleep_own, explainedAns_own, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['[Own PCA] Animal: ' thisAnimal], ...
           'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal ' � PCA from own Wake data']);
    plotPCAFigure(pcaData_own, freqHzAll, clustWake, clustSleep, ...
                  nClustersWake, nClustersSleep, clustersNames, ...
                  clusterColorsWake, clusterColorsSleep);

    % ----------------------------------------------------------
    % APPROACH 2: Project onto global (all-animals) Wake PCA
    % ----------------------------------------------------------
    scoresA_wake_g   = (pxxWake  - meanA_global) * coeffA_global;
    % scoresA_hun_g    = (pxxHun   - meanA_global) * coeffA_global;
    scoresA_sleep_g  = (pxxSleep - meanA_global) * coeffA_global;

    scoresAnt_wake_g  = (nPxxWake  - meanAnt_global) * coeffAnt_global;
    % scoresAnt_hun_g   = (nPxxHun   - meanAnt_global) * coeffAnt_global;
    scoresAnt_sleep_g = (nPxxSleep - meanAnt_global) * coeffAnt_global;

    scoresAns_wake_g  = (pxxWake_norm  - meanAns_global) * coeffAns_global;
    % scoresAns_hun_g   = (pxxHun_norm   - meanAns_global) * coeffAns_global;
    scoresAns_sleep_g = (pxxSleep_norm - meanAns_global) * coeffAns_global;

    pcaData_global = {
        coeffA_global,   scoresA_wake_g,      scoresA_sleep_g,   explainedA_global,   'Raw',              dx1_raw, dx2_raw, pcaXlim_raw, pcaYlim_raw;
        coeffAnt_global, scoresAnt_wake_g, scoresAnt_sleep_g, explainedAnt_global, 'Norm. to trial',   dx1_nt,  dx2_nt,  pcaXlim_nt,  pcaYlim_nt;
        coeffAns_global, scoresAns_wake_g, scoresAns_sleep_g, explainedAns_global, 'Norm. to segment', dx1_ns,  dx2_ns,  pcaXlim_ns,  pcaYlim_ns;
    };

    figure('Name', ['[Global PCA] Animal: ' thisAnimal], ...
           'Units','normalized', 'Position',[0 0 1 1]);
    sgtitle(['Animal: ' thisAnimal ' � projected onto all-animals Wake PCA']);
    plotPCAFigure(pcaData_global, freqHzAll, clustWake, clustSleep, ...
                  nClustersWake, nClustersSleep, clustersNames, ...
                  clusterColorsWake, clusterColorsSleep);
end
