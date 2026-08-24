%% Figure 2 - Wake vs Sleep figrue
%% get the data
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';

%% Panel A - examples of sleep and wake


%% Panel B - correlogram both sleep and wake
% wake - PV157,Wake3,hour 0-1
% sleep - PV157, Night1, hour 6-7

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

WA.setCurrentRecording('Animal=PV157,recNames=Night2');
sleepFreqData = WA.getFreqBandDetection(...
    'fMax',          fMax, ...
    'segmentLength', segmentLength, ...
    'binDuration',   binDuration, ...
    'WelchOL',       WelchOL, ...
    'dftPoints',     dftPoints, ...
    'overwrite',     1, ...
    'tStart',        7*60*60*1000, ...
    'win', win ...
    );

    WA.plotFreqBandDetection;

% plot them together with the dendro tree:

f = figure('Position',[100 100 798 601]);

% Outer margin around the whole layout
margin = 0.05;

% Reserve exact space so nothing overflows
leftWidth  = 0.18;
cbarWidth  = 0.025;
cbarGap    = 0.015;

usableWidth  = 1 - 2*margin;
usableHeight = 1 - 2*margin;
ax3Width = usableWidth - leftWidth - cbarWidth - cbarGap;

ax1 = axes('Parent', f, 'Position', [margin, margin + usableHeight/2, leftWidth, usableHeight/2]);
ax2 = axes('Parent', f, 'Position', [margin, margin,                 leftWidth, usableHeight/2]);
ax3 = axes('Parent', f, 'Position', [margin + leftWidth, margin, ax3Width, usableHeight]);

% Lock positions so axis square / dendrogram plotting can't move them
set(ax1, 'PositionConstraint', 'innerposition');
set(ax2, 'PositionConstraint', 'innerposition');
set(ax3, 'PositionConstraint', 'innerposition');

maxDendroClusters = 2;

% wake tree into ax1
wakeC = corrcoef(wakeFreqData.normsPxx);
[~, orderW, clustersW] = DendrogramMatrix(wakeC,'linkMetric','euclidean','linkMethod','ward', ...
                        'maxClusters',maxDendroClusters, 'figureHandle',f, ...
                        'toPlotBinaryTree',1,'hDendro',[ax1,ax3],'plotOrderLabels',0);

% sleep tree into ax2
sleepC = corrcoef(sleepFreqData.normsPxx);
[~, orderS, clusterS] = DendrogramMatrix(sleepC,'linkMetric','euclidean','linkMethod','ward', ...
                        'maxClusters',maxDendroClusters, 'figureHandle',f, ...
                        'toPlotBinaryTree',1,'hDendro',[ax2,ax3],'plotOrderLabels',0);

% Combined ordered matrix into ax3
normsPXXboth = [wakeFreqData.normsPxx, sleepFreqData.normsPxx];
corrMatboth = corrcoef(normsPXXboth);
orderBoth = [orderW, orderS + length(orderW)];
reordered_corr = corrMatboth(orderBoth, orderBoth);

delete(ax3)
ax3 = axes('Parent', f, 'Position', [margin + leftWidth, margin, ax3Width, usableHeight]);
set(ax3, 'PositionConstraint', 'innerposition');
axes(ax3)
imagesc(reordered_corr);
axis square;
colormap(parula);
clim([-1 1]);

% Explicitly position the colorbar in the reserved slice � won't overflow
cb = colorbar(ax3);
cb.Position = [margin + leftWidth + ax3Width + cbarGap*0.3, margin + 0.15*usableHeight, cbarWidth, 0.7*usableHeight];

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
load([analysisFolderWake filesep 'wakeSleepData.mat'])

%% pData.animalC = categorical(pData.Animal);
% plot AC
figure;
violinplot(pData,"type",'P2V');
hold on;
swarmchart(pData.type, pData.P2V, 15, 'k', 'filled', 'XJitterWidth', 0.3);
hold off;

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

%% plot bouts - sleep vs wake


allBouts.group = categorical(strcat(string(allBouts.type), "-", string(allBouts.label)));

% figure;
% violinplot(allBouts.group,allBouts.duration_s);
% % legend(unique(allBouts.Animal),'Location','best');
% ylabel('Bout duration (s)');
% xlabel('State');
% title('Bout duration distribution');
% %%
groupList = categories(allBouts.group);  % should be: Sleep_HighDB, Sleep_LowDB, Wake_HighDB, Wake_LowDB

% base colors: blue-family for HighDB (1,3), orange-family for LowDB (2,4)
colorMat = [
    0.85 0.33 0.10;   % 1: Sleep_HighDB  - orange (dark)
    0.00 0.45 0.74;   % 2: Sleep_LowDB   - blue (dark)
    0.93 0.60 0.40;   % 3: Wake_HighDB   - orange (light)
    0.40 0.70 0.90    % 4: Wake_LowDB    - blue (light)
];

[~, idx] = ismember(allBouts.group, groupList);
pointColors = colorMat(idx, :);

figure;
swarmchart(allBouts.group, allBouts.duration_s, 5, pointColors, 'filled');
ylabel('Bout duration (s)');
xlabel('State');
title('Bout duration distribution');

% per animals
allBouts.Animal = categorical(allBouts.Animal);
animalList = categories(allBouts.Animal);
groupList = categories(allBouts.group);
cmap = lines(numel(animalList));

nAnimals = numel(animalList);
offsets = linspace(-0.3, 0.3, nAnimals);  % spread animals within each group's x-slot

figure;
hold on;
h = gobjects(nAnimals,1);
for i = 1:nAnimals
    idx = allBouts.Animal == animalList{i};
    
    % numeric x position = group code + this animal's offset
    xNum = double(allBouts.group(idx)) + offsets(i);
    
    h(i) = swarmchart(xNum, allBouts.duration_s(idx), 10, cmap(i,:), 'filled', ...
        'XJitterWidth', 0.15);
end

% set x-ticks to group names at the center of each cluster
xticks(1:numel(groupList));
xticklabels(groupList);
xlim([0.5, numel(groupList)+0.5]);

legend(h, animalList, 'Location', 'best');
ylabel('Bout duration (s)');
xlabel('State');
title('Bout duration distribution, per animal');allBouts.Animal = categorical(allBouts.Animal);
animalList = categories(allBouts.Animal);
groupList = categories(allBouts.group);
cmap = lines(numel(animalList));

%% check Standart deviation for each state
stdgroups= [];
for i = 1:length(groupList)
    group = groupList{i};
    curstd = std(allBouts.duration_s(allBouts.group==group));
    stdgroups= [stdgroups;curstd];
end
STDgruopT = table(groupList,stdgroups,'VariableNames',["Group","STD"]);

%% check long bouts
l = allBouts.duration_s>500;

longBouts=allBouts(l,:);

[found, o] = ismember(longBouts(:, {'Animal','recNames'}), ...
                         pData(:, {'Animal','recNames'}), 'rows');
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
