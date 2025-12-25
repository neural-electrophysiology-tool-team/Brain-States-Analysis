%% Supplemental Materials
% initiating stimTable and other essentials. 

SA=sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTableAll.mat'])
load([analysisFolder filesep 'LMDataAll.mat'])
animalsColors = [
    0/255,124/255,143/255; % Hex: 00798C- dark blue  - PV106
    255/255, 142/255, 71/255;% HEX:  FF8E47 - orange  - PV126
    % 121/255, 67/255, 157/255;% HEX: 79439D- dark purple -PV142
    156/255,44/255,58/255; % HEX: 9C2C3A red-ish - PV143
    28/255, 144/255, 217/255;  % HEX: 1C90D9 - blue - PV149
    237/255,174/255,73/255; % HEX:edae49 - yellow-ish - PV153
    % 77/255, 147/255, 57/255; %HEX: 4D9339 dark green- PV88
    148/255, 120/255, 186/255; % HEX: 9478BA - purple - PV157
    217/255, 62/255, 67/255; % HEX: D93E43 - red - PV159
    255/255, 202/255, 58/255; % HEX: FFCA3A - yellow -  PV161
    97/255, 184/255, 68/255;  % HEX:61B844 - Green -PV162
    ];

uniqueAnimals = unique(stimTable.Animal);

%% Supplementary Figure 1
%% all stimualtions raw data one night.
SA.setCurrentRecording('Animal=PV153,recNames=Night15');
stims = SA.getStimTriggers;
ch = SA.recTable.defaulLFPCh(SA.currentPRec);
pre = 20*1000;
tStart = stims(1)-pre; win = stims(end)-tStart+50*1000;
[LFP,LFP_t] = SA.currentDataObj.getData(ch,tStart,win);
LFP = squeeze(LFP);
dsLFP = downsample(LFP,200);
dsLFP_t = downsample(LFP_t,200);

%% plot
chunkDur_s = 158.2;                 % seconds per row
chunkDur_ms = chunkDur_s * 1000;  % ms
gap = 1200;                       % vertical spacing

t = dsLFP_t;   % time in ms
y = dsLFP;

dt = median(diff(t));                 % ms
chunkSamples = round(chunkDur_ms / dt);
nChunks = ceil(length(y) / chunkSamples);

% height of stim marker within each row (in signal units)
stimHalfHeight = gap * 0.5;          % adjust if you want taller/shorter lines
%
f=figure; hold on

for k = 1:nChunks
    idx1 = (k-1)*chunkSamples + 1;
    idx2 = min(k*chunkSamples, length(y));

    % time within chunk (ms)
    t0 = t(idx1);
    tChunk = t(idx1:idx2) - t0;

    % vertical offset (top -> bottom)
    yOffset = (nChunks - k) * gap;

    % plot chunk
    yChunk = y(idx1:idx2) + yOffset;
    plot(tChunk/1000, yChunk, 'k')  % x-axis in seconds

    % ---- stim markers ONLY for this row ----
    stimRel = (stims - tStart);   % ms in same reference as t
    stimIdx = stimRel >= t(idx1) & stimRel <= t(idx2);

    if any(stimIdx)
        x = (stimRel(stimIdx) - t0) / 1000;   % seconds within chunk

        % draw short vertical segments centered around this row's offset
        y0 = yOffset - stimHalfHeight;
        y1 = yOffset + stimHalfHeight;

        for xi = x(:)'
            line([xi xi], [y0 y1], 'Color', 'r', 'LineWidth', 1.5);
        end
    end
end

xlabel('Time (s)')
set(gca,'YTick',[]); xlim([0,chunkDur_s])
box on

% save figue
set(f,'PaperPosition', [1 1 7 9]);
% set(f,"PaperPositionMode","auto")
fileName=[analysisFolder filesep 'AllNightRawLFPPV153N12'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Supplementary Figure 1
%% Sup 1 C, D,E
% Traces of stimulations with virus:

i = 6;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;

% stimulations timings:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch};
stimStartT = stimTable.stimStartT(i);
stimEndT = stimTable.stimEndT(i);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+200;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=100000;
ch = 17;
% get the raw data and LFP:
% for j = 10:40
 j = 39;
 VipTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
 VidbTmp=DB.bufferedDelta2BetaRatio(VipTmp);
 [vlfp,vlfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    

  f = figure;
    h1 = subplot(2,1,1);
    % set(f, 'Position', [100, 100, 1200, 400]);
    % yyaxis left
    plot(vlfp_t/1000,squeeze(vlfp),'k'); hold on;
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims/1000,'r','LineWidth',1.5)
    % yyaxis right
    plot(VidbTmp,'Color','b','LineWidth',2)
    title('With virus');
    % sgtitle(sprintf('Trial num: %i',j))
    hold off;
% end

% Traces of stimulations without virus:

i = 9;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;

% stimulations timings:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch};
stimStartT = stimTable.stimStartT(i);
stimEndT = stimTable.stimEndT(i);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+400;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=100000;
ch = 17;

% get the raw data and LFP:
% for j = 10:40
 j = 28;
 pTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
 dbTmp=DB.bufferedDelta2BetaRatio(pTmp);
 [lfp,lfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    
  % f = figure;
    h2 = subplot(2,1,2);
    % set(f, 'Position', [100, 100, 1200, 400]);
    % yyaxis left
    plot(lfp_t/1000,squeeze(lfp),'k'); hold on;
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims/1000,'r','LineWidth',1.5)
    % yyaxis right
    plot(dbTmp,'Color','b','LineWidth',2)
    title('Without Virus');
    % sgtitle(sprintf('W/O Trial num:%i',j))
    hold off;
%end


%save figue
set(f,'PaperPosition', [1 1 4 2.5]);
fileName=[analysisFolder filesep 'tracesVirusNovirus'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% Sup 1F

viAnimals = ["161","149"];
curType = "47";

f=figure;

%plot the data
% plot nights without :
nonViTrials = contains(stimTable.Remarks,curType) &...
    ~contains(stimTable.Remarks,'Ex') & ...
    any(~contains(stimTable.Animal,viAnimals),2);
nNon = sum(nonViTrials);
Nnon= length(unique(stimTable.Animal(nonViTrials)));
ViTrials = contains(stimTable.Remarks,curType) &...
    ~contains(stimTable.Remarks,'Ex') & ...
    any(contains(stimTable.Animal,viAnimals),2);
nVi = sum(ViTrials);
Nvi = length(unique(stimTable.Animal(ViTrials)));
 
% Statistics:
[pMW_sham,h_sham] = ranksum(stimTable.dbDiffShamM(ViTrials), stimTable.dbDiffShamM(nonViTrials));
    fprintf('Mann-Whitney Signed-Rank Test p-value - sham: %.4f\n', pMW_sham);
    fprintf('h = %i\n',h_sham)

[pMW_stim, h_stim] = ranksum(stimTable.dbDiffStimM(ViTrials), stimTable.dbDiffStimM(nonViTrials));
    fprintf('Mann-Whitney Signed-Rank Test p-value - stim: %.4f\n', pMW_stim);
    fprintf('h = %i\n',h_stim)

%plot
x = 1:2;
h1 = plot(x,[stimTable.dbDiffShamM(nonViTrials) stimTable.dbDiffStimM(nonViTrials)], ...
    'Marker','.','LineStyle','-','Color',[1 0.65 0]);
hold on
h2 = plot(x,[stimTable.dbDiffShamM(ViTrials) stimTable.dbDiffStimM(ViTrials)], ...
    'Marker','.','LineStyle','-','Color','g');
xlim([0.9 2.1])
xticks([x]);xticklabels(["Sham","Stim"])

legend([h1(1), h2(1)],"Animals not injected","Animals injected with ChR2")


% save figue

set(f,'PaperPosition',[1 1 2 3]);
fileName=[analysisFolder filesep 'virusChangeStimShamBlue'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);





%% Supplementary Figure 2 - SWR
%% Sup 2A
i = 63;%(17 before). not PV106,N20
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
stims = SA.getStimTriggers; %stimulations timimng in ms
firstTrig = stims(1:8:end);

j=2;
pre = 5000;
ch = SA.recTable.defaulLFPCh(SA.currentPRec);
tStart = firstTrig(j)-pre;
win = pre+40*1000;
[rawLFP,rawLFP_t] = SA.currentDataObj.getData(ch,tStart,win);
%plot
figure;
plot((rawLFP_t-pre)/1000,squeeze(rawLFP),'k');hold on;
xline((stims(8*j-7:8*j)-stims(8*j-7))/1000,'r');
%save
set(gcf ,'PaperPosition',[1 2 4.3 1.7]);
fileName=[analysisFolder filesep 'SWRraw_rev'];
% print(fileName,'-depsc','-vector');
print(fileName, '-dpdf', '-r300');

%% Sup 2B
% plot all raw traces, 1 nights:

i = 63;%(17 before). not PV106,N20
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
    
% get the stimualtions time stamps:
% t_ch = stimTable.StimTrighCh(i);
% T=SA.getDigitalTriggers;
stims = SA.getStimTriggers; %stimulations timimng in ms
    
% get the raw data for the second around each stim:
pre = 200; 
post = 1000;
[mV, mV_t] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),stims-pre,(pre+post));
% [mAn, mAn_t] = SA.currentDataObj.getAnalogData(SA.recTable.stimDiodeCh(SA.currentPRec),stims-pre,(pre+post));

% plot
fraw = figure;
x = linspace(-pre,post,length(mV_t));
mV1 = squeeze(mV);
[d,h] = hist2(repmat(x,height(mV1),1),mV1,'dX1',3,'dX2',5);
hold on
plot(downsample(x,100),downsample(mean(mV1,1),100),'Color','r','LineWidth',2)
ylims =[-600 400];
ylim(ylims)
% title('Raw traces for 1 night')
x_shade = [0 200 200 0];  % X-coordinates of the shaded region
y_shade = [ylims(1) ylims(1) ylims(2) ylims(2)]; % Y-coordinates covering the full y-range
patch(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('Time(ms)');ylabel('Voltage[uV]')

% subplot(2,1,2)
% x = linspace(-pre,post,length(mV_t));
% mAn1 = squeeze(mAn);
% [d,h] = hist2(repmat(x,height(mAn1),1),mAn1,'dX1',3,'dX2',5);
% hold on
% plot(downsample(x,100),downsample(mean(mAn1,1),100),'Color','r','LineWidth',2)
% ylims =[-1000 700];
% ylim(ylims)
% % title('Raw traces for 1 night')
% x_shade = [0 200 200 0];  % X-coordinates of the shaded region
% y_shade = [ylims(1) ylims(1) ylims(2) ylims(2)]; % Y-coordinates covering the full y-range
% patch(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% xlabel('Time(ms)');ylabel('Voltage[uV]')
% 

% save figure
set(fraw ,'PaperPosition',[1 2 4.3 1.7]);
fileName=[analysisFolder filesep 'SWRhist_rev'];
% print(fileName,'-depsc','-vector');
print(fileName, '-dpdf', '-r300');

%% Sup 2C+D
% plot cros corr 1 night:
% detect SWR using SWR detection method.

% for 1 night first:
i = 63;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
SA.getSharpWaves('detectOnlyDuringSWS',false)
SW = SA.getSharpWaves('detectOnlyDuringSWS',false);
swrT = SW.tSW'; %SWR timings in ms
swrN = length(swrT); 



% get the stimualtions time stamps:
% t_ch = stimTable.StimTrighCh(i);
% T=SA.getDigitalTriggers;
% stims = T.tTrig{t_ch}; %stimulations timimng in ms
stims= SA.getStimTriggers;
first_stims = stims(1:8:end-2);
% cross cor: Short, one night 
% Convert Timestamps to Binary Time Series:
bin_size = 200; % ms 
winShort = 10*1000; % ten seconds

ylims = [-0.05 0.1];
t_min = min([stims swrT]);
t_max = max([stims swrT]);
time_bins = t_min:bin_size:t_max;

% Convert timestamps to binary spike trains
swrBin = histcounts(swrT, [time_bins, time_bins(end)+bin_size]);
stimsBin = histcounts(stims, [time_bins, time_bins(end)+bin_size]);
max_lag = winShort/bin_size;
[c, lags,bounds] = crosscorr(stimsBin,swrBin,"NumLags",max_lag); % Normalized cross-correlation. remove 'coeff' for not-normalized

% plt:
lag_times = (lags * bin_size)/1000; % Convert lag indices to time (and moving from ms to s
f = figure;
subplot (1,2,1)
plot(lag_times, c,'Color','k','Marker','.');
yline(bounds,'Color',[0.7 0.7 0.7],'LineStyle','--')
yline(0,'Color',[0.7 0.7 0.7])
xlabel(['Lag (s)']);
ylim(ylims)
ylabel('Cross-correlation');
title('Cross-Correlation SWR and stimulations');

% cross cor: Long, one night 
% Convert Timestamps to Binary Time Series:
bin_size = 5000; % ms 
winLong = 200*1000; % ten seconds

t_min = min([stims swrT]);
t_max = max([stims swrT]);
time_bins = t_min:bin_size:t_max;

% Convert timestamps to binary spike trains
swrBin = histcounts(swrT, [time_bins, time_bins(end)+bin_size]);
% stimsBin = histcounts(stims, [time_bins, time_bins(end)+bin_size]);
firstStimsBin = histcounts(first_stims, [time_bins, time_bins(end)+bin_size]);
max_lag = winLong/bin_size;
[c, lags,bounds] = crosscorr(firstStimsBin,swrBin,"NumLags",max_lag); % Normalized cross-correlation. remove 'coeff' for not-normalized

% plt:
lag_times = (lags * bin_size)/1000; % Convert lag indices to time (and moving from ms to s
% figure;
subplot(1,2,2)
plot(lag_times, c,'Color','k','Marker','.');
yline(bounds,'Color',[0.7 0.7 0.7],'LineStyle','--')
yline(0,'Color',[0.7 0.7 0.7])
xlabel(['Lag (s)']);
ylim(ylims)
% ylabel('Cross-correlation');
title('Cross-Correlation SWR and stimulations');

% save figure:
set(f,'PaperPosition',[1 2 7 2]);
fileName=[analysisFolder filesep 'SWRcrosPV161N18'];
print(fileName,'-dpdf','-r300');
 
%% Sup 2E+F
fName = [analysisFolder filesep 'ShWcrossCorrDataAll.mat'];
load(fName)

wavelength = 'white';
curTrials = (contains(stimTable.Remarks,wavelength) |  contains(stimTable.Remarks,"DayTime"))...
            & ~contains(stimTable.Remarks,'Ex');
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
curDataShort = cShort(curTrials,:);
curDataLong = cLong(curTrials,:);

fr = figure;
% plot Short according to animaL
subplot(1,2,1)
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curDataShort)
    plot(lagTimeS, curDataShort(i,:),'Color',curColorMat(i,:))%,'Marker','.','MarkerSize',4)
end
plot(lagTimeS, mean(curDataShort,1)','Color','k','LineWidth',1.5)%,'Marker','.','MarkerSize',5)
yline(mean(boundsShort(curTrials,:)),'Color',[0.4 0.4 0.4],'LineStyle','--')
yline(0,'Color',[0.4 0.4 0.4])
ylims = [-0.02 .12];
ylim(ylims)

% plot Long according to animaL
subplot(1,2,2)
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curDataLong)
    plot(lagTimeL, curDataLong(i,:),'Color',curColorMat(i,:))%,'Marker','.','MarkerSize',4)
end
plot(lagTimeL, mean(curDataLong,1)','Color','k','LineWidth',1.5)%,'Marker','.','MarkerSize',5)
yline(mean(boundsLong(curTrials,:)),'Color',[0.4 0.4 0.4],'LineStyle','--')
yline(0,'Color',[0.4 0.4 0.4])
ylim(ylims)

% save figure:
set(fr,'PaperPosition',[1 2 7 2]);
fileName=[analysisFolder filesep 'SWRcrosRedNights'];
print(fileName,'-dpdf','-r300');

%% Supplementary Figure 3
%% Sup 3A
% plot one night: Beta during SW cycle

i =22 ;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);

dbStimData = {1./stimTable.betaSW(i).Pre 1./stimTable.betaSW(i).Stim 1./stimTable.betaSW(i).Post};

fdbDec = figure;
groupNames = ["Pre", "Stim","Post"];
colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

yswarm = [dbStimData{1} dbStimData{2} dbStimData{3}];
xswarm=[1*ones(1,length(dbStimData{1})), 2*ones(1,length(dbStimData{2})), 3*ones(1,length(dbStimData{3}))];
colorswarm = [repmat([0.2, 0.6, 0.8;], length(dbStimData{1}), 1);  % Red for Array 1
    repmat([0.9, 0.4, 0.3], length(dbStimData{2}), 1);  % Green for Array 2
    repmat([0.5, 0.8, 0.5], length(dbStimData{3}), 1)]; % Blue for Array 3
swarmchart(xswarm,yswarm,10,colorswarm,"filled",'o','XJitterWidth',0.5)
ylabel('1/beta power during SWS')
% ylim([ylims])
set(gca, 'XTick', 1:numel(dbStimData), 'XTickLabel', groupNames);
title('Beta power during SWS - one night')

% statistics:

% Group names
groupNames = {'Pre', 'During', 'After'};

% Combine data into a single vector and create a grouping variable
allData = [dbStimData{1},dbStimData{2},dbStimData{3}]'; % Concatenate all data
groupLabels = [...
    repmat({'Pre'}, numel(dbStimData{1}), 1); ...
    repmat({'During'}, numel(dbStimData{2}), 1); ...
    repmat({'After'}, numel(dbStimData{3}), 1)]';

% Kruskal-Wallis test
[pKruskal, tblKruskal, statsKruskal] = kruskalwallis(allData, groupLabels, 'off');
fprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal);

if pKruskal < 0.05
    disp('Significant differences found in Kruskal-Wallis test. Proceeding with pairwise comparisons...');

    if isstring(groupLabels)
        groupLabels = cellstr(groupLabels);
    end

    % Get unique group names
    [groupNames, ~, groupIdx] = unique(groupLabels);
    numGroups = numel(groupNames);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = allData(groupIdx == i);
            data_j = allData(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [groupNames{i} ' vs ' groupNames{j}];
            raw_pvals(idx,1) = p;
            idx = idx + 1;
        end
    end

    % Bonferroni correction
    corrected_pvals = min(raw_pvals * length(raw_pvals), 1);

    % Display results
    fprintf('\nPairwise Wilcoxon Rank-Sum Test (Mann-Whitney U):\n');
    for i = 1:length(raw_pvals)
        fprintf('%s:\t raw p = %.4f,\t Bonferroni-corrected p = %.4f\n', ...
            comparisons{i}, raw_pvals(i), corrected_pvals(i));
    end
else
    disp('No significant differences found in Kruskal-Wallis test.');
end


%savefigure
set(fdbDec,'PaperPosition',[1 1 2.7 2.3]);
fileName=[analysisFolder filesep 'BetaSWOneNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% Sup 3B
% plot all nights: Beta SW- only red

wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength)...
        & ~contains(stimTable.Remarks,'Ex') ...
        & ~any(isnan(stimTable.betaSWMeans),2); 
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curData = 1./stimTable.betaSWMeans(curTrials,:);


% statistics:

[p, tbl, stats] = friedman(curData, 1,'off'); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
if p<0.05
    before = curData(:,1);
    during = curData(:,2);
    after = curData(:,3);

    % Pairwise Wilcoxon signed-rank tests
    [p_before_during, ~, stats_before_during] = signrank(before, during);
    [p_during_after, ~, stats_during_after] = signrank(during, after);
    [p_after_before, ~, stats_after_before] = signrank(after, before);
    raw_pvals = [p_before_during,p_during_after,p_after_before];
    num_comparisons = length(raw_pvals);

    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Before vs During: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('During vs After: p-value = %.4f\n', corrected_pvals_bonferroni(2));
    fprintf('After vs Before: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
end
    %plot
fdb = figure;
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curData)
    plot(curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10)
end
plot(mean(curData,1),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
xlim([0.5, 3.5])
xticks(1:3)
xticklabels(groupNames)
ylabel('1/Beta means during SWS ')

% savefigure
set(fdb,'PaperPosition',[1 1 2.7 2.3]);
fileName=[analysisFolder filesep 'BetachangeSWS'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);



%% Supplementary Figure 4
%% Sup 4A
% plot one night: D/B decrease full cycle

i =22 ;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);

dbStimData = {stimTable.dbFull(i).Pre stimTable.dbFull(i).Stim stimTable.dbFull(i).Post};

fdbDec = figure;
groupNames = ["Pre", "Stim","Post"];
colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

yswarm = [dbStimData{1} dbStimData{2} dbStimData{3}];
xswarm=[1*ones(1,length(dbStimData{1})), 2*ones(1,length(dbStimData{2})), 3*ones(1,length(dbStimData{3}))];
colorswarm = [repmat([0.2, 0.6, 0.8;], length(dbStimData{1}), 1);  % Red for Array 1
    repmat([0.9, 0.4, 0.3], length(dbStimData{2}), 1);  % Green for Array 2
    repmat([0.5, 0.8, 0.5], length(dbStimData{3}), 1)]; % Blue for Array 3
swarmchart(xswarm,yswarm,10,colorswarm,"filled",'o','XJitterWidth',0.5)
ylabel('D/B power during full cycle')
ylim([0 300])
set(gca, 'XTick', 1:numel(dbStimData), 'XTickLabel', groupNames);
title('D/B power during full cycles - one night')

% statistics:

% Group names
groupNames = {'Pre', 'During', 'After'};

% Combine data into a single vector and create a grouping variable
allData = [dbStimData{1},dbStimData{2},dbStimData{3}]'; % Concatenate all data
groupLabels = [...
    repmat({'Pre'}, numel(dbStimData{1}), 1); ...
    repmat({'During'}, numel(dbStimData{2}), 1); ...
    repmat({'After'}, numel(dbStimData{3}), 1)]';

% Kruskal-Wallis test
[pKruskal, tblKruskal, statsKruskal] = kruskalwallis(allData, groupLabels, 'off');


if pKruskal < 0.05
    disp('Significant differences found in Kruskal-Wallis test. Proceeding with pairwise comparisons...');

    if isstring(groupLabels)
        groupLabels = cellstr(groupLabels);
    end

    % Get unique group names
    [groupNames, ~, groupIdx] = unique(groupLabels);
    numGroups = numel(groupNames);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = allData(groupIdx == i);
            data_j = allData(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [groupNames{i} ' vs ' groupNames{j}];
            raw_pvals(idx,1) = p;
            idx = idx + 1;
        end
    end

    % Bonferroni correction
    corrected_pvals = min(raw_pvals * length(raw_pvals), 1);

    % Display results
    fprintf('\nPairwise Wilcoxon Rank-Sum Test (Mann-Whitney U):\n');
    for i = 1:length(raw_pvals)
        fprintf('%s:\t raw p = %.4f,\t Bonferroni-corrected p = %.4f\n', ...
            comparisons{i}, raw_pvals(i), corrected_pvals(i));
    end
    
else
    disp('No significant differences found in Kruskal-Wallis test.');
end


%savefigure

set(fdbDec,'PaperPosition',[1,1,2.7,2.3]);
fileName=[analysisFolder filesep 'DBfullcycOneNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Sup 4B 
% plot all nights: D/B full cycle - only red
type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength)...
        & ~contains(stimTable.Remarks,'Ex') ...
        & ~any(isnan(stimTable.dbMeans),2); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curData = stimTable.dbMeans(curTrials,:);


% statistics:
[p, tbl, stats] = friedman(curData, 1,'off'); 
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
if p<0.05
    before = stimTable.dbMeans(curTrials,1);
    during = stimTable.dbMeans(curTrials,2);
    after = stimTable.dbMeans(curTrials,3);


    % Pairwise Wilcoxon signed-rank tests
    [p_before_during, ~, stats_before_during] = signrank(before, during);
    [p_during_after, ~, stats_during_after] = signrank(during, after);
    [p_after_before, ~, stats_after_before] = signrank(after, before);
    raw_pvals = [p_before_during,p_during_after,p_after_before];
    num_comparisons = length(raw_pvals);
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    %disp:
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Before vs During: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('During vs After: p-value = %.4f\n', corrected_pvals_bonferroni(2));
    fprintf('After vs Before: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
end

%plot
fdb = figure;
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curData)
    plot(curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10)
end
plot(mean(stimTable.dbMeans(curTrials,:),1),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
xlim([0.5, 3.5])
ylim([0 350])
xticks(1:3)
xticklabels(groupNames)
ylabel('D/B means during full cycle')

% savefigure
set(fdb,'PaperPosition',[1 1 2.7 2.3]);
fileName=[analysisFolder filesep 'DBchangeFullCycle'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Supplementary Figure 5

%% Sup 5A
 p = 1000*60*45;
 cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
     'headtagse_cali_recs/calibration_results.mat']).cali_result;

recName = 'Animal=PV159,recNames=Night4';
SA.setCurrentRecording(recName);
curHS = SA.recTable.Headstage{SA.currentPRec};%stimTable.Headstage(i,:);
headstageAmpCalib = cali_result.(curHS);
cur_sensativity = headstageAmpCalib.sensetivity;
cur_zeroGbias = headstageAmpCalib.zeroGbais;

SA.getLizardMovements('sensitivity', cur_sensativity(1,:)', 'zeroGBias', ...
    cur_zeroGbias(1,:)')
LM = SA.getLizardMovements;
pitch = -LM.angles(2,:);
[angleF, angleF_t] = getHeadLifts(pitch,LM.t_static_ms,100,5);

%get sleep times:
SA.getDelta2BetaRatio;
SA.getDelta2BetaAC;
AC = SA.getDelta2BetaAC;
curSleepStartT = AC.tStartSleep;
curSleepEndT = AC.tEndSleep;

fHA = figure;
plot(angleF_t/(1000*60*60), angleF,'k');
xlabel('Time (hours)'); ylabel('Head Angle')
xline(curSleepStartT/(1000*60*60),'k')
xline((curSleepStartT+p)/(1000*60*60),'g')
xline(curSleepEndT/(1000*60*60),'b')
yline(0,'--','Color',[0.4 0.4 0.4])

%save Figure
set(fHA,'PaperPositionMode','auto');
fileName=[SA.currentPlotFolder filesep 'headAngleNoStim'];
print(fileName,'-dpdf','-r300');



%% Sup 5B
i = 17;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;%('overwrite',1);
pitchAngles = -LM.angles(2,:);
[angleF, angleF_t] = getHeadLifts(pitchAngles,LM.t_static_ms,100,5);

wakeEnd = 1000*60*60*1;
if wakeEnd>stimTable.sleepStartT(i)
    wakeEnd = stimTable.sleepStartT(i);
end

parts ={[0, wakeEnd],[stimTable.sleepStartT(i),stimTable.stimStartT(i)],...
    [stimTable.stimStartT(i),stimTable.stimEndT(i)], ...
    [stimTable.stimEndT(i),stimTable.sleepEndT(i)]};
numParts = numel(parts);

% plot
f = figure;
plot(angleF_t/(1000*60*60), angleF,'k');
xlabel('Time (hours)'); ylabel('Head Angle')
xline(stimTable.sleepStartT(i)/(1000*60*60),'b')
xline(stimTable.stimStartT(i)/(1000*60*60),'r')
xline(stimTable.stimEndT(i)/(1000*60*60),'Color','r');
xline(stimTable.sleepEndT(i)/(1000*60*60),'b')
yline(0,'--','color',[0.5 0.5 0.5]); ylim([-30 70])
legend({'';'Start Sleep';'Start Stimulations';'End Stimulations';'End Sleep'})

%save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'headAnglesPV159N34'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% intensity supplamentry:

%% figure 3J - intensity shifts: delta2beta in SW change
clearvars -except stimTable SA LMData animalsColors uniqueAnimals analysisFolder stimWaveL Animals
stimWaveL = ["blue","strongBlue","white","DayTime"];
Animals = unique(stimTable.Animal);

couples = {"blue","white"};

dbSWdata = [];
groupNum = [];
dbSWMeans = [];
colorMat =[];
psfromZero =[];
ns = [];
Ns = [];

raw_pvals = zeros(1,length(couples));

for i = 1:length(couples)
    for j = 1:2
        type = (i-1)*2 +j;
        curType = stimWaveL(type);
        curTrials = contains(stimTable.Remarks,curType) & ...
            ~contains(stimTable.Remarks,'Ex') &...
            all(~isnan(stimTable.dbSWMeans(:,1:2)), 2) & ...
            ismember(stimTable.Animal,Animals);
                        
        n = sum(curTrials);
        N = length(unique(stimTable.Animal(curTrials)));
        ns = [ns; n];
        Ns = [Ns; N];

        curDBSWdata = stimTable.dbSWMeans(curTrials,2)-stimTable.dbSWMeans(curTrials,1);
        dbSWdata = [dbSWdata;curDBSWdata];
        groupNum = [groupNum; repmat(type,height(curDBSWdata),1)];
        dbSWMeans = [dbSWMeans, mean(curDBSWdata,'omitmissing')];
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :);
        colorMat = [colorMat; curColorMat];
        pfromZero = signrank(curDBSWdata,0);
        psfromZero = [psfromZero,pfromZero];
        
    end
    %statistics
    curCoupleNum = [(i*2-1),(i*2)];
    data1 = dbSWdata(groupNum==curCoupleNum(1));
    data2 = dbSWdata(groupNum==curCoupleNum(2));
    p = ranksum(data1,data2);
    raw_pvals(i) = p;


end
    num_comparisons = length(couples);
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    fprintf('DBSW: Wilcoxon rank-sum (Mann-Whitney) test results with Bonferroni correction:\n');
    fprintf('blueCouple: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('whiteCouple: p-value = %.4f', corrected_pvals_bonferroni(2));
    fprintf('\n ns:');fprintf('%d,',ns); fprintf('\n Ns:');fprintf('%d,',Ns); fprintf('\n')
    fprintf('diff from 0:');fprintf('%.4f,',psfromZero);fprintf('\n');
%plot:
f=figure;
swarmchart(groupNum,dbSWdata,20,colorMat,'filled','XJitterWidth',0.6);
hold on; scatter(1:type,dbSWMeans,'k','Marker','+')
labels = {'Blue','Strong Blue','White','Strong White'};
xticks(1:length(stimWaveL));xticklabels(labels)
ylabel('DBdiff in SW Stim-Pre')
yline(0,'--', 'Color',[0.4 0.4 0.4])

% savefigure
set(f,'PaperPosition',[1 5 4 3]);
fileName=[analysisFolder filesep 'IntensityDBSWdiffStimPre'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% figure 3K - intensity shifts: head angle SD:
clearvars -except stimTable SA LMData animalsColors uniqueAnimals analysisFolder stimWaveL Animals
headAngleSD = LMData.headAngleSD;

couples = {"blue","white"};

headAngleSDdata = [];
groupNum = [];
headAngleMeans = [];
colorMat =[];
ns = [];
Ns = [];

raw_pvals = zeros(1,length(couples));
psfromZero = [];
for i = 1:length(couples)
    for j = 1:2
        type = (i-1)*2 +j;
        curType = stimWaveL(type);
        curTrials = contains(stimTable.Remarks,curType) & ...
            ~contains(stimTable.Remarks,'Ex') &...
            ismember(stimTable.Animal,Animals);
                        
        n = sum(curTrials);
        N = length(unique(stimTable.Animal(curTrials)));
        ns = [ns; n];
        Ns = [Ns; N];

        curHeadAngleSD = headAngleSD(curTrials,3)-headAngleSD(curTrials,2);
        headAngleSDdata = [headAngleSDdata;curHeadAngleSD];
        groupNum = [groupNum; repmat(type,height(curHeadAngleSD),1)];
        headAngleMeans = [headAngleMeans, mean(curHeadAngleSD,'omitmissing')];
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :);
        colorMat = [colorMat; curColorMat];
        pFromZero = signrank(curHeadAngleSD,0);
        psfromZero = [psfromZero,pFromZero];
    end
    %statistics
    curCoupleNum = [(i*2-1),(i*2)];
    data1 = headAngleSDdata(groupNum==curCoupleNum(1));
    data2 = headAngleSDdata(groupNum==curCoupleNum(2));
    p = ranksum(data1,data2);
    raw_pvals(i) = p;


end
    num_comparisons = length(couples);
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    fprintf('headAngle: Wilcoxon rank-sum (Mann-Whitney) test results with Bonferroni correction:\n');
    fprintf('blueCouple: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('whiteCouple: p-value = %.4f', corrected_pvals_bonferroni(2));
    fprintf('\n ns:');fprintf('%d,',ns); fprintf('\n Ns:');fprintf('%d,',Ns); fprintf('\n')
    fprintf('diff from 0:');fprintf('%.4f,',psfromZero);fprintf('\n');

%plot:
f=figure;
swarmchart(groupNum,headAngleSDdata,20,colorMat,'filled','XJitterWidth',0.6);
hold on; scatter(1:type,headAngleMeans,'k','Marker','+')
labels = {'Blue','Strong Blue','White','Strong White'};
xticks(1:length(stimWaveL));xticklabels(labels)
ylabel('delta headAngle SD Stim-Pre')
yline(0,'--', 'Color',[0.4 0.4 0.4])

% savefigure
set(f,'PaperPosition',[1 5 4 3]);
fileName=[analysisFolder filesep 'IntensityHeadAngleSD'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);



%% figure 3L - intensity shifts: head mov:
clearvars -except stimTable SA LMData animalsColors uniqueAnimals analysisFolder stimWaveL Animals

couples = {"blue","white"};

headMovDdata = [];
groupNum = [];
headMovMeans = [];
colorMat =[];
ns = [];
Ns = [];

LMstimbinM = cell2mat(LMData.LMstimbin); % takes out the nan val
raw_pvals = zeros(1,length(couples));
psfromZero = [];
for i = 1:length(couples)
    for j = 1:2
        type = (i-1)*2 +j;
        curType = stimWaveL(type);
        curTrials = contains(stimTable.Remarks,curType) & ...
            ~contains(stimTable.Remarks,'Ex') &...
            ismember(stimTable.Animal,Animals);
                        
        n = sum(curTrials);
        N = length(unique(stimTable.Animal(curTrials)));
        ns = [ns; n];
        Ns = [Ns; N];
        LMpre = LMData.LMpre(curTrials);
        LMstimbintrialM = mean(LMstimbinM(curTrials,:),2); % mean for each night
        curDiffStimPre = LMstimbintrialM-LMpre;
        headMovDdata = [headMovDdata;curDiffStimPre];
        groupNum = [groupNum; repmat(type,height(curDiffStimPre),1)];
        headMovMeans = [headMovMeans, mean(curDiffStimPre,'omitmissing')];
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :);
        colorMat = [colorMat; curColorMat];
        pFromZero = signrank(curDiffStimPre,0);
        psfromZero= [psfromZero,pFromZero];
    end
    %statistics
    curCoupleNum = [(i*2-1),(i*2)];
    data1 = headMovDdata(groupNum==curCoupleNum(1));
    data2 = headMovDdata(groupNum==curCoupleNum(2));
    p = ranksum(data1,data2);
    raw_pvals(i) = p;
end
    num_comparisons = length(couples);
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    fprintf('headMov: Wilcoxon rank-sum (Mann-Whitney) test results with Bonferroni correction:\n');
    fprintf('blueCouple: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('whiteCouple: p-value = %.4f', corrected_pvals_bonferroni(2));
    fprintf('\n ns:');fprintf('%d,',ns); fprintf('\n Ns:');fprintf('%d,',Ns); fprintf('\n')
    fprintf('diff from 0:');fprintf('%.4f,',psfromZero);fprintf('\n');

%plot:
f=figure;
swarmchart(groupNum,headMovDdata,20,colorMat,'filled','XJitterWidth',0.6);
hold on; scatter(1:type,headMovMeans,'k','Marker','+')
labels = {'Blue','Strong Blue','White','Strong White'};
xticks(1:length(stimWaveL));xticklabels(labels)
ylabel('delta movement Stim-Pre')
yline(0,'--', 'Color',[0.4 0.4 0.4])

% savefigure
set(f,'PaperPosition',[1 5 4 3]);
fileName=[analysisFolder filesep 'IntensityHeadMov'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% eyelid supplamentry:
EyelidFolder = '/media/sil1/Data/Nitzan/Experiments/Eyelids_obsorption/TAMS Scan_PV143';
eyelidfilename1 =[EyelidFolder filesep 'eyelid1_pv143.Sample.Raw.csv'];
spectro1= readtable(eyelidfilename1);
eyelidfilename2 =[EyelidFolder filesep 'eyelid2_pv143.Sample.Raw.csv'];
spectro2= readtable(eyelidfilename2);
means = mean([spectro1.x_T, spectro2.x_T]);
nspectro1 = spectro1.x_T/means(1);
nspectro2 = spectro2.x_T/means(2);
avgTrans = mean([spectro1.x_T, spectro2.x_T],2);
nAvgTrans = mean([nspectro1, nspectro2],2);
%% plot
% figure; plot(spectro1.nm,spectro1.x_T,color=[0.5 0.5 0.5])
% hold on; plot(spectro2.nm,spectro2.x_T,color=[0.5 0.5 0.5]);
% plot(spectro2.nm,avgTrans,'k')

figure; plot(spectro1.nm,nspectro1,color=[0.8 0.8 0.8])
hold on; plot(spectro2.nm,nspectro2,color=[0.5 0.5 0.5]);
plot(spectro2.nm,nAvgTrans,'k')


xlabel('Wavelength (nm)'); ylabel('% Transmission');

% X regions to shade
xRanges = [ ...
    435   485;
    500   550;
    592   667;
    662  737];

% Colors (RGBA)
colors = [ ...
    0 0 1;   % blue
    0 1 0;   % green
    1 0.5 0;   % orange
    1 0 0]; % red

alphaVal = 0.15;

% Get y limits
yl = ylim;

for i = 1:size(xRanges,1)
    patch( ...
        [xRanges(i,1) xRanges(i,2) xRanges(i,2) xRanges(i,1)], ...
        [yl(1) yl(1) yl(2) yl(2)], ...
        colors(i,:), ...
        'FaceAlpha',alphaVal, ...
        'EdgeColor','none');
end

uistack(findobj(gca,'Type','line'),'top')  % keep line on top
% savefigure;
set(gcf,'PaperPosition',[1 5 4 3]);
fileName=[analysisFolder filesep 'NspectrometerEyelidSup'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% Third eye Supplamentary figure:
animalsThird = ["PV157","PV161"];

i = 25;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
plotStimSham(SA);
SA.plotDelta2BetaSlidingAC(stim=1)
%% plot the data
% plot nights without :
isThird = ismember(string(stimTable.Animal), animalsThird);   % Nx1, always
withEyeTrials   = stimTable.Mani==1 & isThird;
withoutETrials  = stimTable.Mani==2 & isThird;
nwith = sum(withEyeTrials);
Nwith= length(unique(stimTable.Animal(withEyeTrials)));
nWO = sum(withoutETrials);
NWO = length(unique(stimTable.Animal(withoutETrials)));
 
% Statistics:
[pMW_sham,h_sham] = ranksum(stimTable.dbDiffShamM(withoutETrials), stimTable.dbDiffShamM(withEyeTrials));
    fprintf('Mann-Whitney Signed-Rank Test p-value - sham: %.4f\n', pMW_sham);
    fprintf('h = %i\n',h_sham)

[pMW_stim, h_stim] = ranksum(stimTable.dbDiffStimM(withoutETrials), stimTable.dbDiffStimM(withEyeTrials));
    fprintf('Mann-Whitney Signed-Rank Test p-value - stim: %.4f\n', pMW_stim);
    fprintf('h = %i\n',h_stim)

%plot
x = 1:2;
f=figure;
h1 = plot(x,[stimTable.dbDiffShamM(withEyeTrials) stimTable.dbDiffStimM(withEyeTrials)], ...
    'Marker','.','LineStyle','-','Color',[1 0.65 0]);
hold on
h2 = plot(x,[stimTable.dbDiffShamM(withoutETrials) stimTable.dbDiffStimM(withoutETrials)], ...
    'Marker','.','LineStyle','-','Color',[11/255 102/255 35/255]);
xlim([0.9 2.1])
xticks([x]);xticklabels(["Sham","Stim"])

legend([h1(1), h2(1)],["Trials with Third Eye","Trials without Third Eye"],'Location','northwest')


% save figue
% 
set(f,'PaperPosition',[1 1 2 3]);
fileName=[analysisFolder filesep 'ThirdEye'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


