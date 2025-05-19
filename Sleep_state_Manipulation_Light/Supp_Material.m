%% Supplemental Materials
% initiating stimTable and other essentials. 

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTable.mat'])
% load([analysisFolder filesep 'LMdata.mat'])
animalsColors = [
    255/255, 142/255, 71/255;% HEX:  FF8E47 - orange  - PV126
    28/255, 144/255, 217/255;  % HEX: 1C90D9 - blue - PV149
    148/255, 120/255, 186/255; % HEX: 9478BA - perpule - PV157
    217/255, 62/255, 67/255; % HEX: D93E43 - red - PV159
    255/255, 202/255, 58/255; % HEX: FFCA3A - yellow -  PV161
    97/255, 184/255, 68/255;  % HEX:61B844 - Green -PV162
];
uniqueAnimals = unique(stimTable.Animal);


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





%% Supplementary Figure 2
%% Sup 2B
%% plot all raw traces, 1 nights:

i = 17;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
    
% get the stimualtions time stamps:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch}; %stimulations timimng in ms
    
% get the raw data for the second around each stim:
pre = 200; 
post = 1000;
[mV, mV_t] = SA.currentDataObj.getData(SA.recTable.defaulLFPCh(SA.currentPRec),stims-pre,(pre+post));

% plot
fraw = figure;
x = linspace(-pre,post,length(mV_t));
mV1 = squeeze(mV);
[d,h] = hist2(repmat(x,height(mV1),1),mV1,'dX1',3,'dX2',5);
hold on
plot(downsample(x,100),downsample(mean(mV1,1),100),'Color','r','LineWidth',2)
ylims =[-750 500];
ylim(ylims)
% title('Raw traces for 1 night')
x_shade = [0 200 200 0];  % X-coordinates of the shaded region
y_shade = [ylims(1) ylims(1) ylims(2) ylims(2)]; % Y-coordinates covering the full y-range
patch(x_shade, y_shade, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
xlabel('Time(ms)');ylabel('Voltage[uV]')

% save figure
set(fraw ,'PaperPosition',[1 2 4.3 1.7]);
fileName=[analysisFolder filesep 'SWRrawPV159N34'];
% print(fileName,'-depsc','-vector');
print(fileName, '-dpdf', '-r300');

%% Sup 2C+D, plot cros corr 1 night:
% detect SWR using SWR detection method.

% for 1 night first:
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
SA.getSharpWaves('detectOnlyDuringSWS',false)
SW = SA.getSharpWaves('detectOnlyDuringSWS',false);
swrT = SW.tSW'; %SWR timings in ms
swrN = length(swrT); 

% get the stimualtions time stamps:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch}; %stimulations timimng in ms
first_stims = stims(1:8:end-2);
% cross cor: Short, one night 
% Convert Timestamps to Binary Time Series:
bin_size = 200; % ms 
winShort = 10*1000; % ten seconds

ylims = [-0.04 0.1];
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
 
%% Supp 2E+F
fName = [analysisFolder filesep 'ShWcrossCorrData.mat'];
load(fName)

wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) ...
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


%% Supplementary Figure 4
%% Sup 4A: plot one night: D/B decrease full cycle

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

%% Supp 4B

%% plot all nights: D/B full cycle - only red
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
[p, tbl, stats] = friedman(curData, 1); % Here, 1 indicates within-subjects design
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
% plot(stimTable.dbSWMeans(curTrials,:)','Color',[0.5 0.5 0.5],'Marker','.','MarkerSize',10)
%create color code:

[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
hold on; 
for i = 1:height(curData)
    plot(curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10)
end
plot(mean(stimTable.dbMeans(curTrials,:),1,'omitnan'),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
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

%% Supplementary Figure 6


