%% Figure 2:
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

%% Figure 2A - one night: D/B decrease

i =74; % PV106,
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);

dbStimData = {stimTable.dbSW(i).Pre stimTable.dbSW(i).Stim stimTable.dbSW(i).Post};

fdbDec = figure;
groupNames = ["Pre", "Stim","Post"];
colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];

yswarm = [stimTable.dbSW(i).Pre stimTable.dbSW(i).Stim stimTable.dbSW(i).Post]; %same data in one array
xswarm=[1*ones(1,length(stimTable.dbSW(i).Pre)), 2*ones(1,length(stimTable.dbSW(i).Stim)), 3*ones(1,length(stimTable.dbSW(i).Post))];
colorswarm = [repmat([0.2, 0.6, 0.8;], length(stimTable.dbSW(i).Pre), 1);  % Red for Array 1
    repmat([0.9, 0.4, 0.3], length(stimTable.dbSW(i).Stim), 1);  % Green for Array 2
    repmat([0.5, 0.8, 0.5], length(stimTable.dbSW(i).Post), 1)]; % Blue for Array 3
swarmchart(xswarm,yswarm,20,colorswarm,"filled",'o','XJitterWidth',0.5)
ylabel('D/B power in SWS bouts')
% ylim([0 450])
set(gca, 'XTick', 1:numel(dbStimData), 'XTickLabel', groupNames);
title('D/B power during SWS bouts - one night')

% statistics:

% Combine data into a single vector and create a grouping variable
allData = [dbStimData{1},dbStimData{2},dbStimData{3}]'; % Concatenate all data
groupLabels = [...
    repmat({'Pre'}, numel(dbStimData{1}), 1); ...
    repmat({'During'}, numel(dbStimData{2}), 1); ...
    repmat({'After'}, numel(dbStimData{3}), 1)]';

% Kruskal-Wallis test
[pKruskal, tblKruskal, statsKruskal] = kruskalwallis(allData, groupLabels, 'off');
fprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal);
annotation('textbox', [0.1, 0.8, 0.2, 0.1], 'String', ...
    sprintf('Kruskal-Wallis test p-value: %.4f\n', pKruskal), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


if pKruskal < 0.05
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
    corrected_pvals = min(raw_pvals * length(raw_pvals),1);

    % Display results
    fprintf('\nPairwise Wilcoxon Rank-Sum Test (Mann-Whitney U):\n');
    for i = 1:length(raw_pvals)
        fprintf('%s:\t raw p = %.4f,\t Bonferroni-corrected p = %.4f\n', ...
            comparisons{i}, raw_pvals(i), corrected_pvals(i));
    end


end

%savefigure
set(fdbDec,'PaperPosition',[1 5 2 2]);
fileName=[analysisFolder filesep 'DBSWSpreStimPostOneNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% Figure 2B - D/B decrease all white nights
% plot all red nights 

type = 'white';
wavelength = 'white'; 
curTrials = (contains(stimTable.Remarks,wavelength)|contains(stimTable.Remarks,'DayTime') ) &...
    ~contains(stimTable.Remarks,'Ex') ...
    & ~any(isnan(stimTable.dbSWMeans),2); 
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Stim-Pre', 'Stim-Post'};
diffStimPre = stimTable.dbSWMeans(curTrials,2)-stimTable.dbSWMeans(curTrials,1);
diffStimPost = stimTable.dbSWMeans(curTrials,2)-stimTable.dbSWMeans(curTrials,3);
plotdata = [diffStimPre;diffStimPost];
x=[repmat(1,[n,1]);repmat(2,[n,1])];

%plot
fdb = figure;
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
plotdata = [diffStimPre;diffStimPost];
x=[repmat(1,[n,1]);repmat(2,[n,1])];
curColorMat = repmat(animalsColors(animalIndices, :),[length(groupNames),1]); 
swarmchart(x,plotdata,30,curColorMat,'filled','XJitter','randn','XJitterWidth',0.8)
hold on; scatter([1:length(groupNames)],[mean(plotdata(x==1)),mean(plotdata(x==2))],50,'k','Marker','+');
yline(0,'Color',[0.4 0.4 0.4],'LineStyle','--')
xlim([0.5, 2.5])
% ylim([-300 100])
xticks(1:2);xticklabels(groupNames); ylabel('D/B means during SWS')

%statistics 
% check if different from zero, 
[pfromZeroPre, h] = signrank(diffStimPre, 0);
[pfromZeroPost, h] = signrank(diffStimPost, 0);

fprintf('stim-pre diff:\t p = %.4f\n', pfromZeroPre)
fprintf('stim-post diff:\t p = %.4f\n', pfromZeroPost)

% savefigure
set(fdb,'PaperPosition',[1 5 2 2]);
fileName=[analysisFolder filesep 'DBSWSwhiteNightsDiffs'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

% groupNames = {'Pre', 'During', 'After'};
% curData = stimTable.dbSWMeans(curTrials,:);
% 
% 
% % statistics:
% [p, tbl, stats] = friedman(curData, 1,'off'); % paired data
% fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% if p<0.05
%     % data for three groups
%     before = stimTable.dbSWMeans(curTrials,1);
%     during = stimTable.dbSWMeans(curTrials,2);
%     after = stimTable.dbSWMeans(curTrials,3);
% 
%     % Pairwise Wilcoxon signed-rank tests
%     [p_before_during, ~, stats_before_during] = signrank(before, during);
%     [p_during_after, ~, stats_during_after] = signrank(during, after);
%     [p_after_before, ~, stats_after_before] = signrank(after, before);
% 
%     raw_pvals = [p_before_during,p_during_after,p_after_before];
%     num_comparisons = 3;
%     corrected_pvals_bonferroni = min(raw_pvals * num_comparisons,1);
%     fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
%     fprintf('Before vs During: p-value = %.4f \n', corrected_pvals_bonferroni(1));
%     fprintf('During vs After: p-value = %.4f\n', corrected_pvals_bonferroni(2));
%     fprintf('After vs Before: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
% end 
% %plot
% fdb = figure;
% [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
% curColorMat = animalsColors(animalIndices, :); 
% hold on; 
% for i = 1:height(curData)
%     plot(curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10)
% end
% plot(mean(stimTable.dbSWMeans(curTrials,:),1,'omitnan'),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
% xlim([0.5, 3.5])
% % ylim([0 450])
% xticks(1:3)
% xticklabels(groupNames)
% ylabel('D/B means during SWS')
% 
% annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
%     sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');
% 
% annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
%     sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');
% 
% % savefigure
% set(fdb,'PaperPosition',[1 5 2 2]);
% fileName=[analysisFolder filesep 'DBSWSwhiteNights'];
% print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Figure 2C+D - exmple for head movement (clean)
% plot the full movement data for a night
close all
clearvars -except stimTable SA LMData animalsColors uniqueAnimals analysisFolder
% load([analysisFolder filesep 'LMDataAll.mat'])
% for one night:

% recInd = (contains(stimTable.Remarks,"DayTime")|contains(stimTable.Remarks,"white"));
% animals =stimTable.Animal(recInd);recnames = stimTable.recNames(recInd);
% stimTableInd = find(recInd);
% recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
%  for i = 1:length(recList)
    % i=13;

% SA.setCurrentRecording(recList{i})
recName = 'Animal=PV153,recNames=Night15';
SA.setCurrentRecording(recName)
DB = SA.getDelta2BetaRatio;
% ind = stimTableInd(i);
ind = 69;
curLMwake = LMData.LMwake(ind);
curLMpre = LMData.LMpre(ind);
curLMstim = LMData.LMstimbin{ind};
curLMpost = LMData.LMpost(ind);

figure;

% Total number of horizontal slots = 5
% Subplot widths:
w_small = 0.12;  % Width for small subplots (1/5)
w_large = 0.39;  % Width for the large subplot (2/5)

h = 0.8;  % Height of all subplots (80% of the figure height)
bottom = 0.1;  % Distance from the bottom edge of the figure

% Subplot 1: Position manually
left1 = 0.05;  % Left edge for subplot 1
a1 = subplot('Position', [left1, bottom, w_small, h]);
plot(1, curLMwake, '.', 'Color', 'black', 'MarkerSize', 20);
xticks(1); xticklabels('Wake');
ylabel('Mov/s');

% Subplot 2: Position manually
left2 = left1 + w_small + 0.05;  % Space after subplot 1
a2 = subplot('Position', [left2, bottom, w_small, h]);
plot(1, curLMpre, '.', 'Color', 'black', 'MarkerSize', 20);
xticks(1); xticklabels('Sleep before');

% Subplot 3: Larger width
left3 = left2 + w_small + 0.05;  % Space after subplot 2
a3 = subplot('Position', [left3, bottom, w_large, h]);
xdur = (1:length(curLMstim))*10;
plot(xdur, curLMstim, '-o', 'Color', 'black', 'Marker','.','MarkerSize',20);
xlabel('Stimulations Avg.: Time from start trial');
% ylabel('Trials Avg.');
hold on;
xline(0, 'Color', 'r', 'LineWidth', 2);
xline(37, 'Color', 'r', 'LineWidth', 2);
hold off;

% Subplot 4: Position manually
left4 = left3 + w_large + 0.05;  % Space after subplot 3
a4 = subplot('Position', [left4, bottom, w_small, h]);
plot(1, curLMpost, '.', 'Color', 'black', 'MarkerSize', 20);
xticks(1); xticklabels('Sleep After');

% Link y-axes and set limits
linkaxes([a1, a2, a3, a4], 'y'); %ylim([0 6.5]);

% Add 4 shared title
% sgtitle(['Mean movement during stimulation,' recList{i}]);
sgtitle(['Mean movement during stimulation,' recName]);
% end
% % save fig
set(gcf,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'lizMovWholeNightPV106N20'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);



%% Figure 2E - movement rate difference between stim times vs no stim
wavelength = 'white';
curTrials = (contains(stimTable.Remarks,wavelength)|contains(stimTable.Remarks,'DayTime')) & ...
    ~contains(stimTable.Remarks,'Ex');

n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

LMstim = cell2mat(LMData.LMstimbin(curTrials));
stimMov = mean(LMstim(:,1:4),2);
notStimMov = mean(LMstim(:,5:end),2);
plotdata = [stimMov,notStimMov];

f=figure;
Groups = ["Stim","After Stim"];
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :);
x = 1:2;
hold on 
for i = 1:length(plotdata)
    plot(x,plotdata(i,:),'Color',curColorMat(i,:),'Marker','.')
end
plot(x, mean(plotdata),'Color','k','Marker','.')
xticks(x), xticklabels(Groups); xlim([0.7 2.3]);
ylabel('Mov/S')
% ylim([0 37]);

%statistics:
[pWilcoxon, ~, statsWilcoxon] = signrank(stimMov,notStimMov);
fprintf('Wilcoxon Signed-Rank Test p-value: %.4f\n', pWilcoxon);
disp(statsWilcoxon)


% savefigure
set(f,'PaperPosition',[1 1 3 2]);
fileName=[analysisFolder filesep 'stim-notstimcyclesWhite'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Figure 2F
% plot head movements - white nights!
close all
clearvars -except stimTable SA LMData animalsColors uniqueAnimals analysisFolder

% load([analysisFolder filesep 'LMDataAll.mat'])
wavelength = 'white';
curTrials = (contains(stimTable.Remarks,wavelength)|contains(stimTable.Remarks,'DayTime')) & ...
    ~contains(stimTable.Remarks,'Ex');

n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

LMpre = LMData.LMpre(curTrials, :);
LMwake = LMData.LMwake(curTrials, :);
LMstimbinM = cell2mat(LMData.LMstimbin); % takes out the nan val
LMstimbintrialM = mean(LMstimbinM(curTrials,:),2); % mean for each night

LMplotData = [LMwake, LMpre, LMstimbintrialM];

% check the statistics:
[p, tbl, stats] = friedman(LMplotData, 1,'off'); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
if p<0.05
    % Pairwise Wilcoxon signed-rank tests
    [p_wake_pre, ~, stats_wake_pre] = signrank(LMwake, LMpre);
    [p_wake_during, ~, stats_wake_during] = signrank(LMwake, LMstimbintrialM);
    [p_pre_during, ~, stats_pre_during] = signrank(LMstimbintrialM,LMpre);

    raw_pvals = [p_wake_pre,p_pre_during,p_wake_during];
    num_comparisons = length(raw_pvals); 
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    % Display results with Bonferroni correction
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Wake vs pre: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('pre vs During: p-value = %.4f\n', corrected_pvals_bonferroni(2));
    fprintf('Wake vs During: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
end


fLMr = figure;
Groups = ["Wake","Pre","During Stim"];

[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :);
x = 1:3;
hold on 
for i = 1:length(LMplotData)
    plot(x,LMplotData(i,:),'Color',curColorMat(i,:),'Marker','.')
end
plot(x, mean(LMplotData),'Color','k','Marker','.')
xticks(x), xticklabels(Groups); xlim([0.7 3.3]);
ylabel('Mov/S')
ylim([0 37]);

% savefigure
set(fLMr,'PaperPosition',[1 1 3 2]);
fileName=[analysisFolder filesep 'LMwhiteNights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Figure 2G - Head angle SD
% load([analysisFolder filesep 'LMData.mat'])
headAngleSD = LMData.headAngleSD;

wavelength = 'white';
curTrials = (contains(stimTable.Remarks,wavelength)|contains(stimTable.Remarks,'DayTime')) & ...
    ~contains(stimTable.Remarks,'Ex');
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

groupNames = {'Wake-Stim', 'Post-Stim'};
diffWakeStim = headAngleSD(curTrials,1)-headAngleSD(curTrials,3);
diffStimPre = headAngleSD(curTrials,2)-headAngleSD(curTrials,3);

%plot
fdb = figure;
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
plotdata = [diffWakeStim;diffStimPre];
x=[repmat(1,[n,1]);repmat(2,[n,1])];
curColorMat = repmat(animalsColors(animalIndices, :),[length(groupNames),1]); 
swarmchart(x,plotdata,30,curColorMat,'filled','XJitter','randn','XJitterWidth',0.8)
hold on; scatter([1:length(groupNames)],[mean(plotdata(x==1)),mean(plotdata(x==2))],50,'k','Marker','+');
yline(0,'Color',[0.4 0.4 0.4],'LineStyle','--')
xlim([0.5, 2.5])
% ylim([0 450])
xticks(1:2);xticklabels(groupNames); ylabel('D/B means during SWS')

%statistics 
% check if different from zero, 
[pfromZeroWake, h] = signrank(diffWakeStim, 0);
[pfromZeroPre, h] = signrank(diffStimPre, 0);

fprintf('Wake-stim diff:\t p = %.4f\n', pfromZeroWake)
fprintf('pre-stim diff:\t p = %.4f\n', pfromZeroPre)

% savefigure
set(fdb,'PaperPosition',[1 5 2 2]);
fileName=[analysisFolder filesep 'HeadAngleSDdiffsColorsExt'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);




% groupNames = {'Pre', 'During', 'After'};
% curHeadSD = headAngleSD(curTrials,1:3);
% 
% figure;
% x1 = 1:width(curHeadSD);
% [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
% curColorMat = animalsColors(animalIndices, :); 
% for i= 1:height(curHeadSD)
%     plot(x1,curHeadSD(i,:),'Color',curColorMat(i,:),'Marker','.'); hold on;
% end
% plot(x1, mean(curHeadSD), 'Color','k','Marker','.','LineWidth',1.5)
% xlim([0.7,3.2]); xticks(x1); xticklabels(["Wake","Pre","Stim"])
% ylabel('Avg Head SD')
% title('Head Angle SD - white nights')
% 
% % stats:
% [p, tbl, stats] = friedman(curHeadSD, 1,'off'); % Here, 1 indicates within-subjects design
% fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% % p-valure is very low, post hoc:
% if p<0.05
%     wakeSD = curHeadSD(:,1);
%     beforeSD = curHeadSD(:,2);
%     duringSD = curHeadSD(:,3);
% 
%     % Pairwise Wilcoxon signed-rank tests
%     [p_wake_before, ~, stats_wake_before] = signrank(wakeSD, beforeSD);
%     [p_before_during, ~, stats_before_during] = signrank(beforeSD, duringSD);
%     [p_wake_during, ~, stats_wake_during] = signrank(wakeSD, duringSD);
% 
%     raw_pvals = [p_wake_before,p_before_during,p_wake_during];
%     num_comparisons = length(raw_pvals);
%     corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
%     % Display results with Bonferroni correction
%     fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
%     fprintf('Wake vs pre: p-value = %.4f \n', corrected_pvals_bonferroni(1));
%     fprintf('pre vs During: p-value = %.4f\n', corrected_pvals_bonferroni(2));
%     fprintf('Wake vs During: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
% end
% % savefigure
% set(gcf,'PaperPosition',[1 4 2.2 1.6])
% saveas (gcf, [analysisFolder filesep 'HeadLiftswhiteNightsSD.pdf']);

%% Eye Movement: Figure 2H - single night
%% one night
i = 63;%
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
SA.plotSyncedDBEyeMovements;

%% Figure2I - eyemovements - all white nights:
uniqueAnimals = unique(stimTable.Animal);
relativePhaseEye = stimTable.eyeMovPh-stimTable.eyeMovPhDB;


wavelength = 'white';
trials = stimTable.eyeMov ==1 &...
    ~contains(stimTable.Remarks,'Ex') &...
    (contains(stimTable.Remarks,wavelength)|contains(stimTable.Remarks,'DayTime'));
n=sum(trials);
N = length(unique(stimTable.Animal(trials)));
figure;
h1=polaraxes;hold on;
title('eye Mov. white nights')
Rlim=0.5;
% stats:
relativeMean = circ_mean(relativePhaseEye(trials));
dir = circ_ang2rad(180);
    [pval, v] = circ_vtest(relativePhaseEye(trials), dir);
fprintf('v test p-value: %f\n', pval);

% [p, z] = circ_rtest(relativePhaseEye(pVal));
% fprintf('Rayleigh''s test p-value: %f\n', p);
% if p < 0.05
%     fprintf('The data is significantly non-uniform.\n');
% else
%     fprintf('The data uniformly distributed.\n');
% end


hP={};
for i=1:numel(uniqueAnimals)
    p=find(trials & strcmp(stimTable.Animal,uniqueAnimals(i)));
    hP{i}=polarplot([relativePhaseEye(p)';relativePhaseEye(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
end
hold on;
hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);
hP3=polarplot([relativeMean relativeMean],[0,Rlim],'color','k','LineWidth',3);

% hRose=polarhistogram(h1,relativePhaseEye(trials),12,'Normalization','probability');
% hRose.FaceColor=[0.7 0.7 0.7];
% hRose.FaceAlpha=0.5;

text(0.2, Rlim/2, '\delta/\beta');
h1.ThetaTick=[0:90:330];
h1.RTick=[0.1:0.1:0.4];
l1=legend([hP{1}(1),hP3,hRose],{'singleNight','\delta/\beta'},'box','off');
l1.Position=[0.7386    0.8238    0.2125    0.1190];

% save figure:
fileName=[analysisFolder filesep 'EyeMovWhite2'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);% clearvars -except stimTable SA analysisFolder LMdata


%% Figure 3I
% %load data:
% load([analysisFolder filesep 'polarhistoAllNights.mat'])
% 
% % plot polar histogram - red nights
% uniqueAnimals = unique(stimTable.Animal);
% relativePhasePre = (mPhasePre.movs -mPhasePre.DBs);
% relativePhaseStim = (mPhaseStim.movs -mPhaseStim.DBs);
% relativePhasePost = (mPhasePost.movs -mPhasePost.DBs);
% 
% wavelength = 'white';
% pVal = mPhaseStim.movs ~= 0 & ...
%     (contains(stimTable.Remarks,wavelength)|contains(stimTable.Remarks,'DayTime'))...
%     & ~contains(stimTable.Remarks,'Ex');
% n=sum(pVal);
% N=length(unique(stimTable.Animal(pVal)));
% fMOVdbred=figure;
% h1=subplot(1,3,1,polaraxes);hold on;
% title('Pre Stimulations')
% Rlim=0.5;
% relativePreMean = circ_mean(relativePhasePre(pVal));
% 
% hP={};
% for i=1:numel(uniqueAnimals)
%     p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
%     hP{i}=polarplot([relativePhasePre(p)';relativePhasePre(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
% end
% hold on;
% hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);
% hP3=polarplot([relativePreMean relativePreMean],[0,Rlim],'color','k','LineWidth',3);
% 
% hRose=polarhistogram(h1,relativePhasePre(pVal),12,'Normalization','probability');
% hRose.FaceColor=[0.7 0.7 0.7];
% hRose.FaceAlpha=0.5;
% 
% text(0.2, Rlim/2, '\delta/\beta');
% h1.ThetaTick=[0:90:330];
% h1.RTick=[0.1:0.1:0.4];
% l1=legend([hP{2}(1),hP3,hRose],{'singleNight','\delta/\beta','Prob.'},'box','off');
% l1.Position=[0.7386    0.8238    0.2125    0.1190];
% 
% %figure 2 : during stim:
% h2=subplot(1,3,2,polaraxes);hold on;% Stimulation time
% title('During Stimulations')
% Rlim=0.5;
% relativeStimMean = circ_mean(relativePhaseStim(pVal));
% hP={};
% for i=1:numel(uniqueAnimals)
%     p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
%     hP{i}=polarplot([relativePhaseStim(p)';relativePhaseStim(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
% end
% hold on;
% hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);
% hP3=polarplot([relativeStimMean relativeStimMean],[0,Rlim],'color','k','LineWidth',3);
% 
% hRose=polarhistogram(h2,relativePhaseStim(pVal),12,'Normalization','probability');
% hRose.FaceColor=[0.7 0.7 0.7];
% hRose.FaceAlpha=0.5;
% 
% text(0.2, Rlim/2, '\delta/\beta');
% h2.ThetaTick=[0:90:330];
% h2.RTick=[0.1:0.1:0.4];
% 
% %figure 3 : post stim:
% h3=subplot(1,3,3,polaraxes);hold on;% Stimulation time
% title('Post Stimulations')
% Rlim=0.5;
% relativePostMean = circ_mean(relativePhasePost(pVal));
% 
% hP={};
% for i=1:numel(uniqueAnimals)
%     p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
%     hP{i}=polarplot([relativePhasePost(p)';relativePhasePost(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
% end
% hold on;
% hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);
% hP3=polarplot([relativePostMean relativePostMean],[0,Rlim],'color','k','LineWidth',3);
% 
% hRose=polarhistogram(h3,relativePhasePost(pVal),12,'Normalization','probability');
% hRose.FaceColor=[0.7 0.7 0.7];
% hRose.FaceAlpha=0.5;
% 
% text(0.2, Rlim/2, '\delta/\beta');
% h3.ThetaTick=[0:90:330];
% h3.RTick=[0.1:0.1:0.4];
% 
% 
% % savefigure
% set(gcf, 'PaperUnits', 'inches');         % Set paper units to inches
% set(gcf, 'PaperSize', [6,4]);            % Set the paper size (width x height in inches)
% set(gcf, 'PaperPosition', [0, 0, 6,4]);  % Set position on paper to match size exactly
% 
% % Print to PDF
% % set(f,'PaperPositionMode','auto');
% fileName=[analysisFolder filesep 'PolarMovDBwhiteNightold'];
% print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);% clearvars -except stimTable SA analysisFolder LMdata
% 