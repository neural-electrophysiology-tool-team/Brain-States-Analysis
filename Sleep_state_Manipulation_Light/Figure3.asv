%% Figure 3
% initiating stimTable and other essentials. 

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTable.mat'])
load([analysisFolder filesep 'LMdata.mat'])
animalsColors = [
    255/255, 142/255, 71/255;% HEX:  FF8E47 - orange  - PV126
    28/255, 144/255, 217/255;  % HEX: 1C90D9 - blue - PV149
    148/255, 120/255, 186/255; % HEX: 9478BA - perpule - PV157
    217/255, 62/255, 67/255; % HEX: D93E43 - red - PV159
    255/255, 202/255, 58/255; % HEX: FFCA3A - yellow -  PV161
    97/255, 184/255, 68/255;  % HEX:61B844 - Green -PV162
];
uniqueAnimals = unique(stimTable.Animal);

%% Figure 3A - Sham-stim all nights:
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
statsStimSham = struct();
diffData = [];
groupNames = [];
colorMat =[];
ns = [];
Ns = [];
diffmeans = [];
diffSDs = [];


f=figure;
x=1:2;
for type = 1:numType
    h = subplot(1,numType,type);
    curType = stimWaveL(type);
    curName = stimType(type);
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,'Ex') & ...
                all(~isnan(stimTable.dbDiffStimM),2) & ...
                all(~isnan(stimTable.dbDiffShamM),2); 
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    ns = [ns; n];
    Ns = [Ns; N];
    curData = [stimTable.dbDiffShamM(curTrials), stimTable.dbDiffStimM(curTrials)];
    curDatadiff = [stimTable.dbDiffStimM(curTrials)-stimTable.dbDiffShamM(curTrials)];
    diffData = [diffData; curDatadiff];
    diffmeans = [diffmeans; mean(curDatadiff)];
    diffSDs = [diffSDs;std(curDatadiff)];
    groupNames = [groupNames; repmat(type, length(curDatadiff), 1)];
    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :);
    colorMat = [colorMat; curColorMat];
    hold on
    for i = 1:height(curData)
        plot(x, curData(i,:), 'Color',curColorMat(i,:), 'Marker','.', 'MarkerSize',10, 'LineWidth',1)
    end
    plot(x,mean(curData),'Color','k', 'Marker','.', 'MarkerSize',10, 'LineWidth',1.5)
    ylim([-50 200])

    xticks(1:2); xticklabels(["Sham", "Stim"])
    xlim([0.5, 2.5]);
    %wilcoxon:
    p = signrank(curData(:,1), curData(:,2));
    fprintf('%s: p-value = %.5f\n',stimType(type), p);
end

% savefigure
set(f,'PaperPosition',[1 1 4 2]);
fileName=[analysisFolder filesep 'meanNormBDStimSham'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

% add titles. labels...
ylabel('Diff in D2B power')

% statistics: between groups diffs:
% 1. kruskal wallas:
[pkruskal, tbl, stats] = kruskalwallis(diffData,groupNames,'off');

if pkruskal < 0.05
    if isstring(groupNames)
        groupNames = cellstr(groupNames);
    end

    % Get unique group names
    [groupName, ~, groupIdx] = unique(groupNames);
    numGroups = numel(groupName);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = diffData(groupIdx == i);
            data_j = diffData(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [num2str(groupName(i)) ' vs ' num2str(groupName(j))];
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
end 



%% Figure 3B - AC all colors
% plot AC period time changes - According to color and animal
close all
clearvars -except stimTable SA LMdata animalsColors uniqueAnimals analysisFolder

stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
stimData = [];
groupNum=[];
colorMat = [];
Ns = [];
ns = [];
x=1:3;


f=figure;
set(f, 'Position', [100, 100, 800, 400]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curName = stimType(type);
    curTrials = contains(stimTable.Remarks,curType) & ...
                ~contains(stimTable.Remarks,'Ex') &...
                all(~isnan(stimTable.ACcomPer),2) &...
                all(stimTable.ACcomP2V > 0.15,2);
    n = sum(curTrials);
    N = numel(unique(stimTable.Animal(curTrials)));
    ns = [ns; n];
    Ns = [Ns; N];
    curData = stimTable.ACcomPer(curTrials,:);
    curMean = mean(curData,1,'omitnan');
    stimData = [stimData; curData(:,2)];
    groupNum = [groupNum; repmat(type, length(curData), 1)];

%statistics:

    [p, tbl, stats] = friedman(curData, 1,'off'); % Here, 1 indicates within-subjects design
    fprintf('p-value for freidman ANOVA test: %.5f\n',p)
    % % p-valure is very low, post hoc:
    if p<0.05
        % data for the four groups
        beforeAC = curData(:,1);
        duringAC = curData(:,2);
        afterAC = curData(:,3);

        % Pairwise Wilcoxon signed-rank tests
        [p_pre_during, ~, stats_before_during] = signrank(beforeAC, duringAC);
        [p_during_post, ~, stats_during_after] = signrank(duringAC, afterAC);
        [p_pre_post, ~, stats_wake_during] = signrank(beforeAC, afterAC);

        raw_pvals = [p_pre_during,p_during_post,p_pre_post];
        num_comparisons = 3;
        corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);

        % Display results with Bonferroni correction
        fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
        fprintf('Before vs During: p-value = %.4f \n', corrected_pvals_bonferroni(1));
        fprintf('During vs After: p-value = %.4f\n', corrected_pvals_bonferroni(2));
        fprintf('After vs Before: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
    
        %save statistics:
        statsAC.(curName).pAnova = p;
        statsAC.(curName).corrected_pvals_bonferroni = corrected_pvals_bonferroni;
    end
    
    if n>0
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :);
        colorMat = [colorMat; curColorMat];
        hold on;
        for j = 1:height(curData)
            plot(x,curData(j,:)/1000,'Color',curColorMat(j,:),'Marker','.')
        end
        plot(x,curMean/1000,'Color','k','LineWidth',2,'Marker','.')
        yline(156,'--', 'Color',[0.4 0.4 0.4])

    end
    if type==1
        ylabel('Time[s]')
    end
    ylim([40 250])
    xticks(1:3); xlim([0.8 3.2])
    xticklabels({'Pre','During','Post'})
  
end

% savefigure
set(gcf,'PaperPosition',[1 1 4.3 1.2]);
fileName=[analysisFolder filesep 'ACperiodstimAll'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
save([analysisFolder filesep 'statsACperiodstimAll.mat'], "statsAC")


%% Figure 3C - Pv2 at 156 s
close all
clearvars -except stimTable SA LMdata animalsColors uniqueAnimals analysisFolder


stimP2V = [];
groupNum = [];
colorMat = [];
means = [];
SDs = [];
Ns = [];
ns = [];

stimType = ["Blue","Green","Red","WhiteEx"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);

hold on
for type = 1:numType
    curType = stimWaveL(type);
    curName = stimType(type);
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,"Ex");
    n = sum(curTrials);
    N = numel(unique(stimTable.Animal(curTrials)));
    ns = [ns;n];
    Ns = [Ns; N];
    curData = stimTable.P2V_156(curTrials,:);
    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :);
    stimP2V = [stimP2V; curData(:,2)];
    means = [means, mean(curData(:,2))];
    SDs = [SDs,std(curData(:,2))];
    groupNum = [groupNum; repmat(type,height(curData),1)];
    colorMat = [colorMat; curColorMat];
    
end
f=figure;
x=1:3;
swarmchart(groupNum,stimP2V,10,colorMat,'filled','XJitterWidth',0.5);
hold on; scatter(1:4,means,'k','Marker','+')
xticks(1:4);xticklabels(stimType)
ylabel('P2V in 156s in Stim')

% 
% % savefigure
set(gcf,'PaperPosition',[1 1 3.5 3]);
fileName=[analysisFolder filesep 'P2V156n'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%

[p, tbl, statsP2Vfs] = kruskalwallis(stimP2V,groupNum,'off');

if p < 0.05
    % figure;
    % cP2Vfs = multcompare(statsP2Vfs, 'CType', 'dunn-sidak');
    % Convert group to cell array of strings if needed
    if isstring(groupNum)
        groupNum = cellstr(groupNum);
    end

    % Get unique group names
    [groupNames, ~, groupIdx] = unique(groupNum);
    numGroups = numel(groupNames);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = stimP2V(groupIdx == i);
            data_j = stimP2V(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [num2str(groupNames(i)) ' vs ' num2str(groupNames(j))];
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
end


%% Figure 3D - D/B decrease
close all
clearvars -except stimTable SA LMdata animalsColors uniqueAnimals analysisFolder

stimType = ["Blue","Green","Red","WhiteEx"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
DBdiff = [];
groupNum = [];
colorMat = [];
ns = [];
Ns = [];
psFromZero = [];
DBMean = [];

for type = 1:numType
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,"Ex") &...
                all(~isnan(stimTable.dbSWMeans(:,1:2)), 2);
    ns =[ns; sum(curTrials)];
    Ns = [Ns; numel(unique(stimTable.Animal(curTrials)))];
    curData = stimTable.dbSWMeans(curTrials,2)-stimTable.dbSWMeans(curTrials,1);
    DBdiff = [DBdiff; curData];
    DBMean = [DBMean; mean(curData)];
    groupNum = [groupNum; repmat(type,length(curData),1)];
    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :);
    colorMat = [colorMat;curColorMat];
    [pfromZero, h] = signrank(curData, 0);
    psFromZero = [psFromZero , pfromZero];
end

%stats:
[p, tbl, statsDBdiff] = kruskalwallis(DBdiff,groupNum,'off');

if p < 0.05
    if isstring(groupNum)
        groupNum = cellstr(groupNum);
    end

    % Get unique group names
    [groupNames, ~, groupIdx] = unique(groupNum);
    numGroups = numel(groupNames);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = DBdiff(groupIdx == i);
            data_j = DBdiff(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [num2str(groupNames(i)) ' vs ' num2str(groupNames(j))];
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
end

%Plot
f=figure;
swarmchart(groupNum,DBdiff,10,colorMat,'filled','XJitterWidth',0.5);
hold on;
scatter(1:4, DBMean,14,'k','Marker','+')
ylabel('d/b diff')
xticks(1:numType); xticklabels(stimType)
yline(0,'Color',[0.4 0.4 0.4],'LineStyle','--')
ylim([-350 150])

% savefigure
set(f,'PaperPosition',[2 1 2.5 1.5]);
fileName=[analysisFolder filesep 'DBdecreaseDiffAllnigthscolors'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% Figure 3E

load([analysisFolder filesep 'LMData.mat'])
headAngleSD = LMData.headAngleSD;
% HeadAngleAvg = LMData.HeadAngleAvg;

% plot SDs
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
statsHeadangSd = struct();

% mean normelized change in D/B.
fHA=figure;
sgtitle('Head Angles over night')
headSDdiff = [];
groupNum = [];
headSDdiffM = [];
colorMat = [];
psFromZero = [];
ns = [];
Ns = [];
for type = 1:numType

    %create data subset:
    curName = stimType(type);
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType)& ~contains(stimTable.Remarks,'Ex' );
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :); 
    curHeadSD = headAngleSD(curTrials,1:3);
    curHeadSDmean = mean(curHeadSD,1);

    %data for swarmplot
    curData = curHeadSD(:,3) - curHeadSD(:,2);
    headSDdiff = [headSDdiff; curData];
    headSDdiffM = [headSDdiffM; mean(curData)];
    groupNum = [groupNum; repmat(type,length(curData),1)];
    colorMat = [colorMat; curColorMat];
    % statistic diff from zero
    [pfromZero, h] = signrank(curData, 0);
    psFromZero = [psFromZero , pfromZero];
    ns = [ns; n];
    Ns = [Ns; N];
    
end

fHL = figure;
swarmchart(groupNum,headSDdiff,10,colorMat,'filled','XJitterWidth',0.5);
hold on;
scatter(1:4, headSDdiffM,14,'k','Marker','+')
ylabel('Head Angle SD Diff -  (Stim-Pre)')
xticks(1:numType); xticklabels(stimType)
yline(0,'Color',[0.4 0.4 0.4],'LineStyle','--')

[pKW, tbl, stats] = kruskalwallis(headSDdiff,groupNum,'off');

if pKW < 0.05
    if isstring(groupNum)
        groupNum = cellstr(groupNum);
    end

    % Get unique group names
    [groupNames, ~, groupIdx] = unique(groupNum);
    numGroups = numel(groupNames);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = headSDdiff(groupIdx == i);
            data_j = headSDdiff(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [num2str(groupNames(i)) ' vs ' num2str(groupNames(j))];
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

end

%save figure #2
set(fHL,'PaperPosition',[1 5 3.5 2]);
fileName=[analysisFolder filesep 'headAngleSDalltypesDiffs'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% Figure 3F - movement diff all nights
close all
clearvars -except stimTable SA LMData animalsColors uniqueAnimals analysisFolder

load([analysisFolder filesep 'LMData.mat'])

stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
Groups = ["Wake","Pre","Stim"];

diffStimPre = [];
groupNum = [];
colorMat = [];
diffStimPreM =[];
ns = [];
Ns = [];

% Head movements:
fLMdiff=figure;
for type = 1:numType
    subplot(2,1,1)
    %plot the data
    curName = stimType(type);
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType) & ...
             ~contains(stimTable.Remarks,'Ex');
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    ns = [ns; n];
    Ns = [Ns; N];
    LMpre = LMData.LMpre(curTrials);
    LMwake = LMData.LMwake(curTrials);
    LMstimbinM = cell2mat(LMData.LMstimbin); % takes out the nan val
    LMstimbintrialM = mean(LMstimbinM(curTrials,:),2); % mean for each night
    
    curDiffStimPre = LMstimbintrialM-LMpre;
    diffStimPre = [diffStimPre;curDiffStimPre];
    diffStimPreM = [diffStimPreM; mean(curDiffStimPre)];
    groupNum = [groupNum; repmat(type,length(curDiffStimPre),1)];
    % curDiffPreWake = LMpre;
    % curDdiffStimWake = LMstimbintrialM;
    curData = [LMpre,LMstimbintrialM];
    % plot normelized to wake:

    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :); 
    colorMat = [colorMat; curColorMat];    

    hold on
    for i = 1:length(curData)
        plot([type*2-1, type*2],curData(i,:),'Color',curColorMat(i,:), 'Marker','.','LineWidth',0.5)
    end
    plot([type*2-1, type*2],mean(curData),'Color','k', 'Marker','.','LineWidth',1)

   
    % check the statistics:
    % Pairwise Wilcoxon signed-rank tests

    p = signrank(curData(:,1), curData(:,2));
    fprintf('%s: p-value = %.5f\n',stimType(type), p);

    % add titles. labels...
    if type==1
        ylabel('Head Movements')
    end
end
xlim([0.7 8.2])
xticks(1:8);xticklabels(repmat(["Pre","Stim"],1,numType))

% plot only stim-pre diff:
subplot(2,1,2)
swarmchart(groupNum,diffStimPre,15,colorMat,'filled','XJitterWidth',0.5);
hold on;
scatter(1:4,diffStimPreM,15,'k','Marker','+')
xticks(1:4), xticklabels(stimType); xlim([0.7 4.3]); %ylim([0 10]);
yline(0,'--','Color',[0.4 0.4 0.4])
[pKW, tbl, stats] = kruskalwallis(diffStimPre,groupNum,'off');

if pKW < 0.05
    if isstring(groupNum)
        groupNum = cellstr(groupNum);
    end

    % Get unique group names
    [groupNames, ~, groupIdx] = unique(groupNum);
    numGroups = numel(groupNames);

    % Initialize
    comparisons = {};
    raw_pvals = [];
    idx = 1;

    % Loop through all group pairs
    for i = 1:numGroups-1
        for j = i+1:numGroups
            % Extract data for group i and j
            data_i = diffStimPre(groupIdx == i);
            data_j = diffStimPre(groupIdx == j);

            % Wilcoxon rank-sum (Mann-Whitney U)
            [p, ~] = ranksum(data_i, data_j);

            % Store results
            comparisons{idx,1} = [num2str(groupNames(i)) ' vs ' num2str(groupNames(j))];
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

end

psFromZero = zeros(1,numType);
for i= 1:numType
    curData = diffStimPre(groupNum==i);
    [pfromZero, h] = signrank(curData, 0);
    psFromZero(1,i) = pfromZero; 
end

% savefigure
set(fLMdiff,'PaperPosition',[1 5 3.5 4]);
fileName=[analysisFolder filesep 'LMdiffs'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


