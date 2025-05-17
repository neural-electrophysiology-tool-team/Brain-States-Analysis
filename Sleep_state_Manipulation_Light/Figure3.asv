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

%% Figure 3D - D/B decrease

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
%plot Period Times:


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



%%
ns= [];
for type = 1:numType
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,"Ex");
    ns =[ns; sum(curTrials)];
end