% head angles : diffs:

%% prepare head angle data for plots:
load([analysisFolder filesep 'HeadAngleAvg.mat'])
load([analysisFolder filesep 'HeadAngleSD.mat'])

headAngDiff = diff(HeadAngleAvg,[],2);
% set the zero to 90 Deg, according to accelerometer data ( this is the z
% axis, when it is 90 the accelerometer is penpendicular to the ground)
HeadAngleAvgP = HeadAngleAvg -90;
%% plot SDs
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
% numAnimal = length(animals);
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
    h = subplot(2,2,type);

    %create data subset:
    curName = stimType(type);
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType)& ~contains(stimTable.Remarks,'Ex' )...
            & (headAngDiff(:,1)>3 | headAngDiff(:,1)<-3); 
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :); 
    % curHeadAvg = HeadAngleAvgP(curTrials,:);
    % curHeadAvgmean = mean(curHeadAvg,1,'omitnan');

    curHeadSD = headAngleSD(curTrials,1:3);
    curHeadSDmean = mean(curHeadSD,1,'omitnan');

    %data for swarmplot
    curData = curHeadSD(:,3) - curHeadSD(:,2);
    headSDdiff = [headSDdiff; curData];
    headSDdiffM = [headSDdiffM; mean(curData)];
    groupNum = [groupNum; repmat(type,length(curData),1)];
    colorMat = [colorMat; curColorMat];
    [pfromZero, h] = signrank(curData, 0);
    psFromZero = [psFromZero , pfromZero];
    ns = [ns; n];
    Ns = [Ns; N];
    %statistics:

    [p, tbl, stats] = friedman(curHeadSD, 1,'off'); % Here, 1 indicates within-subjects design
    fprintf('p-value for freidman ANOVA test: %.5f\n',p)
    % % p-valure is very low, post hoc:
    % data for the four groups
    wakeAng = curHeadSD(:,1);
    beforeAng = curHeadSD(:,2);
    duringAng = curHeadSD(:,3);
    % afterAng = curHeadSD(:,4);

    % Bonferroni-corrected alpha level
    alpha = 0.05 / 3;

    % Pairwise Wilcoxon signed-rank tests
    [p_wake_pre, ~, stats_wake_before] = signrank(wakeAng, beforeAng);
    [p_pre_during, ~, stats_before_during] = signrank(beforeAng, duringAng);
    % [p_during_after, ~, stats_during_after] = signrank(duringAng, afterAng);
    [p_wake_during, ~, stats_wake_during] = signrank(wakeAng, duringAng);

    % Display results with Bonferroni correction
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Wake vs Before: p-value = %.4f (Significant if < %.4f)\n', p_wake_pre, alpha);
    fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_pre_during, alpha);
    % fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
    fprintf('Wake vs During: p-value = %.4f (Significant if < %.4f)\n', p_wake_during, alpha);

    %save statistics:
    statsHeadangSd.(curName).alpha = alpha;
    statsHeadangSd.(curName).pAnova = p;
    statsHeadangSd.(curName).p_wake_pre = p_wake_pre;
    statsHeadangSd.(curName).p_wake_during = p_wake_during;
    statsHeadangSd.(curName).p_pre_during = p_pre_during;
    % statsHeadangSd.(curName).p_during_after = p_during_after;

    %plot the data

    if n>0
        x1 = 1:3;
        for j = 1:height(curHeadSD)
            plot(x1,curHeadSD(j,:),'Color',curColorMat(j,:),'Marker','.'); hold on;
        end
        plot(x1, mean(curHeadSD), 'Color','k','Marker','.','LineWidth',1.5)
        xlim([0.7,3.2]); xticks(x1); xticklabels(["Wake","Sleep","Stim"])
        ylabel('Avg Head Angles SD')
        % yline(0,'--','Headstage penpendicular to floor')
    end

    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   % ylim([-30 50])
    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('head angles from ground')
    title(stimType(type))
 
end
% savefigure
set(fHA,'Position',[50 50 1000 720]);
fileName=[analysisFolder filesep 'headLiftsalltypesSD'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)],'-bestfit');
save([analysisFolder filesep 'headLiftsalltypesSD.mat'], "statsHeadangSd")


fHL = figure;
swarmchart(groupNum,headSDdiff,10,colorMat,'filled','XJitterWidth',0.5);
hold on;
scatter(1:4, headSDdiffM,14,'k','Marker','+')
ylabel('Head Angle SD Diff -  (Stim-Pre)')
xticks(1:numType); xticklabels(stimType)
yline(0,'Color',[0.4 0.4 0.4],'LineStyle','--')
% ylim([-350 150])

%save figure #2
set(fHL,'PaperPosition',[1 5 3.5 2]);
fileName=[analysisFolder filesep 'headAngleSDalltypesDiffs'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% diffs: 
diff32 = HeadAngleAvgP(:,3)-HeadAngleAvgP(:,2);


stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);

diffData1 = [];
groupNum = [];
colorMat = [];
ns = [];
Ns= [];

for type = 1:numType

    %create data subset:
    curName = stimType(type);
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType)& ~contains(stimTable.Remarks,'Ex' )...
            & (headAngDiff(:,1)>3 | headAngDiff(:,1)<-3); 
    curData = diff32(curTrials);
    diffData1 = [diffData1; curData];
    groupNum = [groupNum; repmat(type,length(curData),1)];
    ns =[ns; length(curData)];
    Ns = [Ns;length(unique(stimTable.Animal(curTrials)))];

    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :); 
    colorMat = [colorMat; curColorMat];
end

% statistics: 
% 1. kruskal wallas:
[p, tbl, stats] = kruskalwallis(diffData1,groupNum,'off');

if p < 0.05
    c = multcompare(stats, 'CType', 'dunn-sidak');
end

%plot the data
fHA=figure;
swarmchart(groupNum,diffData1,10,colorMat,"filled")
xticks(1:4); xticklabels(stimType);
ylabel('Head Angle diff ')
title('Head angle diff - (stim-pre)')

% savefigure
set(fHA,'PaperPosition',[1 1 4 2.5]);
fileName=[analysisFolder filesep 'headAngDiffStim_Pre'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)],'-bestfit');

%% diffs - to wake:

diff31 = HeadAngleAvgP(:,3)-HeadAngleAvgP(:,1);
diff21 = HeadAngleAvgP(:,2)-HeadAngleAvgP(:,1);

stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
% 
diffData2 = [];
% groupNum = [];
% colorMat = [];
% ns = [];
% Ns= [];
figure;
for type = 1:numType

    %create data subset:
    curName = stimType(type);
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType)& ~contains(stimTable.Remarks,'Ex' )...
            & (headAngDiff(:,1)>3 | headAngDiff(:,1)<-3); 
    curData = [diff21(curTrials) , diff31(curTrials)];
    diffData2 = [diffData2; curData];
    % groupNum = [groupNum; repmat([type*2-1, type*2],length(curData),1)];
    ns =[ns; length(curData)];
    Ns = [Ns;length(unique(stimTable.Animal(curTrials)))];

    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :); 
    % colorMat = [colorMat; curColorMat];
    hold on
    for i = 1:length(curData)
        plot([type*2-1, type*2],curData(i,:),'Color',curColorMat(i,:), 'Marker','.','LineWidth',0.5)
    end
    plot([type*2-1, type*2],mean(curData),'Color','k', 'Marker','.','LineWidth',1)
    % yline(0,'--','Color',[0.4 0.4 0.4])


% wilcoxon:
p = signrank(curData(:,1), curData(:,2));
fprintf('%s: p-value = %.5f\n',stimType(type), p);


end
ylabel('Normelized head angle (to wake)')
xlim([0.7 8.2])
    xticks(1:8);xticklabels(repmat(["Pre","Stim"],1,numType))

% %plot the data
% fHA=figure;
% swarmchart(groupNum,diffData2,10,colorMat,"filled")
% % xticks(1:4); xticklabels(stimType);
% ylabel('Head Angle diff ')
% title('Head angle diff - (stim-pre)')
% 
% savefigure
set(fHA,'PaperPosition',[1 1 4 2.5]);
fileName=[analysisFolder filesep 'headAngDiffStim_Pre_toWAKE'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
