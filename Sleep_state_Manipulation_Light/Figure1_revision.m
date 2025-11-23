%% Figure 1:
% initiating stimTable and other essentials. 
SA=sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTableAll.mat'])
load([analysisFolder filesep 'LMDataAll.mat'])
animalsColors = [
    0,0.474509803921569,0.549019607843137; % dark blue  - PV106
    255/255, 142/255, 71/255;% HEX:  FF8E47 - orange  - PV126
    0.819607843137255,0.286274509803922,0.356862745098039; %red-ish - PV143
    28/255, 144/255, 217/255;  % HEX: 1C90D9 - blue - PV149
    0.929411764705882,0.682352941176471,0.286274509803922; %yellow-ish - PV153
    148/255, 120/255, 186/255; % HEX: 9478BA - perpule - PV157
    217/255, 62/255, 67/255; % HEX: D93E43 - red - PV159
    255/255, 202/255, 58/255; % HEX: FFCA3A - yellow -  PV161
    97/255, 184/255, 68/255;  % HEX:61B844 - Green -PV162
    
];
uniqueAnimals = unique(stimTable.Animal);
%% Figure 1 A, C, E
% parameters calculations for next plots:
% single traces from one night:

i = 82; % set the recording to PV153,N11
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;

% stimulations timings:
stims = SA.getStimTriggers;%all stim timings
stimStartT = stims(1);
stimEndT = stims(end);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+200;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=130000;
ch = 17;

%% Figure 1A
% plot trace from before stimulation:
p = 2*60*60*1000 +(10*1000);
win = 180*1000;
tStart = firstTrig(1)-p +(590*1000);

pTmpC = find(DB.t_ms>(tStart) & DB.t_ms<(tStart+win));
dbTmpC = DB.bufferedDelta2BetaRatio(pTmpC);
[lfpC,lfp_tC] = SA.currentDataObj.getData(ch,tStart,win);
% plot:
f = figure;
set(f, 'Position', [100, 100, 600, 100]);
hold on
yyaxis left
plot(lfp_tC/1000,squeeze(lfpC),'k'); hold on;
ylabel('microvolts')
yyaxis right
plot(dbTmpC,'Color','b','LineWidth',2)
ylabel('\delta/\beta ')
hold off
xlabel('Time[s]')

ax2 = gca;
ax2.YAxis(1).Color = 'k';
ax2.YAxis(2).Color = 'blue';
ax2.YAxis(1).Limits = [-500 500]; % Get left y-axis limits
ax2.YAxis(2).Limits = [-300, 300];

%save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'cleanTracesPV153N15'];
print(fileName,'-dpdf','-r300');

%% Figure 1C+E
% plot trace + raster:
j = 19; %trial number example:

pTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
dbTmp=DB.bufferedDelta2BetaRatio(pTmp);
[lfp,lfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    
% get rasterdata
curPhyFoler = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort'];
tIc = SA.currentDataObj.convertPhySorting2tIc(curPhyFoler);
bin=100; %ms
startTime = firstTrig(j)-pre;
width = pre+post; %150 sec in ms.
[M]=BuildBurstMatrix(tIc.ic,round(tIc.t/bin),round(startTime/bin),round(width/bin));

% plot trace + raster:
f = figure;
ax1=subplot(2,1,1);
t = linspace(-pre/1000,post/1000,length(lfp_t));
plot(t,squeeze(lfp),'k'); hold on;
curstims = (trial(j,:) -trial(j,1));
xline(curstims/1000,'r','LineWidth',1.5);
tdb = linspace(-pre/1000,post/1000,length(dbTmp));
yyaxis right
plot(tdb,dbTmp,'Color','b','LineWidth',2);
ylim([-300 300]);                   % limit for left axis
ax = gca;
ax.YColor = 'b';
sgtitle(sprintf('Trial num: %i',j));
hold off;

% Plot RASTER!! spikes for each cluster, for each trig time!
ax2=subplot(2,1,2);
hold on;
spikemat = squeeze(M);
%normelization:
Fs = 1000/bin; % sampling rate (150000 samples = 150 seconds ? Fs = 1000 Hz)

nCols = size(spikemat,2);
tEdges = linspace(-pre/1000,post/1000,nCols+1);
h2 = imagesc(tEdges, 1:size(spikemat,1), spikemat); box on;
colormap(flipud(gray));   % flip so 0 is white, max is black
axis xy;
clim([0 max(spikemat(:))]);%box on;
% clim([0 1]);
% colorbar;
set(h2, 'Interpolation','nearest');   % ?? Key line: prevents blending between pixels

xlabel('Time (s)');
ylabel('unit #');
% ploting the stimulatiis
curstims = trial(j,:) -trial(j,1);
xline(curstims/1000,'r','LineWidth',1.5)
linkaxes([ax1,ax2],'x')
xlim([-pre/1000-0.5,post/1000]);
ylim([0.5 size(spikemat,1)+0.7])

% % save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'singleTrialRasterPV153N33t32'];
print(fileName,'-depsc','-vector');

%% Figure 1D
i = 90; %PV153, N31
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
stimShamFile=[SA.currentAnalysisFolder filesep  'stimSham.mat'];
load(stimShamFile); %load data
    
trialsham = height(StimDBSham)-70+1;
    
colorLim=[0 600];
f=figure;
subplot(4,2,[1:2:6]);
imagesc(StimDBSham(trialsham:end,:),colorLim);
ylabel('Trial #');title('Sham');hold on;
set(gca,'XTick',[]);
cb=colorbar('Position',[0.47 0.76 0.013 0.17]);
ylabel(cb,'\delta/\beta');
line([pre/1000 pre/1000],ylim,'color','r');
line([(pre+stimDuration)/1000 (pre+stimDuration)/1000],ylim,'color','r');
subplot(4,2,7);plot(ts-pre/1000,mean(StimDBSham,'omitnan'),'k');% plot mean
xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/2);xlim([-pre/1000 post/1000])
line([0 0],ylim,'color','r');
line([stimDuration/1000 stimDuration/1000],ylim,'color','r');
subplot(4,2,[2:2:6]);
imagesc(StimDB,colorLim);ylabel('Trial #');
title('Stim');set(gca,'XTick',[]);
cb=colorbar('Position',[ 0.91 0.76 0.013 0.17]);ylabel(cb,'\delta/\beta');
line([pre/1000 pre/1000],ylim,'color','r');
line([(pre+stimDuration)/1000 (pre+stimDuration)/1000],ylim,'color','r');
subplot(4,2,8);plot(ts-pre/1000,mean(StimDB,'omitnan'),'k');
xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/2);
line([0 0],ylim,'color','r');
line([stimDuration/1000 stimDuration/1000],ylim,'color','r');
xlim([-pre/1000 post/1000])

%   savefigure
fileName=[analysisFolder filesep 'stim_sham_activationPV153N31.pdf'];
saveas (f, fileName);

%% Figure 1F

% Plot Spike rates avrages for all nights:
% 2 subplots: 1 unit, all units
% night for the single unit:
i = 56;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
tIcPath= [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort'];
tIc = SA.currentDataObj.convertPhySorting2tIc(tIcPath);

% parameters:
bin=1000;
pre = 20*1000; %ms
width = 150*1000; %150 sec in ms.
stims = SA.getStimTriggers;
trial = reshape(stims,[8,length(stims)/8])';
firstTrig=stims(1:8:end-2);
startTimes = round(firstTrig-pre);
curstims = ((trial(1,:) -trial(1,1)))/1000;
[M]=BuildBurstMatrix(tIc.ic,round(tIc.t/bin),round(startTimes/bin),round(width/bin));
meanT = squeeze(mean(M,1,'omitnan'));
% all nights, all units:
meanAllRecsPath = [analysisFolder filesep 'meanSpikeRate1000msbin.mat'];
load(meanAllRecsPath)
%% get unit annotation
labelsAll = [];
spikesInd = stimTable.spikes ==1& (strcmp(stimTable.Remarks,'white')|contains(stimTable.Remarks,'DayTime'));
animals =stimTable.Animal(spikesInd);recnames = stimTable.recNames(spikesInd);
recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);

for i = 1:length(recList)
    SA.setCurrentRecording(recList{i});
    tIcPath = [SA.currentDataObj.recordingDir '/spikeSorting/kilosort'];
    tIc = SA.currentDataObj.convertPhySorting2tIc(tIcPath);
    curLabels=tIc.label;
    labelsAll = [labelsAll; curLabels];
end%

%% plot

f= figure;
%figrue1: one unit from 1 nights/
ax1 = subplot(2,1,1);
t=linspace(-pre/1000,(width-pre)/1000,size(meanAllRecs,2));
plot(t,meanT(4,:),'Color','k', 'LineWidth',1)
hold on
xline(curstims,'r','LineWidth',1)
ylabel('Spike/s')
title('one unit mean across all trials')
hold off


% %figure 3: all units - all nights
% ax2=subplot(2,1,2);
% title('Mean all units, all Red Nights')
% plot(t,mean(meanAllRecs,1,'omitnan'),'Color','k', 'LineWidth',1); hold on;
% xline(curstims,'r','LineWidth',1)
% xlabel('Time[s]'); 
% ylabel('Spike/s')
% n= size(meanAllRecs,1);
% spikesInd = stimTable.spikes ==1& (strcmp(stimTable.Remarks,'white')|contains(stimTable.Remarks,'DayTime'));
% N = sum(spikesInd);
% annotation('textbox', [0.8, 0.2, 0.03, 0.1], 'String', ...
%     sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');
% 
% % 
%figure 3: good units - all nights
ax2=subplot(2,1,2);
title('Mean good units, all white Nights')
data = meanAllRecs(goodUnitsInd,:);
plot(t,mean(data,1,'omitnan'),'Color','k', 'LineWidth',1); hold on;
xline(curstims,'r','LineWidth',1)
xlabel('Time[s]'); 
ylabel('Spike/s')
n= size(data,1);
spikesInd = stimTable.spikes ==1& (strcmp(stimTable.Remarks,'white')|contains(stimTable.Remarks,'DayTime'));
N = sum(spikesInd);
annotation('textbox', [0.8, 0.2, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


linkaxes([ax1,ax2],'x'); xlim([-20 130])
% %savefigure
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'spikeRates-allnights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Figure 1G
%% get the pre-post stim spike rates.
clearvars -except stimTable SA LMData meanAllRecs analysisFolder 
meanAllRecsPath = [analysisFolder filesep 'meanSpikeRate1000msbin.mat'];
load(meanAllRecsPath)

% spiking recs:
spikesInd = stimTable.spikes ==1& (strcmp(stimTable.Remarks,'white')|contains(stimTable.Remarks,'DayTime'));
animals =stimTable.Animal(spikesInd);recnames = stimTable.recNames(spikesInd);
recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
SA.setCurrentRecording(recList{1});
stims = SA.getStimTriggers;
trial = reshape(stims,[8,length(stims)/8])';
curstims = (trial(1,:) -trial(1,1))/1000;
%parameters:
bin=1000; 
pre = 20*1000; %ms
width = 150*1000; %150 sec in ms.
spikeRateT = linspace(-pre/1000,(width-pre)/1000,width/bin);

% timings for pre+post
preStimI = find(spikeRateT> 9.5 & spikeRateT< 19.500);
postStimI = find(spikeRateT> (curstims(end)+2) & spikeRateT< (curstims(end)+12));

% cal data:
preData=mean(meanAllRecs(:,preStimI),2);
postData=mean(meanAllRecs(:,postStimI),2);

% %% plot
% 
% fall = figure;
% subplot(1,2,1)
% IIdata = [preData, postData];
% plot([1,2] ,IIdata,'Color',[0.7 0.7 0.7],'Marker','.','MarkerSize',4);
% hold on;
% plot([1,2] ,mean(IIdata),'Color','k', 'LineWidth',2,'Marker','.','MarkerSize',4);
% xticks([1,2]);xticklabels(["Before Stims", "After Stims"]);xlim([0.7 2.3])
% % grid on;
% ylabel('Spikes/S')
% sgtitle ('all units')
% n = length(IIdata); 
% N= length(recList);
% 
% subplot(1,2,2)
% diff= preData- postData; %ISI = after, ITI =before
% x = [ones(length(IIdata),1)];
% colors = [ 0.5 0.7 0.8];
% swarmchart(x,diff,10,colors,'filled','XJitterWidth',0.8);
% ylabel('Spikes/S')
% xticks(1);xticklabels("Before-After");
% xlim([0.5 1.5])
% yline(0,'--','Color',[0.4 0.4 0.4])
% percent = (sum(diff<0)/length(diff))*100; % how many units are bellow zero in percent
% [pWilcoxon, ~, statsWilcoxon] = signrank(preData, postData);
% 
% annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
%     sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');
% annotation('textbox', [0.8, 0.7, 0.03, 0.1], 'String', ...
%     sprintf('Wilcoxon p-val: %.4f\n', pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');
% % avefigures
%     set(fall,'PaperPosition',[1 1 2.1 2]);
%     fileName=[analysisFolder filesep 'ITIISIallunits'];
%     print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
% 



%% plot only good units:
goodUnitsInd = find(contains(labelsAll,'good'));
preDataG = preData(goodUnitsInd);
postDataG= postData(goodUnitsInd);
fall = figure;
subplot(1,2,1)
IIdata = [preDataG, postDataG];
plot([1,2] ,IIdata,'Color',[0.7 0.7 0.7],'Marker','.','MarkerSize',4);
hold on;
plot([1,2] ,mean(IIdata),'Color','k', 'LineWidth',2,'Marker','.','MarkerSize',4);
xticks([1,2]);xticklabels(["Before Stims", "After Stims"]);xlim([0.7 2.3]);

% grid on;
ylabel('Spikes/S')
sgtitle ('agood units')
n = length(IIdata); 
N= length(recList);

subplot(1,2,2)
diff= preDataG- postDataG; %ISI = after, ITI =before
x = [ones(length(IIdata),1)];
colors = [ 0.5 0.7 0.8];
swarmchart(x,diff,10,colors,'filled','XJitterWidth',0.8);
ylabel('Spikes/S')
xticks(1);xticklabels("Before-After");
xlim([0.5 1.5])
yline(0,'--','Color',[0.4 0.4 0.4])
%stats:
percent = (sum(diff<0)/length(diff))*100; % how many units are bellow zero in percent
[pWilcoxon, ~, statsWilcoxon] = signrank(preDataG, postDataG)

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.8, 0.7, 0.03, 0.1], 'String', ...
    sprintf('Wilcoxon p-val: %.4f\n', pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
% avefigures
    set(fall,'PaperPosition',[1 1 2.1 2]);
    fileName=[analysisFolder filesep 'ITIISIallunits'];
    print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);



%% Figure 1H - stim sham all red nights:

% wavelength =;
curTrials = (contains(stimTable.Remarks,'white')|contains(stimTable.Remarks,'DayTime'))& ...
    ~contains(stimTable.Remarks,'Ex') & ...
    all(~isnan(stimTable.dbDiffStimM), 2) &...
    all(~isnan(stimTable.dbDiffShamM), 2); 
    
n = sum(curTrials);
animals = unique(stimTable.Animal(curTrials));
N = length(animals);

% statistical tests for that figure 
groupSham = stimTable.dbDiffShamM(curTrials);
groupStim = stimTable.dbDiffStimM(curTrials);

% 1. Wilcoxon Signed-Rank Test
[pWilcoxon, ~, statsWilcoxon] = signrank(groupSham, groupStim);
fprintf('Wilcoxon Signed-Rank Test p-value: %.4f\n', pWilcoxon);
disp(statsWilcoxon)

% plot:
figure;
title('Change in mean D/B norm across trials - All nights')
x=[1,2];
curData = [stimTable.dbDiffShamM(curTrials), stimTable.dbDiffStimM(curTrials)];
% color code per animal:
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); % Use the indices to directly fetch corresponding colors
curMeanStim = mean(stimTable.dbDiffStimM(curTrials),1);
curMeanSham = mean(stimTable.dbDiffShamM(curTrials),1);
hold on
for i = 1:height(curData)
    plot(x, curData(i,:), 'Color',curColorMat(i,:), 'Marker','.', 'MarkerSize',10, 'LineWidth',1)
end
plot(x,[curMeanSham,curMeanStim],'color','k','LineWidth',2,'Marker','.','MarkerSize',10)
hold off


% ylim([-40 120])
xlim([0.9 2.1]);
xticks([1, 2]); % Position of the x-ticks
xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks

%save figure:
set(gcf,'PaperPosition',[1 4 1.2 1.6])
saveas (gcf, [analysisFolder filesep 'DBdiffStimShamRedNights.pdf']);

    
%% Figure 1I
% plot sliding AC sith stimulations: (bottom)
%set the recording:
i = 74; %PV106, N32 (74)
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
AC = SA.getDelta2BetaAC;
SA.plotDelta2BetaSlidingAC ('stim',1);

% plot partial AC for each part of the night:
ACpre = stimTable.ACpre{i};
ACstim = stimTable.ACstim{i};
ACpost = stimTable.ACpost{i};

PDPcolors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];
ACstructs = {ACpre,ACstim,ACpost};
labels = {'preStim', 'Stim', 'postStim'};
for j = 1:3
    struct2vars(ACstructs{j});
    fAC = figure;
    %plot:
    lineHandles = plot(autocorrTimes/1000,real(xcf),'Color',PDPcolors(j,:),'LineWidth',4);
    ylim([-0.5 1])
    set(lineHandles(1),'MarkerSize',4);
    % grid('on');
    xlabel('Period [s]');
    ylabel('Auto corr.');
    hold('on');

    plot(period/1000,real(xcf(pPeriod)),'o','MarkerSize',5,'color','k');
    text(period/1000,0.05+real(xcf(pPeriod)),num2str(period/1000));

    a = axis;
    plot([a(1) a(2)],[0 0],'-k');
    hold('off');
    title (labels{j})
   % save fig:
   set(fAC,'PaperPositionMode','auto');
   fileName=[analysisFolder filesep 'dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(round(parDbAutocorr.tStart)) labels{j}];
   print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

end

%% Figure 1J
% plot AC - only white nights:
curTrials = (contains(stimTable.Remarks,'white')|contains(stimTable.Remarks,'DayTime')) & ...
            ~contains(stimTable.Remarks,'Ex') & ...
            all(~isnan(stimTable.ACcomPer),2) &...
            all(stimTable.ACcomP2V > 0.15,2);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

% statistical tests:
% first, we need to check the differences in genral, in Freidman test,
% which is a a-parametrical ANOVA test. then we can use wilcoxon post-hoc
% to check where is the differenc (with benforoni corection)

% Assuming data in columns where each row is a subject and each column is a timepoint
[p, tbl, stats] = friedman(stimTable.ACcomPer(curTrials,:), 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
if p<0.05
    before = stimTable.ACcomPer(curTrials,1);
    during = stimTable.ACcomPer(curTrials,2);
    after = stimTable.ACcomPer(curTrials,3);

    % Pairwise Wilcoxon signed-rank tests
    [p_before_during, ~, stats_before_during] = signrank(before, during);
    [p_during_after, ~, stats_during_after] = signrank(during, after);
    [p_after_before, ~, stats_after_before] = signrank(after, before);

    raw_pvals = [p_before_during,p_during_after,p_after_before];
    num_comparisons = 3;
    % Display results with Bonferroni correction
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Before vs During: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('During vs After: p-value = %.4f\n', corrected_pvals_bonferroni(2));
    fprintf('After vs Before: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
end

%plot
figure;
x = 1:3;
% color code per animal:
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
curMean = mean(stimTable.ACcomPer(curTrials,:)/1000,1);
curData = stimTable.ACcomPer(curTrials,:)/1000;
hold on
for i = 1:height(curData)
    plot(x, curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10,'LineWidth',1)
end
plot(x, curMean,'Color','k','LineWidth',2)
xlim([0.75 3.25])
xticks(x)  % Set ticks after xlim to avoid automatic adjustment
xticklabels({'Pre', 'During', 'Post'})

ylim([50 200])
ylabel('Period Time[s]')
title ('Perios Times changes - all white nights')
yline(156,'--',Color=[0.5 0.5 0.5])

% savefigure
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'ACperiodwhites.pdf']);