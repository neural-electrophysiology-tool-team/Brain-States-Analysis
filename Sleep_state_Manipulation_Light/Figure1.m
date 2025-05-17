%% Figure 1:
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
%% Figure 1 A, C, E
%% single traces from one night:

i = 22; % set the recording to PV161,N18
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;

% stimulations timings:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch}; %all stim timings
stimStartT = stims(1);
stimEndT = stims(end);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+200;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=130000;
ch = 17;

%% Figure 1A: plot trace from before stimulation:
p = 90*60*1000;
win = 180*1000;
tStart = firstTrig(1)-p+(870*1000);

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

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'blue';
ax.YAxis(1).Limits = [-500 500]; % Get left y-axis limits
ax.YAxis(2).Limits = [-800, 800];

%save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'cleanTraces'];
print(fileName,'-dpdf','-r300');

%% Figure 1C+E: plot trace + raster:

j = 16; %trial number example:
pTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
dbTmp=DB.bufferedDelta2BetaRatio(pTmp);
[lfp,lfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    
% add a raster plot according to stimulation
% Load spike times and cluster IDs (adjust file names as needed)
curPhyFoler = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort4'];
spikes=load([curPhyFoler filesep 'spike_data.mat']);

spikeTimes_ms = spikes.spike_times/(SA.currentDataObj.samplingFrequency(ch)/1000); %spike times in ms. (get it from phy in samples)
spikeClusters = spikes.spike_clusters; % cluster ID per spike
cluster_info = readtable([curPhyFoler filesep 'cluster_info.tsv'], 'FileType', 'text', 'Delimiter', '\t');

% Parameters for the plot
allClusters  = cluster_info.cluster_id(find(contains(cluster_info.group,'good')|contains(cluster_info.group,'mua')));
goodClusters = cluster_info.cluster_id(find(contains(cluster_info.group,'good')));
nClusters = length(allClusters);
nClustersg = length(goodClusters);

% plot trace + raster:
f = figure;
h1 = subplot(2,1,1);
plot(lfp_t/1000,squeeze(lfp),'k'); hold on;
curstims = (trial(j,:) -trial(j,1)) +pre ;
xline(curstims/1000,'r','LineWidth',1.5)
plot(dbTmp,'Color','b','LineWidth',2)
sgtitle(sprintf('Trial num: %i',j));
hold off;


% Plot RASTER!! spikes for each cluster, for each trig time!
    h2 = subplot(2,1,2);
    hold on;
    for l = 1:nClusters
        clusterID = allClusters(l);
        
        clusterSpikeTimes = spikeTimes_ms(spikeClusters == clusterID); % times for this unit
        curClusterSpikeT_ms=clusterSpikeTimes(clusterSpikeTimes>(firstTrig(j)-pre) & clusterSpikeTimes<(firstTrig(j)+post));
        TClusterSpikeT_ms = (curClusterSpikeT_ms-(firstTrig(j)-pre));
        % Plot each spike as a tick at y = cluster number
        for k = 1:length(TClusterSpikeT_ms)
            line([TClusterSpikeT_ms(k), TClusterSpikeT_ms(k)], [l - 0.4, l + 0.4], 'Color', 'k'); % tick mark
        end
        
    end
    curstims = trial(j,:) -trial(j,1) +pre ;
   
    xline(curstims,'r','LineWidth',1.5)
    box on
    % Convert x-axis labels from ms to s by setting the x-axis ticks and labels
    xticks = get(gca, 'XTick');          % Get current x-axis tick values in ms
    set(gca, 'XTick', xticks);           % Set the same ticks
    set(gca, 'XTickLabel', xticks / 1000); % Display tick labels in seconds
    xlabel('Time (s)');
    % xlabel('Time (ms)');
    ylabel('Neuron/Unit');
    ylim([0.4,nClusters+0.5])
    title('Raster Plot');
    hold off;

% save figure:
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'singleTrialRasterPV161N18t16'];
print(fileName,'-depsc','-vector');

%% Figure 1F

% Plot Spike rates avrages for all nights:
    % 3 subplots: 1 unit, all units same night, all units
    % night for the single unit:
    i = 22;
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    curPhyFoler = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort4'];
    spikeRateFile = [curPhyFoler filesep 'spikeRateAll.mat'];
    load(spikeRateFile,"spikeRateAll","spikeRateT","allClusters","OL");
    % all nights, all units:
    allnightsfilename = [analysisFolder filesep 'allNightsSpikingRate.mat'];
    load(allnightsfilename,"AllNightsUnits","spikeRecs","xPositions")

unitMean = mean(spikeRateAll,3);
bestU = 48;
bestUind = find(allClusters==bestU);
% stimDiff = mean(mean(diff(trial,[],2)));
% xPositions = (pre + (0:7)*(stimDiff))/OL;  % Example positions for lines

f= figure;
%figrue1: one unit from 1 nights/
h1 = subplot(3,1,1);
plot(spikeRateT/OL,unitMean(bestUind,:),'Color','k', 'LineWidth',1)
hold on
xline(xPositions,'r','LineWidth',1)
ylabel('Spike/s')
title('one unit mean across all trials')
hold off

%figure 2: all units one night: 
h2=subplot(3,1,2);
meanOneNight = mean(unitMean,1);
n = height(unitMean);
title 'Mean all units - One Night')
plot(spikeRateT/OL,meanOneNight,'Color','k', 'LineWidth',1); hold on;
xline(xPositions,'r','LineWidth',1)
ylabel('Spike/s')
hold off
annotation('textbox', [0.8, 0.5, 0.03, 0.1], 'String', ...
    sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


%figure 3: all units - all nights
h3=subplot(3,1,3);
title('Mean all units, all Red Nights')
plot(spikeRateT/OL,mean(AllNightsUnits,1),'Color','k', 'LineWidth',1); hold on;
xline(xPositions,'r','LineWidth',1)
xlabel('Time[s]'); 
ylabel('Spike/s')
n= height(AllNightsUnits);
N = numel(spikeRecs);
annotation('textbox', [0.8, 0.2, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% linkaxes([h1,h2,h3],'x')
%savefigure
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'spikeRates-allnights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Figure 1 G:
% load data
ITIISImatfilename = [analysisFolder filesep 'ITIISIdata.mat'];
load(ITIISImatfilename, "ITIsamples", "ISIsamples","spikeRecs")


fall = figure;
subplot(1,2,1)
IIdata = [ITIsamples, ISIsamples];
plot([1,2] ,IIdata,'Color',[0.7 0.7 0.7],'Marker','.','MarkerSize',4);
hold on;
plot([1,2] ,mean(IIdata),'Color','k', 'LineWidth',2,'Marker','.','MarkerSize',4);
xticks([1,2]);xticklabels(["Before Stims", "After Stims"]);xlim([0.5 2.5])
% grid on;
ylabel('Spikes/S')
sgtitle ('all units')
n = length(IIdata); 
N= length(spikeRecs);

subplot(1,2,2)
diff= ITIsamples- ISIsamples; %ISI = after, ITI =before
x = [ones(length(IIdata),1)];
colors = [ 0.5 0.7 0.8];
swarmchart(x,diff,1,colors,'filled','XJitterWidth',0.8);
ylabel('Spikes/S')
xticks(1);xticklabels("Before-After");
xlim([0.5 1.5])
yline(0,'--','Color',[0.4 0.4 0.4])
percent = (sum(diff<0)/length(diff))*100; % how many units are bellow zero in percent
[pWilcoxon, ~, statsWilcoxon] = signrank(ITIsamples, ISIsamples);

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


    %%