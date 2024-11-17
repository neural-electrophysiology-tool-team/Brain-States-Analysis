%% light manipulation analysis 08/05/2024
% This is the new version. 

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
analysisFolder = '/media/sil3/Data/Pogona_Vitticeps/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTable.mat'])
load([analysisFolder filesep 'LMdata.mat'])

%% analysis folder
% analysisFolder = '/media/sil3/Data/Pogona_Vitticeps/NitzanAnalysisFiles';
% SA.batchProcessData('getDelta2BetaRatio',{})

%% get the stim sham for single rec
SA.setCurrentRecording('Animal=PV162,recNames=Night23');
getStimSham(SA,11);
SA.getDelta2BetaRatio('tStart', 0 ,'win',1000*60*60*11);
SA.getDelta2BetaAC;
SA.plotDelta2BetaRatio
SA.plotDelta2BetaAC;
SA.plotDelta2BetaSlidingAC;
plotStimSham(SA)

%% creating the StimTable:

%% get all the stim sham avg from all recs and put in a new table
% this part goes over all the records in SA.
%  for every recoerd that is tagged (1/2/3..) 
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
% SA.setCurrentRecording('Animal=PV162,recNames=Night27');
maniRecs = SA.recTable.Mani>0; % taking all the rows with manipulation
stimTable = SA.recTable(maniRecs,{'Animal','recNames','Remarks','Mani','LizMov','StimTrighCh'});  % creating new table
stimTable.StimAvg = cell(height(stimTable),1);
stimTable.StimAvgSham = cell(height(stimTable),1);
stimTable.times = cell(height(stimTable),1);
stimTable.stimDuration = zeros(height(stimTable),1);
stimTable.ACpre = cell(height(stimTable),1);
stimTable.ACstim = cell(height(stimTable),1);
stimTable.ACpost = cell(height(stimTable),1);

%% geting the AND AC stim sham data from all records 

for i = 1:height(stimTable)

    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    % run all required analysis:
    SA.getDelta2BetaRatio;
    SA.getDelta2BetaAC;
    SA.getDigitalTriggers
    % get and save the stim avg
    s = getStimSham(SA,stimTable.StimTrighCh(i),1);
    disp('got Stim for this rec')
    stimTable.StimAvg(i) = {mean(s.StimDB,1)};
    stimTable.StimAvgSham(i) = {mean(s.StimDBSham,1)};
    stimTable.times(i) = {s.ts};
    stimTable.stimDuration(i) = s.stimDur;
    disp('stimsham in table')

    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimStartT = T.tTrig{t_ch}(1);
    stimEndT = T.tTrig{t_ch}(end);
    firstTrig=T.tTrig{t_ch}(1:8:end-2);
    endStim=T.tTrig{t_ch}(8:8:end)+400;
%     stimDuration=(endStim(1)-firstTrig(1));

     % calculate the AC and the P2V for each part of the stimulation
    p = 30*60*1000; %some time diff for the cycle to change.
    ACwin = 2*60*60*1000; % 2 hrs in ms

    SA.getDelta2BetaAC('tStart',stimStartT-ACwin, 'win',ACwin,'overwrite',1);
    stimTable.ACpre(i) = {SA.getDelta2BetaAC('tStart',stimStartT-ACwin, 'win',ACwin)};
    SA.getDelta2BetaAC('tStart',stimStartT+p, 'win',ACwin,'overwrite',1);
    stimTable.ACstim(i) = {SA.getDelta2BetaAC('tStart',stimStartT+p, 'win',ACwin)};
    SA.getDelta2BetaAC('tStart',stimEndT+p,'win',ACwin,'overwrite',1);
    stimTable.ACpost(i) = {SA.getDelta2BetaAC('tStart',stimEndT+p,'win',ACwin)};
    disp('AC in table')
    
      
end
clear recName
clear s
%% Getting the liz mov for all animals + add to table
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
stimTable.LM_DBt = cell(height(stimTable),1);
for i = 1:height(stimTable)
    % set te current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    % run all required analysis:
    if stimTable.LizMov(i) ==1
        SA.getLizardMovements
        LM = SA.getLizardMovements;

        DB = SA.getDelta2BetaRatio;
        

        % calculate the number of movements to each DB bin
        LM_DBt = zeros(size(DB.t_ms));
        % Loop through each bin in DB and count the events in LM that fall within each bin
        for j = 1:length(DB.t_ms)-1
            % Count events from LM that fall within the current bin (DB(i) to DB(i+1))
            LM_DBt(j) = sum(LM.t_mov_ms >= DB.t_ms(j)& LM.t_mov_ms < DB.t_ms(j+1));
        end
        % Count any events at the last bin edge
        LM_DBt(end) = sum(LM.t_mov_ms >= DB.t_ms(end));

        %put in stimTable:
        stimTable.LM_DBt(i) = {LM_DBt};
        disp('LM in stimTabl')
    end
end


%% save stimTable
save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');

%% Figure 1 - traces+ spikes

%% look at single traces from one night:

i = 30;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
% recName2 = 'Animal=PV157,recNames=Night37';
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;
% SA.getDelta2BetaRatio
% SA.plotDelta2BetaRatio

% stimulations timings:
t_ch = stimTable.StimTrighCh(i);
T=SA.getDigitalTriggers;
stims = T.tTrig{t_ch};
stimStartT = stims(1);
stimEndT = stims(end);
firstTrig=stims(1:8:end-2);
endStim=stims(8:8:end)+400;
trial = reshape(stims,[8,length(stims)/8])';

pre=20000;
post=100000;
ch = 17;

%% plot plot trace + raster:

 % for j=22:numel(firstTrig)
    j = 22;
    pTmp=find(DB.t_ms>(firstTrig(j)-pre) & DB.t_ms<(firstTrig(j)+post));
    %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
    dbTmp=DB.bufferedDelta2BetaRatio(pTmp);
    [lfp,lfp_t] = SA.currentDataObj.getData(ch,firstTrig(j)-pre,pre+post);
    
    % add a raster plot according to stimulation
    % Load spike times and cluster IDs (adjust file names as needed)
    curPhyFoler = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort4'];
    spikes=load([curPhyFoler filesep 'spike_data.mat']);

    % spikeTimes = readNPY([curPhyFoler filesep 'spike_times.npy'] ); % spike times in samples.
    spikeTimes_ms = spikes.spike_times/ (SA.currentDataObj.samplingFrequency(ch)/1000); %spike times in ms.
    spikeClusters = spikes.spike_clusters; % cluster ID per spike
    cluster_info = readtable([curPhyFoler filesep 'cluster_info.tsv'], 'FileType', 'text', 'Delimiter', '\t');

    % Parameters for the plot
    allClusters  = cluster_info.cluster_id(find(contains(cluster_info.group,'good')|contains(cluster_info.group,'mua')));
    goodClusters = cluster_info.cluster_id(find(contains(cluster_info.group,'good')));
    % uniqueClusters = unique(spikeClusters); % all neuron/unit ID
    nClusters = length(allClusters);
   
%% plot trace + raster:
    figure;
    h1 = subplot(2,1,1);
    % set(f, 'Position', [100, 100, 1200, 400]);
    % yyaxis left
    plot(lfp_t/1000,squeeze(lfp),'k'); hold on;
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims/1000,'r')
    % yyaxis right
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
    xline(curstims,'r')
    box on
    % Convert x-axis labels from ms to s by setting the x-axis ticks and labels
    xticks = get(gca, 'XTick');          % Get current x-axis tick values in ms
    set(gca, 'XTick', xticks);           % Set the same ticks
    set(gca, 'XTickLabel', xticks / 1000); % Display tick labels in seconds
    xlabel('Time (s)');
    % xlabel('Time (ms)');
    ylabel('Neuron/Unit');
    ylim([1,nClusters+1])
    title('Raster Plot');
    hold off;

    % saveas(gcf, [SA.currentPlotFolder filesep sprintf('DBRaterLFPT%i.pdf',j)])
    % waitforbuttonpress;

% end


% saveas (gcf, [analysisFolder filesep sprintf('singleTrialRasterPV161N13t%i.pdf',j)]);


%% spike rates for each recording: 

meanWin = 1000;
win = pre+post;
OL = 100;
load([curPhyFoler filesep 'spikeRateAll.mat'],"spikeRateAll","spikeRateT")

%% spikerate for the rec - in you run the previous block don't run this
spikeRateT = meanWin/2:OL:win-meanWin/2; %time in ms
spikeRateAll = zeros(numel(allClusters),length(spikeRateT),numel(firstTrig));
for j = 1:numel(firstTrig)
    [curSpikeRate] = getSpikeRate(spikes,allClusters,firstTrig(j)-pre,win,meanWin,OL);
    spikeRateAll(:,:,j)= curSpikeRate;

end
% save([curPhyFoler filesep 'spikeRateGoodunits.mat'],"spikeRateGM","spikeRateT",'-mat')
% load([curPhyFoler filesep 'spikeRateTrigers.mat'],"spikeRateM")
save([curPhyFoler filesep 'spikeRateAll.mat'],"spikeRateAll","spikeRateT",'-mat')

%% plot unit avarage:  1 night
unitMean = mean(spikeRateAll,3);

f= figure;
bestU = 48;
bestUind = find(allClusters==bestU);
h1 = subplot(3,1,1);
plot(spikeRateT/OL,unitMean(bestUind,:),'Color','k', 'LineWidth',1)
hold on
stimDiff = mean(mean(diff(trial,[],2)));
xPositions = (pre + (0:7)*(stimDiff))/OL;  % Example positions for lines
xline(xPositions,'r','LineWidth',1)
ylabel('Spike/s')
xlabel('Time[s]')
title('one unit mean across all trials')
hold off

h2=subplot(3,1,2);
[~, Gindices] = ismember(goodClusters, allClusters);
meanGood = mean(unitMean(Gindices,:),1);

% plot(spikeRateT/OL,unitMean(Gindices,:),'Color',[0.7 0.7 0.7], 'LineWidth',1.5)
hold on;
title ('Mean all good units')
plot(spikeRateT/OL,meanGood,'Color','k', 'LineWidth',1)
xline(xPositions,'r','LineWidth',1)
xlabel('Time[s]'); 
ylabel('Spike/s')

h3=subplot(3,1,3);
% plot(spikeRateT/OL,unitMean,'Color',[0.7 0.7 0.7], 'LineWidth',1.5)
hold on;
title ('Mean all units')
plot(spikeRateT/OL,mean(unitMean,1),'Color','k', 'LineWidth',1)
xline(xPositions,'r','LineWidth',1)
xlabel('Time[s]'); 
ylabel('Spike/s')
%savefigure
set(f,'PaperPositionMode','auto');
fileName=[SA.currentPlotFolder filesep 'spikeRates'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% plot spike rate all red nights
% select the relevant recordings:
% the ones that have spikes and are red nights that have any reaction
%spike rate parameters:
spikeRecs = [5,17,19,22];
pre=20000;
post=100000;
meanWin = 1000;
win = pre+post;
OL = 100;
spikeRateT = meanWin/2:OL:win-meanWin/2; %time in ms


%%
for k = 1:length(spikeRecs)
    i= spikeRecs(k);
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    
    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch};
    firstTrig=stims(1:8:end-2);
    endStim=stims(8:8:end)+400;
    % trial = reshape(stims,[8,length(stims)/8])';

    % get the spikes Data
    % Load spike times and cluster IDs (adjust file names as needed)
    curPhyFoler = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort4'];
   
    %RUN THE CODE IN PYHTON TO CHANGE THE SPIKES FILES TO MAT
    spikes=load([curPhyFoler filesep 'spike_data.mat']);
    
    spikeTimes_ms = spikes.spike_times/ (SA.currentDataObj.samplingFrequency(1)/1000); %spike times in ms.
    spikeClusters = spikes.spike_clusters; % cluster ID per spike
    cluster_info = readtable([curPhyFoler filesep 'cluster_info.tsv'], 'FileType', 'text', 'Delimiter', '\t');

    % Parameters for the plot
    allClusters  = cluster_info.cluster_id(find(contains(cluster_info.group,'good')|contains(cluster_info.group,'mua')));
    goodClusters = cluster_info.cluster_id(find(contains(cluster_info.group,'good')));
    % uniqueClusters = unique(spikeClusters); % all neuron/unit ID
    nClusters = length(allClusters);

    %get the spike rate:
    spikeRateFile = [curPhyFoler filesep 'spikeRateAll.mat'];
    if exist(spikeRateFile , 'file') == 2
        disp('Spike Rate File exists.');
    continue % go to the next ittereation 
    else
        
        spikeRateT = meanWin/2:OL:win-meanWin/2; %time in ms
        spikeRateAll = zeros(numel(allClusters),length(spikeRateT),numel(firstTrig));
        for j = 1:numel(firstTrig)
            [curSpikeRate] = getSpikeRate(spikes,allClusters,firstTrig(j)-pre,win,meanWin,OL);
            spikeRateAll(:,:,j)= curSpikeRate;

        end
        % save([curPhyFoler filesep 'spikeRateGoodunits.mat'],"spikeRateGM","spikeRateT",'-mat')
        % load([curPhyFoler filesep 'spikeRateTrigers.mat'],"spikeRateM")
        save(spikeRateFile, "spikeRateAll","spikeRateT",'-mat')
    end
end

%% plot all neurons from all nights traces & baseline changes
% ITI vector for all units and ISI vector
% timings: ITI- interTrial: first 200 s, and last 200 s vs ISI: seconds 2:4 from each stimulation
ITIt = find(spikeRateT<20000 | spikeRateT>110000);

ITIsamples = [];
ISIsamples = [];

AllNightsUnits = [];
goodUnits = [];
% loop on all relevant nights
for k = 1:length(spikeRecs)
    i= spikeRecs(k);
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    % get the stimulation times:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch};
    trial = reshape(stims,[8,length(stims)/8])';
    stimDiff = mean(mean(diff(trial,[],2)));
    xPositions = (pre + (0:7)*(stimDiff))/OL;

    %    from each night : get mean for each unit
    curPhyFoler = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'Kilosort4'];
    spikeRateFile = [curPhyFoler filesep 'spikeRateAll.mat'];
    load(spikeRateFile,"spikeRateAll");

    unitM = mean(spikeRateAll,3);
    AllNightsUnits = [AllNightsUnits;unitM];%add to one matrix for all units, all nights

    % get an array for logical for all the good clusters. 
    cluster_info = readtable([curPhyFoler filesep 'cluster_info.tsv'], ...
        'FileType', 'text', 'Delimiter', '\t');
    allClusters  = cluster_info.cluster_id(find(contains(cluster_info.group,'good')|contains(cluster_info.group,'mua')));
    goodClusters = cluster_info.cluster_id(find(contains(cluster_info.group,'good')));
    goodUnits = [goodUnits ;ismember(allClusters,goodClusters)];

    % caluclate the mean for the ITI and for ISI
    % normelize to seconds
    curITI = mean(unitM(:,ITIt),2);
    
    ISIt = extract_spike_intervals(spikeRateT, xPositions);
    curISI = mean(unitM(:,ISIt),2);
    
    % save in the relevant vector
    ITIsamples = [ITIsamples;curITI];
    ISIsamples = [ISIsamples;curISI];
end
    
%% Plot Traces and ISI vs ITI:
    % 3 subplots: 1 unit, all units same night, all units

unitMean = mean(spikeRateAll,3);
bestU = 48;
bestUind = find(allClusters==bestU);
stimDiff = mean(mean(diff(trial,[],2)));
xPositions = (pre + (0:7)*(stimDiff))/OL;  % Example positions for lines

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


%% PLOT ISI vs ITI: violin plot with swarm plot. 
N = 4;
groupNames = ["ITI"; "ISI"];
colors = [0.6 0.8 0.6; 1 0.7 0.4];

fall = figure; % all units:
subplot(1,2,1)
IIdata = [ITIsamples, ISIsamples];
violin(IIdata, groupNames, 'facecolor',colors,'edgecolor',[0.6 0.6 0.6],'medc',[]);
hold on
x = [ones(length(ITIsamples),1) 2*ones(length(ITIsamples),1)];
swarmchart(x,IIdata,15,colors,'filled');
ylabel('Spikes/S')

subplot(1,2,2)
plot([1,1.2] ,IIdata, '-o','Color',[0.5 0.5 0.5]);
hold on;
plot([1,1.2] ,mean(IIdata), '-o','Color','k', 'LineWidth',2);
xlim([0.7 1.5])
xticks([1,1.2]);xticklabels(["Inter-Trial Interval", "Inter-Stim Interval"]);
ylabel('Spikes/S')
sgtitle('all units')
n = length(ITIsamples);
[pWilcoxon, ~, statsWilcoxon] = signrank(ITIsamples, ISIsamples);

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.8, 0.7, 0.03, 0.1], 'String', ...
    sprintf('Wilcoxon  p-va: %.4f\n', pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
% savefigures
    set(fall,'PaperPositionMode','auto');
    fileName=[analysisFolder filesep 'ITIISIallunits'];
    print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

fgood = figure; % only good units:
subplot(1,2,1)
goodUnitsL = logical(goodUnits);
IIdataGood = [ITIsamples(goodUnitsL), ISIsamples(goodUnitsL)];
violin(IIdataGood, groupNames, 'facecolor',colors,'edgecolor',[0.6 0.6 0.6],'medc',[]);
hold on
x = [ones(length(IIdataGood),1) 2*ones(length(IIdataGood),1)];
swarmchart(x,IIdataGood,15,colors,'filled');
ylabel('Spikes/S')

subplot(1,2,2)
plot([1,1.2] ,IIdataGood, '-o','Color',[0.5 0.5 0.5]);
hold on;
plot([1,1.2] ,mean(IIdataGood), '-o','Color','k', 'LineWidth',2);
xlim([0.7 1.5])
sgtitle ('good units)')
n = sum(goodUnits); 
xticks([1,1.2]);xticklabels(["Inter-Trial Interval", "Inter-Stim Interval"]);
ylabel('Spikes/S')
[pWilcoxon, ~, statsWilcoxon] = signrank(ITIsamples(goodUnitsL), ISIsamples(goodUnitsL));

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.8, 0.7, 0.03, 0.1], 'String', ...
    sprintf('Wilcoxon p-val: %.4f\n', pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
% avefigures
    set(fgood,'PaperPositionMode','auto');
    fileName=[analysisFolder filesep 'ITIISIgoodunits'];
    print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% Delta 2 Beta general decrease: SWS

%% plot D2B general decrease for each part of the stimulation - SWS parts
% get the data:
dbSWMeans = zeros([height(stimTable),3]);

for i = 1:height(stimTable)
% i =22 ;
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    DB = SA.getDelta2BetaRatio;

    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimStartT = T.tTrig{t_ch}(1);
    stimEndT = T.tTrig{t_ch}(end);

    p = 30*60*1000; %some time diff for the cycle to change.
    ACwin = 2*60*60*1000; % 4 hrs in ms
    cycWin = 2*60*60*1000;
    ACStartTimes = [stimStartT-ACwin,stimStartT+p,stimEndT+p];
    partsTimings = [stimStartT-cycWin,stimStartT+p,stimEndT+p];

    dbSW.Pre = [];
    dbSW.Stim = [];
    dbSW.Post = [];
    parts = fieldnames(dbSW);  % Get a cell array of field names
    
    for j=1:numel(parts)
        %calculate AC and SC for this part:
        part = parts{j};
        SA.getDelta2BetaAC('tStart',ACStartTimes(j),'win',ACwin,'overwrite',1);
        
        SA.getSlowCycles('excludeIrregularCycles',0,'overwrite',1)
        SC = SA.getSlowCycles;
        
        curtimings = [partsTimings(j) partsTimings(j)+cycWin];
        pCyc = find(SC.TcycleOnset>= curtimings(1)& SC.TcycleOnset<=curtimings(2));
        curCyclesOns = SC.TcycleOnset(pCyc);
        curCyclesMids = SC.TcycleMid(pCyc);

        for k=1:numel(curCyclesOns)
        %get the sws timings and the DB for them:
            pTmp = find(DB.t_ms>curCyclesOns(k) & DB.t_ms<curCyclesMids(k));
            dbSW.(part)(k) = mean(DB.bufferedDelta2BetaRatio(pTmp));
        end
    end
    
    dbSWMeans(i,:) = [mean(dbSW.Pre,'omitnan') mean(dbSW.Stim,'omitnan') mean(dbSW.Post,'omitnan')];
    %save in stimTable
    stimTable.dbSW(i) = dbSW;
end

stimTable.dbSWMeans = dbSWMeans;
%%
save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');

%% plot one night: D/B decrease
 
% for i = 1:height(stimTable)
i =22 ;
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    
    fdbDec = figure;
    groupNames = ["Pre", "Stim","Post"];
    colors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];
    dbStimData = {stimTable.dbSW(i).Pre stimTable.dbSW(i).Stim stimTable.dbSW(i).Post};
    violin(dbStimData, groupNames, 'facecolor',colors,'edgecolor',[0.6 0.6 0.6],'medc',[]);
    hold on
    % jitterAmount = 0.4;
    % for j = 1:numel(dbStimData)
    %     % Generate x-coordinates with jitter for each group
    %     x = j + (rand(size(dbStimData{j})) - 0.5) * jitterAmount;
    % 
    %     % Plot data points for this group
    %     scatter(x, dbStimData{j}, 20, colors(j,:), 'filled', 'MarkerFaceAlpha', 0.5); % Adjust size and transparency
    % end
    yswarm = [stimTable.dbSW(i).Pre stimTable.dbSW(i).Stim stimTable.dbSW(i).Post];
    xswarm=[1*ones(1,length(stimTable.dbSW(i).Pre)), 2*ones(1,length(stimTable.dbSW(i).Stim)), 3*ones(1,length(stimTable.dbSW(i).Post))];
    colorswarm = [repmat([0.2, 0.6, 0.8;], length(stimTable.dbSW(i).Pre), 1);  % Red for Array 1
          repmat([0.9, 0.4, 0.3], length(stimTable.dbSW(i).Stim), 1);  % Green for Array 2
          repmat([0.5, 0.8, 0.5], length(stimTable.dbSW(i).Post), 1)]; % Blue for Array 3
    swarmchart(xswarm,yswarm,15,colorswarm,"filled","o")

    ylabel('D/B power during SWS bouts')
    
    set(gca, 'XTick', 1:numel(dbStimData), 'XTickLabel', groupNames);
    title('D/B power during SWS bouts - one night')
    
    %savefigure
    set(fdbDec,'PaperPositionMode','auto');
    fileName=[SA.currentPlotFolder filesep 'DBSWSpreStimPost'];
    print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% plot all nights: - only red

type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength); %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));


% statistics:

% Assuming data in columns where each row is a subject and each column is a timepoint
bdSWDataStat = stimTable.dbSWMeans(curTrials,:);
cleanDBSW = bdSWDataStat(~any(isnan(bdSWDataStat), 2), :);
n=height(cleanDBSW);
[p, tbl, stats] = friedman(cleanDBSW, 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
before = stimTable.dbSWMeans(curTrials,1);
during = stimTable.dbSWMeans(curTrials,2);
after = stimTable.dbSWMeans(curTrials,3);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_before_during, ~, stats_before_during] = signrank(before, during);
[p_during_after, ~, stats_during_after] = signrank(during, after);
[p_after_before, ~, stats_after_before] = signrank(after, before);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('After vs Before: p-value = %.4f (Significant if < %.4f)\n', p_after_before, alpha);

%plot
figure;
plot(stimTable.dbSWMeans(curTrials,:)','Color',[0.5 0.5 0.5],'Marker','.','MarkerSize',10)
hold on; 
plot(mean(stimTable.dbSWMeans(curTrials,:),1,'omitnan'),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
xlim([0.5, 3.5])
grid on;
ylim([0 450])
xticks(1:3)
xticklabels(groupNames)
ylabel('D/B means during SWS')

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
    sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.15, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f (Significant if < %.4f)', p_before_during, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_after_before), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% savefigure
set(gcf,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'DBSWS'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
 %% plot D/B decrease - all nights

animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","WhiteEx"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);
markers = {'o', 's', 'd', '^', 'v','x'};
x=1:3;
%plot Period Times:
f=figure;
set(f, 'Position', [100, 100, 800, 400]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType)&~contains(stimTable.Remarks,"Ex");
    n = sum(curTrials);
    curCol = plotColors{type};
    curMean = mean(stimTable.ACcomPer(curTrials,:)/1000,1,'omitnan');
    
    if n>0

        plot(x,curMean,'color',curCol,'LineWidth',4)
        hold on
        for animal = 1:numAnimal
            curAni = animals{animal};
            curcurT = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
            if any(curcurT)
            plot(x, stimTable.ACcomPer(curcurT,:)/1000,'Marker',markers{animal},'Color',curCol ...
                ,'MarkerFaceColor',curCol)
            end
        end
        hold off

        annotation('textbox', [.05 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('Time[s]')
    ylim([40 250])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('Perios Times According to stim wavelangth and animal')

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'ACperiodstimColor.pdf']);


%% AC changes 

%% plot sliding AC sith stimulations:
%set the recording:
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
AC = SA.getDelta2BetaAC;
SA.plotDelta2BetaSlidingAC ('stim',1,'stimCh',stimTable.StimTrighCh(i));

%% AC - plot specific before during and after AC - for a specific rec/
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
ACpre = stimTable.ACpre{i};
ACstim = stimTable.ACstim{i};
ACpost = stimTable.ACpost{i};

ACstructs = {ACpre,ACstim,ACpost};
labels = {'preStim', 'Stim', 'postStim'};

for j = 1:3
    struct2vars(ACstructs{j});
    fAC = figure;
    h = axes;
    %plot:
    lineHandles = stem(autocorrTimes/1000,real(xcf),'filled','r-o');
    ylim([-0.4 1])
    set(lineHandles(1),'MarkerSize',4);
    grid('on');
    xlabel('Period [s]');
    ylabel('Auto corr.');
    hold('on');

    plot(period/1000,real(xcf(pPeriod)),'o','MarkerSize',5,'color','k');
    text(period/1000,0.05+real(xcf(pPeriod)),num2str(period/1000));

    a = axis;
    plot([a(1) a(1); a(2) a(2)],[xcf_bounds([1 1]) xcf_bounds([2 2])],'-b');
    plot([a(1) a(2)],[0 0],'-k');
    hold('off');
    
    % save fig:

    set(fAC,'PaperPositionMode','auto');
    fileName=[SA.currentPlotFolder filesep 'dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win) labels{j} '.pdf'];
    saveas(fAC,fileName)
    % print(fileName,'-pdf',['-r' num2str(SA.figResJPG)]);
    % if printLocalCopy
    %     fileName=[cd filesep obj.recTable.Animal{obj.currentPRec} '_Rec' num2str(obj.currentPRec) '_dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(parDbAutocorr.tStart) '_w' num2str(parDbAutocorr.win)];
    %     print(fileName,'-pdf',['-r' num2str(obj.figResJPG)]);
    % end

end

%% AC - get the Data 

ACprePer = zeros(height(stimTable),1);
ACpreP2V = zeros(height(stimTable),1);
ACstimPer = zeros(height(stimTable),1);
ACstimP2V = zeros(height(stimTable),1);
ACpostPer = zeros(height(stimTable),1);
ACpostP2V = zeros(height(stimTable),1);

for i = 1:height(stimTable)
   curACpre = stimTable.ACpre{i};
   curACstim = stimTable.ACstim{i};
   curACpost = stimTable.ACpost{i};
   ACprePer(i) = curACpre.period;
   ACpreP2V(i) = curACpre.peak2VallyDiff;
   ACstimPer(i) = curACstim.period;
   ACstimP2V(i) = curACstim.peak2VallyDiff;
   ACpostPer(i) = curACpost.period;
   ACpostP2V(i) = curACpost.peak2VallyDiff;
   
end


%% plot AC - only Red nights:
type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength); %& contains(stimTable.Animal,curAni);
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
% Example data for three groups
before = stimTable.ACcomPer(curTrials,1);
during = stimTable.ACcomPer(curTrials,2);
after = stimTable.ACcomPer(curTrials,3);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_before_during, ~, stats_before_during] = signrank(before, during);
[p_during_after, ~, stats_during_after] = signrank(during, after);
[p_after_before, ~, stats_after_before] = signrank(after, before);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('After vs Before: p-value = %.4f (Significant if < %.4f)\n', p_after_before, alpha);

%plot
figure;
x = 1:3;
curMean = mean(stimTable.ACcomPer(curTrials,:)/1000,1,'omitnan');
plot(x, stimTable.ACcomPer(curTrials,:)/1000,'Color',[0.5 0.5 0.5],'Marker','.','MarkerSize',10)
hold on
plot(x, curMean,'Color','k','LineWidth',2)
xlim([0.75 3.25])
xticks(x)  % Set ticks after xlim to avoid automatic adjustment
xticklabels({'Pre', 'During', 'Post'})
grid on
ylim([50 300])
ylabel('Period Time[s]')
title ('Perios Times changes - all red nights')
annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
    sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.15, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f (Significant if < %.4f)', p_before_during, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('p-value = %.4f', p_after_before), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


% savefigure
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'ACperiodReds.pdf']);
%% plot AC - According to color and animal

animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);
markers = {'o', 's', 'd', '^', 'v','x'};
x=1:3;
%plot Period Times:
f=figure;
set(f, 'Position', [100, 100, 800, 400]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType);
    n = sum(curTrials);
    curCol = plotColors{type};
    curMean = mean(stimTable.ACcomPer(curTrials,:)/1000,1,'omitnan');
    
    if n>0

        plot(x,curMean,'color',curCol,'LineWidth',4)
        hold on
        for animal = 1:numAnimal
            curAni = animals{animal};
            curcurT = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
            if any(curcurT)
            plot(x, stimTable.ACcomPer(curcurT,:)/1000,'Marker',markers{animal},'Color',curCol ...
                ,'MarkerFaceColor',curCol)
            end
        end
        hold off

        annotation('textbox', [.05 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('Time[s]')
    ylim([40 250])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('Perios Times According to stim wavelangth and animal')

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'ACperiodstimColor.pdf']);

%% plot P2V Times:
f=figure;
set(f, 'Position', [100, 200, 800, 300]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType);
    n = sum(curTrials);
    curCol = plotColors{type};
    curMean = mean(stimTable.ACcomP2V(curTrials,:)/1000,1,'omitnan');
    
    if n>0

        plot(x,curMean,'color',curCol,'LineWidth',4)
        hold on
        for animal = 1:numAnimal
            curAni = animals{animal};
            curcurT = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
            if any(curcurT)
            plot(x, stimTable.ACcomP2V(curcurT,:)/1000,'Marker',markers{animal},'Color',curCol ...
                ,'MarkerFaceColor',curCol)
            end
        end
        hold off

        annotation('textbox', [.07 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    ylabel('P2V')
    ylim([0 0.0015])
    xticklabels({'Pre','During','Post'})
    xticks(1:3); xlim([0.5 3.5])
end

sgtitle ('P2V According to stim wavelangth and animal')

% savefigure
set(gcf,'PaperPosition',[.25 3 6 2])
saveas (gcf, [analysisFolder filesep 'ACP2VstimColor.pdf']);
% clearvars -except stimTable SA analysisFolder
%% StimSham analysis:

%% plot stim according to animal and color
% assuming the table is in the workspace
% load([analysisFolder filesep 'stimTable.mat'])

animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);

times = stimTable.times{1};
pre=50000;
% STIM DUTRATION NEEDS A THINK!!!!!!!
mStimDur = mean(stimTable.stimDuration(3:end));
post=100000;

f=figure;
set(f, 'Position', [100, 100, 1200, 800]);

ylims = [0 300];
for animal = 1:numAnimal
    for type = 1:numType
        h= subplot(numAnimal,numType, (animal-1)*numType+type);
        
        %plot the data
        curAni = animals{animal};
        curType = stimWaveL(type);
        curTrials = contains(stimTable.Remarks,curType) & contains(stimTable.Animal,curAni);
        n = sum(curTrials);
        curCol = plotColors{type};
        curMean = mean(cell2mat(stimTable.StimAvg(curTrials)),1,'omitnan');
        if n>0
            plot(times,curMean,'color',curCol,'LineWidth',4)
            hold on
            for i = 1:height(stimTable)
                if curTrials(i) && stimTable.Mani(i) == 1
                    plot(times,stimTable.StimAvg{i},'Color',curCol)
                end
                if curTrials(i) && stimTable.Mani(i) == 2
                    plot(times,stimTable.StimAvg{i},'Color',curCol,'LineStyle','--')
                end

            end

            xline(pre/1000, 'LineWidth',1,'Color','black')
            xline((pre+mStimDur)/1000,'LineWidth',1,'Color','black')
            hold off
            
            annotation('textbox', [.095 + 0.195*type, 0.98 - 0.142*animal, 0.03, 0.1], 'String', ...
                sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
                'right', 'VerticalAlignment', 'middle');

        end
        if n==0
            plot(0,0)
        end

        % add titles. labels...
        ylim(ylims);xlabel('Time[S]'), ylabel('D2B power')
        if animal == 1
            title(stimType(type))
        end
%         legend()
      
        pos = get(h, 'Position');
        pos(1) = pos(1) + 0.03; % Shift subplot right to make space for row titles
        set(h, 'Position', pos);
        
        % Add row titles using text function
        if type == 1
            ypos = (ylims(1) + ylims(2)) / 2; % Center vertically

            annotation('textbox', [0.1, 0.98 - 0.142*animal, 0.03, 0.1], 'String', animals{animal}, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
                'FontWeight', 'bold');
        
        end
    end 
end


% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'animal_type.pdf']);
clearvars -except stimTable SA analysisFolder

%% plot the statistics:
% i want a bar plot (with points) of the max values of the D/B during
% stim/Sham
% also, plot the max slop, i think that will show better results.

% assuming the table is in the workspace
% load([analysisFolder filesep 'stimTable.mat'])

% calculate max value and max slop
% max:
stimTable.maxStim = cellfun(@max,stimTable.StimAvg);
stimTable.maxSham = cellfun(@max,stimTable.StimAvgSham);

% stimtime minus pre stim - change in D/B:
% take the last 30 sec of stim and substract the 30 sec before the start of
% stim, ans avrage that. 
stimTable.dbDiffStim = cell(height(stimTable),1);
stimTable.dbDiffSham = cell(height(stimTable),1);
firstStimInd = 50000/1000;
win = 30;
for i=1:height(stimTable)
    CurStimDur = round(stimTable.stimDuration(i)/1000);
    % calc substracted stimulation from baseline
    befStim = stimTable.StimAvg{i}(firstStimInd-win:firstStimInd-1);
    durStim = stimTable.StimAvg{i}(firstStimInd+CurStimDur-win:firstStimInd+CurStimDur-1);
    stimTable.dbDiffStim(i) = {durStim-befStim};
    % calc substracted stimulation from baseline
    befSham = stimTable.StimAvgSham{i}(firstStimInd-win:firstStimInd-1);
    durSham = stimTable.StimAvgSham{i}(firstStimInd+CurStimDur-win:firstStimInd+CurStimDur-1);
    stimTable.dbDiffSham(i) = {durSham-befSham};
end
stimTable.dbDiffStimM = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffStim);
stimTable.dbDiffShamM = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffSham);

save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');
 
%% plot D/B diff stimSham bar plot - only Red nights
% animals = unique(stimTable.Animal);
% numAnimal = length(animals);
type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ...
    ~contains(stimTable.Remarks,'Ex') & ...
    (~any(isnan(stimTable.dbDiffShamM), 2) & ~any(isnan(stimTable.dbDiffStimM), 2));% &...
    % ~contains(stimTable.Animal,'157');
    
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

% statistical tests for that figure - diff StimSham only red nights:
% trying different methods: 
groupSham = stimTable.dbDiffShamM(curTrials);
groupStim = stimTable.dbDiffStimM(curTrials);



% 1. Wilcoxon Signed-Rank Test
[pWilcoxon, ~, statsWilcoxon] = signrank(groupSham, groupStim);
fprintf('Wilcoxon Signed-Rank Test p-value: %.4f\n', pWilcoxon);
disp(statsWilcoxon)

% 2. Paired t-test (check normality first)
% Normality test (Shapiro-Wilk) - used
[hNorm1, pNorm1] = adtest(groupSham);
[hNorm2, pNorm2] = adtest(groupStim);

if pNorm1 > 0.05 && pNorm2 > 0.05
    [hTtest, pTtest] = ttest(groupSham, groupStim);
    fprintf('Paired t-test p-value: %.4f\n', pTtest);
else
    fprintf('Data is not normally distributed; Paired t-test may not be appropriate.\n');
end


%plot:
figure;
title('Change in mean D/B norm across trials - All nights')
x=[1,2];
plot(x,[stimTable.dbDiffShamM(curTrials), stimTable.dbDiffStimM(curTrials)] ...
    ,'Color',[0.5, 0.5, 0.5],'Marker','.','MarkerSize',10)
hold on
curMeanStim = mean(stimTable.dbDiffStimM(curTrials),1,'omitnan');
curMeanSham = mean(stimTable.dbDiffShamM(curTrials),1,'omitnan');
plot(x,[curMeanSham,curMeanStim],'color','k','LineWidth',2,'Marker','.','MarkerSize',10)
hold off

grid on
xlim([0.5 2.5]);
xticks([1, 2]); % Position of the x-ticks
xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
annotation('textbox', [0.85, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.15, 0.8, 0.3, 0.2], 'String', ...
    sprintf('Wilcoxon test p =%.3f',pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.15, 0.75, 0.3, 0.2], 'String', ...
    sprintf('T-test p =%.3f',pTtest), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle')

set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'DBdiffStimShamRedNights.pdf']);


clearvars -except stimTable SA analysisFolder

%% plot the bar plot - D/B change - stim sham -  all colors
animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {"blue","green","red", [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);


f=figure;
set(f, 'Position', [100, 100, 1200, 400]);
sgtitle('Max average D/B')
for type = 1:numType
    h = subplot(1,numType,type);
    %plot the data
    %curAni = animals{animal};
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType); %& contains(stimTable.Animal,curAni);
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    curCol = plotColors{type};
    curMeanStim = mean(stimTable.maxStim(curTrials),1,'omitnan');
    curMeanSham = mean(stimTable.maxSham(curTrials),1,'omitnan');
    if n>0
        x=[1,2];
        plot(x,[curMeanSham,curMeanStim],'-o','color',curCol,'LineWidth',3)
        hold on
        plot(x,[stimTable.maxSham(curTrials), stimTable.maxStim(curTrials)] ...
            ,'-o','Color',curCol)
          hold off
    end
    xticks([1, 2]); % Position of the x-ticks
    xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
    xlim([0.5, 2.5]);
    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   
    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('D2B power')
    title(stimType(type))
    ylim([0 450])
end

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'maxStimShamAll1.pdf']);
%% mean normelized change in D/B.
f=figure;
set(f, 'Position', [100, 100, 1200, 400]);
sgtitle('mean normelized D/B')
for type = 1:numType
    h = subplot(1,numType,type);
    %plot the data
    %curAni = animals{animal};
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType); %& contains(stimTable.Animal,curAni);
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    curCol = plotColors{type};
    curMeanNdbStim = mean(stimTable.dbDiffStimM(curTrials),1,'omitnan');
    curMeanNdbSham = mean(stimTable.dbDiffShamM(curTrials),1,'omitnan');
    if n>0
        x=[1,2];
        plot(x,[curMeanNdbSham,curMeanNdbStim],'-o','color',curCol,'LineWidth',3)
        hold on
        plot(x,[stimTable.dbDiffShamM(curTrials), stimTable.dbDiffStimM(curTrials)] ...
            ,'-o','Color',curCol)
          hold off
    end
    xticks([1, 2]); % Position of the x-ticks
    xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
    xlim([0.5, 2.5]);
    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   ylim([-40 160])
    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('D2B power')
    title(stimType(type))
 
end

% savefigure
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'meanNormBDStimSham.pdf']);

% clearvars -except stimTable SA analysisFolder


%% Lizard Movement Data analysis:

%% Get LM DATA

LMdata = getLMData(SA, stimTable,analysisFolder);

%% check the correlation of mevement with th D/B
i = 19;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
LM_DBt= stimTable.LM_DBt{i};
DB = SA.getDelta2BetaRatio;

%%
figure;
yyaxis left
plot(DB.bufferedDelta2BetaRatio)
hold on;
yyaxis right
plot(LM_DBt,'k')
%% plot the full movement data for a night

% for one night:
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
% SA.setCurrentRecording(recName);
% DB = SA.getDelta2BetaRatio;

curLMwake = LMdata.LMwake(i);
curLMpre = LMdata.LMpre(i);
curLMstim = LMdata.LMstimbin{i};
curLMpost = LMdata.LMpost(i);

figure;
hold on;
xWak = 1;
xBef = 2;
xdur = linspace(3,5,length(curLMstim));
xaft = 6;

plot(xWak,curLMwake,'.','Color','black', 'MarkerSize',20);
plot(xBef,curLMpre,'.','Color','black','MarkerSize',20);
plot(xdur,curLMstim,'-o','Color','black','MarkerFaceColor','black');
plot(xaft,curLMpost,'.','Color','black','MarkerSize',20);
xlim([0.5,xaft+0.5])
xticklabels({'Wake','sleep Before','during stimulation','sleep After'})
xticks([xWak,xBef,xdur(round(length(curLMstim)/2)),xaft])
title('Mean movement during stimulation, PV161,Night18')
ylabel('Mov/s')
hold on;
xline(xdur(1),'Color','r','LineWidth',2);
xline(xdur(4),'Color','r','LineWidth',2)
hold off;
% save fig
set(gcf,'PaperPosition',[1,5,3.5,4])
saveas (gcf, [analysisFolder filesep 'lizMovWholeNightPV161N18.pdf']);



%% PLOT - Movement during stimulation - All Nights 
% change LMbin to mat:
LMstimbinM = cell2mat(LMdata.LMstimbin);
LMstimbinMean = mean(LMstimbinM,1);
n = height(LMstimbinM);

figure;
hold on;
xWak = 1;
xBef = 2;
xdur = linspace(3,5,length(LMstimbinMean));
xaft = 6;

plot(xWak,LMdata.LMwake,'.','Color','black', 'MarkerSize',20);
plot(xBef,LMdata.LMpre,'.','Color','black','MarkerSize',20);
plot(xdur,LMstimbinM','-','Color',[0.75, 0.75, 0.75]);
plot(xdur,LMstimbinMean,'-','Color','black','LineWidth',2);
plot(xaft,LMdata.LMpost,'.','Color','black','MarkerSize',20);

xlim([0.5,xaft+0.5])
% ylim([0 100])
Groups = {'Wake','Sleep Before','during stimulation','sleep After'};
xticklabels(Groups)
xticks([xWak,xBef,xdur(round(length(LMstimbinMean)/2)),xaft])
ylabel('Mov/s')

title('Mean movement during stimulation, All Nights, norm')
annotation('textbox', [0.5, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i',n), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

xline(xdur(1),'Color','r','LineWidth',2);
xline(xdur(4),'Color','r','LineWidth',2)
    
% save fig
set(gcf,'PaperPosition',[.25 3 8 6])
saveas (gcf, [analysisFolder filesep 'lizMovAllNights_norm.pdf']);

% clearvars -except stimTable SA analysisFolder
%% plot with mov with violin :
% get data from table - no NAN!
LMpre = LMdata.LMpre(~any(isnan(LMdata.LMpre), 2), :);
LMwake = LMdata.LMwake(~any(isnan( LMdata.LMwake), 2), :);
LMpost = LMdata.LMpost(~any(isnan(LMdata.LMpost), 2), :);
LMstimbinM = cell2mat(LMdata.LMstimbin); % takes out the nan val
LMstimbintrialM = mean(LMstimbinM,2); % mean for eact night

LMvioData = [LMwake, LMpre, LMstimbintrialM, LMpost];
n = height(LMvioData);

% check the statistics:
% Assuming data in columns where each row is a subject and each column is a timepoint

[p, tbl, stats] = friedman(LMvioData, 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Bonferroni-corrected alpha level
alpha = 0.05 / 6;
% Pairwise Wilcoxon signed-rank tests
[p_wake_pre, ~, stats_wake_pre] = signrank(LMwake, LMpre);
[p_wake_during, ~, stats_wake_during] = signrank(LMwake, LMstimbintrialM);
[p_wake_after, ~, stats_wake_after] = signrank(LMwake, LMpost);
[p_pre_during, ~, stats_pre_during] = signrank(LMstimbintrialM,LMpre);
[p_during_after, ~, stats_during_after] = signrank(LMstimbintrialM, LMpost);
[p_pre_after, ~, stats_pre_after] = signrank(LMpre, LMpost);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Wake vs Pre: p-value = %.5f (Significant if < %.4f)\n', p_wake_pre, alpha);
fprintf('Wake vs During: p-value = %.5f (Significant if < %.4f)\n', p_wake_during, alpha);
fprintf('Wake vs Post: p-value = %.5f (Significant if < %.4f)\n', p_wake_after, alpha);
fprintf('During vs Pre: p-value = %.5f (Significant if < %.4f)\n', p_pre_during, alpha);
fprintf('During vs Post: p-value = %.5f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('Pre vs Post: p-value = %.5f (Significant if < %.4f)\n', p_pre_after, alpha);

fileName = [analysisFolder filesep 'postHocPvalLMviolinAllNights.mat'];
save(fileName, 'alpha','p_wake_pre','p_wake_during','p_wake_after','p_pre_during', ...
    'p_during_after','p_pre_after')


fLM = figure;
set(fLM, 'PaperPositionMode','auto')
Groups = ["Wake","Sleep Before","During Stim","Sleep After"];
colors = [0.5 0.5 0.5;0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];
violin(LMvioData, 'xlabel',Groups, 'facecolor',colors,'edgecolor',[0.6 0.6 0.6],'medc',[]);

% savefigure
set(fLM,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'LMviolinallnights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

clearvars -except stimTable SA analysisFolder LMdata
%% plot polar histogram 
% all nights:
mPhasePre.movs = zeros(height(stimTable),1);
mPhasePre.DBs = zeros(height(stimTable),1);
mPhaseStim.movs = zeros(height(stimTable),1);
mPhaseStim.DBs = zeros(height(stimTable),1);
mPhasePost.movs = zeros(height(stimTable),1);
mPhasePost.DBs = zeros(height(stimTable),1);

%%
% for i = 34:height(stimTable)
i = 22;
    if isempty(stimTable.LM_DBt{i})
        disp('run Lizard movement on this recording. Moving to next rec.');
        mPhasePre.movs(i) = NaN;
        mPhasePre.DBs(i) = NaN;
        mPhaseStim.movs(i) = NaN;
        mPhaseStim.DBs(i) = NaN;
        % continue; % Skip to the next iteration;
    end

    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimStartT = T.tTrig{t_ch}(1);
    stimEndT = T.tTrig{t_ch}(end);

    ACwin = 2*60*60*1000;
    p = 30*60*1000;
    % get AC for the pre:
    SA.getDelta2BetaAC('tStart',stimStartT-ACwin ,'win',ACwin , 'overwrite', 1);
    SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);
    % SA.plotSlowCycles
    %for the beggining of the stimulations
    hOutPre = SA.plotLizardMovementDB('stim',1 ,'part',1,'tStartStim', ...
        stimStartT,'tEndStim',stimEndT);
    mPhasePre.movs(i) = hOutPre.mPhaseMov;
    mPhasePre.DBs(i) = hOutPre.mPhaseDB;
    % get AC for the stimulation:
    SA.getDelta2BetaAC('tStart',stimStartT,'win',ACwin,'overwrite', 1);
    SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);

    %for the beggining of the stimulations
    hOutStim = SA.plotLizardMovementDB('stim',1 ,'part',2,'tStartStim', ...
        stimStartT,'tEndStim',stimEndT);
    mPhaseStim.movs(i) = hOutStim.mPhaseMov;
    mPhaseStim.DBs(i) = hOutStim.mPhaseDB;

        % get AC for the Post:
    SA.getDelta2BetaAC('tStart',stimEndT+p,'win',ACwin,'overwrite', 1);
    SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);

    %for the beggining of the stimulations
    hOutPost = SA.plotLizardMovementDB('stim',1 ,'part',3,'tStartStim', ...
        stimStartT,'tEndStim',stimEndT);
    mPhasePost.movs(i) = hOutPost.mPhaseMov;
    mPhasePost.DBs(i) = hOutPost.mPhaseDB;

    % close all
    % end

%%
fileName = [analysisFolder filesep 'polarhistoAllNights.mat'];
save(fileName,"mPhasePre" ,"mPhaseStim","mPhasePost")


%% plot polser histogram - All nights
uniqueAnimals = unique(stimTable.Animal);
relativePhasePre = mPhasePre.movs -mPhasePre.DBs;
relativePhaseStim = mPhaseStim.movs -mPhaseStim.DBs;
relativePhasePost = mPhasePost.movs -mPhasePost.DBs;

% pVal = cellfun(@(x) ~isempty(x),stimTable.LM_DBt);
pVal = mPhasePost.movs ~= 0;
% wavelength = '635';
% pVal = mPhasePost.movs ~= 0 & contains(stimTable.Remarks,wavelength);

animalColor = [
    0.2, 0.6, 0.8;   % soft blue
    0.8, 0.4, 0.2;   % warm orange
    0.4, 0.7, 0.4;   % gentle green
    0.6, 0.3, 0.6;   % muted purple
    0.7, 0.5, 0.3;   % earthy brown
    0.3, 0.6, 0.5    % teal
];

fMOVdb=figure;
h1=subplot(1,3,1,polaraxes);hold on;
title('Pre Stimulations')
Rlim=0.5;

hP={};
for i=1:numel(uniqueAnimals)
    p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
    hP{i}=polarplot([relativePhasePre(p)';relativePhasePre(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalColor(i,:),'LineWidth',1);
end
hold on;
hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);

hRose=polarhistogram(h1,relativePhasePre(pVal),12,'Normalization','probability');
hRose.FaceColor=[0.7 0.7 0.7];
hRose.FaceAlpha=0.5;

text(0.2, Rlim/2, '\delta/\beta');
h1.ThetaTick=[0:90:330];
h1.RTick=[0.1:0.1:0.4];
%h2.ThetaTickLabels([2 3 5 6 8 9 11 12])=cell(size([2 3 5 6 8 9 11 12]));

l1=legend([hP{2}(1),hP3,hRose],{'singleNight','\delta/\beta','Prob.'},'box','off');
l1.Position=[0.7386    0.8238    0.2125    0.1190];

%figure 2 : during stim:
h2=subplot(1,3,2,polaraxes);hold on;% Stimulation time
title('During Stimulations')

Rlim=0.5;

hP={};
for i=1:numel(uniqueAnimals)
    p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
    hP{i}=polarplot([relativePhaseStim(p)';relativePhaseStim(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalColor(i,:),'LineWidth',1);
end
hold on;
hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);

hRose=polarhistogram(h2,relativePhaseStim(pVal),12,'Normalization','probability');
hRose.FaceColor=[0.7 0.7 0.7];
hRose.FaceAlpha=0.5;

text(0.2, Rlim/2, '\delta/\beta');
h2.ThetaTick=[0:90:330];
h2.RTick=[0.1:0.1:0.4];

%figure 3 : post stim:
h3=subplot(1,3,3,polaraxes);hold on;% Stimulation time
title('Post Stimulations')

Rlim=0.5;

hP={};
for i=1:numel(uniqueAnimals)
    p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
    hP{i}=polarplot([relativePhasePost(p)';relativePhasePost(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalColor(i,:),'LineWidth',1);
end
hold on;
hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);

hRose=polarhistogram(h3,relativePhasePost(pVal),12,'Normalization','probability');
hRose.FaceColor=[0.7 0.7 0.7];
hRose.FaceAlpha=0.5;

text(0.2, Rlim/2, '\delta/\beta');
h3.ThetaTick=[0:90:330];
h3.RTick=[0.1:0.1:0.4];


% savefigure

set(gcf, 'PaperUnits', 'inches');         % Set paper units to inches
set(gcf, 'PaperSize', [6,4]);            % Set the paper size (width x height in inches)
set(gcf, 'PaperPosition', [0, 0, 6,4]);  % Set position on paper to match size exactly

% Print to PDF
print(gcf, 'PolarMovDBAllnights.pdf', '-dpdf', '-r300');  % '-r300' sets resolution to 300 DPI
% clearvars -except stimTable SA analysisFolder LMdata