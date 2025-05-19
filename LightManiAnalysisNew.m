%% light manipulation analysis 08/05/2024
% This is the new version. 

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

%% creating the StimTable:

%% get all the stim sham avg from all recs and put in a new table
% this part goes over all the records in SA.
%  for every recoerd that is tagged (1/2/3..) 
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWakeTest.xlsx');
% SA.setCurrentRecording('Animal=PV162,recNames=Night27');
maniRecs = SA.recTable.Mani>0; % taking all the rows with manipulation
stimTable = SA.recTable(maniRecs,{'Animal','recNames','Remarks','Mani','LizMov','StimTrighCh','Headstage'});  % creating new table
stimTable.stimStartT = zeros(height(stimTable),1);
stimTable.stimEndT = zeros(height(stimTable),1);
stimTable.sleepStartT = zeros(height(stimTable),1);
stimTable.sleepEndT = zeros(height(stimTable),1);
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
    SA.getDelta2BetaAC('tStart',0,'overwrite',1);
    SA.getDigitalTriggers;
    
    
    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimTable.stimStartT(i) = T.tTrig{t_ch}(1);
    stimTable.stimEndT(i) = T.tTrig{t_ch}(end);
    
    ACfull =  SA.getDelta2BetaAC;
    stimTable.sleepStartT(i) = ACfull.tStartSleep;
    stimTable.sleepEndT(i) = ACfull.tEndSleep;
     
    % firstTrig=T.tTrig{t_ch}(1:8:end-2);
    % endStim=T.tTrig{t_ch}(8:8:end)+400;
    % stimDuration=(endStim(1)-firstTrig(1));
 
    
    % get and save the stim avg
    s = getStimSham_og(SA,stimTable.StimTrighCh(i),1);
    disp('got Stim for this rec')
    stimTable.StimAvg(i) = {mean(s.StimDB,1,"omitnan")};
    stimTable.StimAvgSham(i) = {mean(s.StimDBSham,1,"omitnan")};
    stimTable.times(i) = {s.ts};
    stimTable.stimDuration(i) = s.stimDur;
    disp('stimsham in table')



     % calculate the AC and the P2V for each part of the stimulation
    p = 30*60*1000; %some time diff for the cycle to change.
    preWin = stimTable.stimStartT(i) - stimTable.sleepStartT(i); %all the time before stim start
    stimWin = stimTable.stimEndT(i) - (stimTable.stimStartT(i)+p); % stimulation period, not including first 30 min
    postWin = stimTable.sleepEndT(i) - (stimTable.stimEndT(i)+p); %post stimulations sleep, not including first 30 mins.

    SA.getDelta2BetaAC('tStart',stimTable.sleepStartT(i), 'win',preWin,'overwrite',1);
    stimTable.ACpre(i) = {SA.getDelta2BetaAC('tStart',stimTable.sleepStartT(i), 'win',preWin)};

    SA.getDelta2BetaAC('tStart',(stimTable.stimStartT(i)+p), 'win',stimWin,'overwrite',1);
    stimTable.ACstim(i) = {SA.getDelta2BetaAC('tStart',(stimTable.stimStartT(i)+p), 'win',stimWin)};
    
    SA.getDelta2BetaAC('tStart',(stimTable.stimEndT(i)+p),'win',postWin,'overwrite',1);
    stimTable.ACpost(i) = {SA.getDelta2BetaAC('tStart',(stimTable.stimEndT(i)+p),'win',postWin)};
    disp('AC in table')
    
      
end
clear recName
clear s
%% calculate stim sham diff:
stimTable.dbDiffStim = cell(height(stimTable),1);
stimTable.dbDiffSham = cell(height(stimTable),1);
stimTable.dbDiffStimM = zeros(height(stimTable),1);
stimTable.dbDiffShamM = zeros(height(stimTable),1);
firstStimInd = 50000/1000;
win = 30;
for i=1:height(stimTable)
    CurStimDur = round(stimTable.stimDuration(i)/1000);
    % calc substracted stimulation from baseline
    befStim = stimTable.StimAvg{i}(firstStimInd-win:firstStimInd-1);
    durStim = stimTable.StimAvg{i}(firstStimInd+CurStimDur-win:firstStimInd+CurStimDur-1);
    stimTable.dbDiffStim(i) = {durStim-befStim};
    stimTable.dbDiffStimM(i) = mean(durStim,'omitnan')-mean(befStim,'omitnan');
    % calc substracted stimulation from baseline
    befSham = stimTable.StimAvgSham{i}(firstStimInd-win:firstStimInd-1);
    durSham = stimTable.StimAvgSham{i}(firstStimInd+CurStimDur-win:firstStimInd+CurStimDur-1);
    stimTable.dbDiffSham(i) = {durSham-befSham};
    stimTable.dbDiffShamM(i) = mean(durSham,'omitnan')-mean(befSham,'omitnan');

end
stimTable.dbDiffStimM1 = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffStim);
stimTable.dbDiffShamM1 = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffSham);

%% AC - get the Data 

ACcomPer = zeros(height(stimTable),3);
ACcomP2V = zeros(height(stimTable),3);
% ACstimPer = zeros(height(stimTable),1);
% ACstimP2V = zeros(height(stimTable),1);
% ACpostPer = zeros(height(stimTable),1);
% ACpostP2V = zeros(height(stimTable),1);

for i = 1:height(stimTable)
   curACpre= stimTable.ACpre{i};
   curACstim = stimTable.ACstim{i};
   curACpost = stimTable.ACpost{i};
   ACcomPer(i,:) = [curACpre.period, curACstim.period, curACpost.period];
   ACcomP2V(i,:) = [curACpre.peak2VallyDiff, curACstim.peak2VallyDiff, curACpost.peak2VallyDiff];

end

stimTable.ACcomPer = ACcomPer;
stimTable.ACcomP2V = ACcomP2V;
%% Getting the liz mov for all animals + add to table
% SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
% stimTable.LM_DBt = cell(height(stimTable),1);
% for i = 1:height(stimTable)
%     % set te current rec:
%     recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
%     SA.setCurrentRecording(recName);
%     % run all required analysis:
%     if stimTable.LizMov(i) ==1
%         SA.getLizardMovements
%         LM = SA.getLizardMovements;
% 
%         DB = SA.getDelta2BetaRatio;
% 
% 
%         % calculate the number of movements to each DB bin
%         LM_DBt = zeros(size(DB.t_ms));
%         % Loop through each bin in DB and count the events in LM that fall within each bin
%         for j = 1:length(DB.t_ms)-1
%             % Count events from LM that fall within the current bin (DB(i) to DB(i+1))
%             LM_DBt(j) = sum(LM.t_mov_ms >= DB.t_ms(j)& LM.t_mov_ms < DB.t_ms(j+1));
%         end
%         % Count any events at the last bin edge
%         LM_DBt(end) = sum(LM.t_mov_ms >= DB.t_ms(end));
% 
%         %put in stimTable:
%         stimTable.LM_DBt(i) = {LM_DBt};
%         disp('LM in stimTabl')
%     end
% end

%% get the d/b during sws:
dbSWMeans = zeros([height(stimTable),3]);
dbMeans = zeros([height(stimTable),3]);
deltaSWMeans = zeros([height(stimTable),3]);
betaSWMeans = zeros([height(stimTable),3]);

for i = 1:height(stimTable)
% i =22 ;
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    DB = SA.getDelta2BetaRatio;
    % 
    stimStartT = stimTable.stimStartT(i);
    stimEndT = stimTable.stimEndT(i);

    p = 30*60*1000; %some time diff for the cycle to change.
    ACwin = 2*60*60*1000; % 4 hrs in ms
    cycWin = 2*60*60*1000;
    ACStartTimes = [stimStartT-ACwin,stimStartT+p,stimEndT+p];
    partsTimings = [stimStartT-cycWin,stimStartT+p,stimEndT+p];

    dbSW = struct('Pre',[],'Stim',[],'Post',[]);
    dbFull = struct('Pre',[],'Stim',[],'Post',[]);
    deltaSW = struct('Pre',[],'Stim',[],'Post',[]);
    betaSW = struct('Pre',[],'Stim',[],'Post',[]);

    parts = fieldnames(dbSW);  % Get a cell array of field names
    
    for j=1:numel(parts)
        %calculate AC and SC for this part:
        part = parts{j};
        SA.getDelta2BetaAC('tStart',ACStartTimes(j),'win',ACwin,'overwrite',1);
        
        SA.getSlowCycles('excludeIrregularCycles',0,'overwrite',1)
        SC = SA.getSlowCycles;
        
        curtimings = [partsTimings(j) partsTimings(j)+cycWin];
        pCyc = find(SC.TcycleOnset>= curtimings(1)& SC.TcycleOnset<=curtimings(2));
        % for SW part of the cycle:
        curCyclesOns = SC.TcycleOnset(pCyc);
        curCyclesMids = SC.TcycleMid(pCyc);
        % for Full cycles:
        curCyclrsOffs = SC.TcycleOffset(pCyc);

        for k=1:numel(curCyclesOns)
        %get the sws timings and the DB for them:
            pTmp = find(DB.t_ms>curCyclesOns(k) & DB.t_ms<curCyclesMids(k));
            dbSW.(part)(k) = mean(DB.bufferedDelta2BetaRatio(pTmp));
            betaSW.(part)(k) = mean(DB.bufferedBetaRatio(pTmp));
            deltaSW.(part)(k) = mean(DB.bufferedDeltaRatio(pTmp));

        %get the full cycles for db and betta
            pTmpF = find(DB.t_ms>curCyclesOns(k) & DB.t_ms<curCyclrsOffs(k));
            dbFull.(part)(k) = mean(DB.bufferedDelta2BetaRatio(pTmpF));
        end
        
        % get the data for regStimSham:
        if j == 1
            DB=SA.getDelta2BetaRatio;
            pre=50000;
            post=100000;
            curSleepDB = zeros([length(curCyclesOns),150]);
            for k = 1:length(curCyclesOns)
                pTmpSS = find(DB.t_ms>(curCyclesOns(k)-pre) & DB.t_ms<=(curCyclesOns(k)+post));
                curSleepDB(k,:) = DB.bufferedDelta2BetaRatio(pTmpSS);
            end

            stimTable.sleepAvg(i) = {mean(curSleepDB,1)};
        end
    end

    
    dbSWMeans(i,:) = [mean(dbSW.Pre,'omitnan') mean(dbSW.Stim,'omitnan') mean(dbSW.Post,'omitnan')];
    dbMeans(i,:) = [mean(dbFull.Pre,'omitnan') mean(dbFull.Stim,'omitnan') mean(dbFull.Post,'omitnan')];
    deltaSWMeans(i,:) = [mean(deltaSW.Pre,'omitnan') mean(deltaSW.Stim,'omitnan') mean(deltaSW.Post,'omitnan')];
    betaSWMeans(i,:) = [mean(betaSW.Pre,'omitnan') mean(betaSW.Stim,'omitnan') mean(betaSW.Post,'omitnan')];
    %save in stimTable
    stimTable.dbSW(i) = dbSW;
    stimTable.dbFull(i) = dbFull;
    stimTable.deltaSW(i) = deltaSW;
    stimTable.betaSW(i) = betaSW;
end

stimTable.dbSWMeans = dbSWMeans;
stimTable.dbMeans = dbMeans;
stimTable.deltaSWMeans = deltaSWMeans;
stimTable.betaSWMeans = betaSWMeans;
%%
stimTable.dbDiffSleep = cell(height(stimTable),1);
stimTable.dbDiffSleepM = zeros(height(stimTable),1);

firstStimInd = 50000/1000;
win = 30;
for i=1:height(stimTable)
    CurStimDur = round(stimTable.stimDuration(i)/1000);
    % calc substracted stimulation from baseline
    befSleep = stimTable.sleepAvg{i}(firstStimInd-win:firstStimInd-1);
    durSleep = stimTable.sleepAvg{i}(firstStimInd+CurStimDur-win:firstStimInd+CurStimDur-1);
    stimTable.dbDiffSleep(i) = {durSleep-befSleep};
    stimTable.dbDiffSleepM(i) = mean(durSleep,'omitnan')-mean(befSleep,'omitnan');
    
end
%% calculate the P2V for 156 sec (tril time)
fs = 156*1000;
P2Vfs = zeros(height(stimTable),3);

for i=1:height(stimTable)
    preAC = stimTable.ACpre{i};
    P2Vfs(i,1) = preAC.xcf(preAC.autocorrTimes ==fs)-preAC.xcf(preAC.autocorrTimes ==preAC.vallyPeriod);
    stimAC = stimTable.ACstim{i};
    P2Vfs(i,2) = stimAC.xcf(stimAC.autocorrTimes ==fs)-stimAC.xcf(stimAC.autocorrTimes ==stimAC.vallyPeriod);
    postAC = stimTable.ACpost{i};
    P2Vfs(i,3) = postAC.xcf(postAC.autocorrTimes ==fs)-postAC.xcf(postAC.autocorrTimes ==postAC.vallyPeriod);
end
stimTable.P2V_156 = P2Vfs;

%% save stimTable
save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');

%% Figure 1 - traces+ spikes

%% look at single traces from one night:

i = 22;
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
post=130000;
ch = 17;

%% plot trace from before stimulation:
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

% saveas (gcf, [analysisFolder filesep sprintf('singleTrialRasterPV161N13t%i.pdf',j)]);
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'cleanTraces'];
print(fileName,'-dpdf','-r300');
% saveas(f,fileName)

%% plot plot trace + raster:

 % for j=16:numel(firstTrig)
    j = 16;
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
    nClustersg = length(goodClusters);
% plot trace + raster:
    f = figure;
    h1 = subplot(2,1,1);
    % set(f, 'Position', [100, 100, 1200, 400]);
    % yyaxis left
    plot(lfp_t/1000,squeeze(lfp),'k'); hold on;
    curstims = trial(j,:) -trial(j,1) +pre ;
    xline(curstims/1000,'r','LineWidth',1.5)
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

    % saveas(gcf, [SA.currentPlotFolder filesep sprintf('DBRaterLFPT%i.pdf',j)])
    % waitforbuttonpress;

% end


% saveas (gcf, [analysisFolder filesep sprintf('singleTrialRasterPV161N13t%i.pdf',j)]);
set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'singleTrialRasterPV161N18t16'];
print(fileName,'-depsc','-vector');
% saveas(f,fileName)
%% spike rates for this recording - load data 

meanWin = 1000;
win = pre+post;
OL = 100;
load([curPhyFoler filesep 'spikeRateAll.mat'],"spikeRateAll","spikeRateT")

%% spikerate for the rec: calculate. if you run the previous block don't run this
spikeRateT = meanWin/2:OL:win-meanWin/2; %time in ms
spikeRateAll = zeros(numel(allClusters),length(spikeRateT),numel(firstTrig));
for j = 1:numel(firstTrig)
    [curSpikeRate] = getSpikeRate(spikes,allClusters,firstTrig(j)-pre,win,meanWin,OL);
    spikeRateAll(:,:,j)= curSpikeRate;

end
% save([curPhyFoler filesep 'spikeRateGoodunits.mat'],"spikeRateGM","spikeRateT",'-mat')
% load([curPhyFoler filesep 'spikeRateTrigers.mat'],"spikeRateM")
save([curPhyFoler filesep 'spikeRateAll.mat'],"spikeRateAll","spikeRateT",'-mat')

%% plot units avarage:  same night from previuos 
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
%% spike rate all red nights - initial parameters
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


%% check for or create the spike rate for all red recs (if not allready done)
for k = 1:length(spikeRecs)
    i= spikeRecs(k);
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    
    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch};
    firstTrig=stims(1:8:end-2);
    endStim=stims(8:8:end)+200;
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
        save(spikeRateFile, "spikeRateAll","spikeRateT","allClusters","goodClusters","OL",'-mat')
    end
end

%% get data for all neurons from all nights traces & baseline changes
% plot the change in spike rate betweem right before stimulations and right
% after. 
% ITI - right before - 10 seconds before ( second 10-20)
% ISI - right after - 10 seconds after 

% load(spikeRateFile, "spikeRateAll","spikeRateT")
ITIt = find(spikeRateT> 9500 & spikeRateT< 19500);

ITIsamples = [];
ISIsamples = [];
ITIsamplesN = [];
ISIsamplesN = [];

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

    unitM = mean(spikeRateAll,3); % avarage across trials.
    AllNightsUnits = [AllNightsUnits;unitM]; %add to one matrix for all units, all nights
    
    % avarage spike rate for each unit:
    unitsBaseline = mean(mean(spikeRateAll,2),3);

    % get an array for logical for all the good clusters. 
    cluster_info = readtable([curPhyFoler filesep 'cluster_info.tsv'], ...
        'FileType', 'text', 'Delimiter', '\t');
    allClusters  = cluster_info.cluster_id(find(contains(cluster_info.group,'good')|contains(cluster_info.group,'mua')));
    goodClusters = cluster_info.cluster_id(find(contains(cluster_info.group,'good')));
    goodUnits = [goodUnits ;ismember(allClusters,goodClusters)]; %indecies for the good units

    % caluclate the mean for the ITI and for ISI
    % normelize to seconds

    curITI = mean(unitM(:,ITIt),2);
    ISIt = find(spikeRateT> xPositions(end)*100+2000 & spikeRateT< xPositions(end)*100+12000);
    curISI = mean(unitM(:,ISIt),2);

    % save in the relevant vector
    ITIsamples = [ITIsamples;curITI];
    ISIsamples = [ISIsamples;curISI];
    
    % normelized: 
    curITIn = curITI./unitsBaseline; % devide each result in the avarge of the average trial 
    curISIn = curISI./unitsBaseline; 
    ITIsamplesN = [ITIsamplesN;curITIn];
    ISIsamplesN = [ISIsamplesN ; curISIn];

end

% save data:
allnightsfilename = [analysisFolder filesep 'allNightsSpikingRate.mat'];
save(allnightsfilename, "AllNightsUnits","spikeRecs","xPositions",'-mat')

ITIISImatfilename = [analysisFolder filesep 'ITIISIdata.mat'];
save(ITIISImatfilename, "ITIsamples", "ISIsamples","ITIsamplesN","ISIsamplesN","goodUnits","spikeRecs",'-mat')

%% Plot Spike rates avrages for all nights:
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


%% plot before and after (ITI/ISI) spiking rate, only good neurons. - not in paper.
fgood = figure; % only good units:
subplot(1,2,1)
goodUnitsL = logical(goodUnits);
IIdataGood = [ITIsamples(goodUnitsL), ISIsamples(goodUnitsL)];
x = [ones(length(IIdataGood),1) 2*ones(length(IIdataGood),1)];
colors = [0.6 0.49 0.97; 1 0.65 0];
swarmchart(x,IIdataGood,15,colors,'filled','XJitterWidth',0.1);
ylabel('Spikes/S')
xticks([1,2]);xticklabels(["Before Stims", "After Stims"]);
xlim([0.5 2.5])

subplot(1,2,2)
plot([1,2] ,IIdataGood,'Color',[0.7 0.7 0.7],'Marker','.','MarkerSize',8);
hold on;
plot([1,2] ,mean(IIdataGood),'Color','k', 'LineWidth',2,'Marker','.','MarkerSize',8);
xticks([1,2]);xticklabels(["Before Stims", "After Stims"]);xlim([0.5 2.5])
% grid on;
ylabel('Spikes/S')
sgtitle ('good units)')
n = sum(goodUnits); 
[pWilcoxon, ~, statsWilcoxon] = signrank(ITIsamples(goodUnitsL), ISIsamples(goodUnitsL));

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.8, 0.7, 0.03, 0.1], 'String', ...
    sprintf('Wilcoxon p-val: %.4f\n', pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
% avefigures
    set(fgood,'PaperPositionMode','auto');
    fileName=[analysisFolder filesep 'ITIISIgoodunitss'];
    print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% plot before and after (ITI/ISI) spiking rate, all units. - Fig 1G
% load data
ITIISImatfilename = [analysisFolder filesep 'ITIISIdata.mat'];
load(ITIISImatfilename, "ITIsamples", "ISIsamples","spikeRecs")


fall = figure;
subplot(1,2,1)
% colors = [0.6 0.49 0.97; 1 0.65 0];%  0.5 0.7 0.8];
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
% goodUnitsL = logical(goodUnits);
diff= ITIsamples- ISIsamples; %ISI = after, ITI =before
x = [ones(length(IIdata),1)];
colors = [ 0.5 0.7 0.8];
swarmchart(x,diff,1,colors,'filled','XJitterWidth',0.8);
ylabel('Spikes/S')
xticks([1]);xticklabels(["Before-After"]);
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

%% plot before and after (ITI/ISI) spiking rate, only good neurons. NORMELIZED: - not in paper

fgoodN = figure; % only good units:
subplot(1,2,1)
goodUnitsL = logical(goodUnits);
IIdataGoodN = [ITIsamplesN(goodUnitsL), ISIsamplesN(goodUnitsL)];
x = [ones(length(IIdataGoodN),1) 2*ones(length(IIdataGoodN),1)];
swarmchart(x,IIdataGoodN,15,colors,'filled','XJitterWidth',0.1);
ylabel('Spikes/S')
xticks([1,2]);xticklabels(["Before Stims", "After Stims"]);
xlim([0.5 2.5])

subplot(1,2,2)
plot([1,2] ,IIdataGoodN,'Color',[0.5 0.5 0.5],'Marker','.');
hold on;
plot([1,2] ,mean(IIdataGoodN),'Color','k', 'LineWidth',2,'Marker','.');
xticks([1,2]);xticklabels(["Before Stims", "After Stims"]);xlim([0.5 2.5])
% grid on;
ylabel('Spikes/S')
sgtitle ('good units - Normelized')
n = sum(goodUnits); 
[pWilcoxon, ~, statsWilcoxon] = signrank(ITIsamplesN(goodUnitsL), ISIsamplesN(goodUnitsL));

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.8, 0.7, 0.03, 0.1], 'String', ...
    sprintf('Wilcoxon p-val: %.4f\n', pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
% avefigures
    set(fgood,'PaperPositionMode','auto');
    fileName=[analysisFolder filesep 'ITIISIgoodunitsN'];
    print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% Delta 2 Beta general decrease: SWS

%% plot D2B general decrease for each part of the stimulation - SWS parts

%% plot one night: D/B decrease - Figure 2 B

i =22 ;
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
swarmchart(xswarm,yswarm,20,colorswarm,"filled",'o','XJitterWidth',0.3)
ylabel('D/B power in SWS bouts')
ylim([0 450])
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
    corrected_pvals = raw_pvals * length(raw_pvals);

    % Display results
    fprintf('\nPairwise Wilcoxon Rank-Sum Test (Mann-Whitney U):\n');
    for i = 1:length(raw_pvals)
        fprintf('%s:\t raw p = %.4f,\t Bonferroni-corrected p = %.4f\n', ...
            comparisons{i}, raw_pvals(i), corrected_pvals(i));
    end


end


%savefigure
set(fdbDec,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'DBSWSpreStimPostOneNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% plot all red nights - figure 2c

type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ~contains(stimTable.Remarks,'Ex') ...
    & ~any(isnan(stimTable.dbSWMeans),2); 
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curData = stimTable.dbSWMeans(curTrials,:);


% statistics:
[p, tbl, stats] = friedman(curData, 1,'off'); % paired data
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
if p<0.05
    % data for three groups
    before = stimTable.dbSWMeans(curTrials,1);
    during = stimTable.dbSWMeans(curTrials,2);
    after = stimTable.dbSWMeans(curTrials,3);

    % Pairwise Wilcoxon signed-rank tests
    [p_before_during, ~, stats_before_during] = signrank(before, during);
    [p_during_after, ~, stats_during_after] = signrank(during, after);
    [p_after_before, ~, stats_after_before] = signrank(after, before);
    
    raw_pvals = [p_before_during,p_during_after,p_after_before];
    num_comparisons = 3;
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons,1);
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
plot(mean(stimTable.dbSWMeans(curTrials,:),1,'omitnan'),'Color','k','LineWidth',2,'Marker','.','MarkerSize',10)
xlim([0.5, 3.5])
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

% savefigure
set(fdb,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'DBSWSredNights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% plot D/B decrease - all nights - not in paper

stimType = ["Blue","Green","Red","WhiteEx"];
stimWaveL = ["47","532","635","LED"];
% plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27], [0.5 0.5 0.5]};
numType = length(stimType);
x=1:2;
%plot Period Times:
f=figure;
set(f, 'Position', [100, 100, 800, 400]);
hold on
for type = 1:numType
    %plot the data
    subplot(1,4,type)
    curType = stimWaveL(type);
    curName = stimType(type);
    curTrials = contains(stimTable.Remarks,curType) &...
                ~contains(stimTable.Remarks,"Ex") &...
                all(~isnan(stimTable.dbSWMeans), 2);
    n = sum(curTrials);
    N = numel(unique(stimTable.Animal(curTrials)));
    % curCol = plotColors{type};
    curData = stimTable.dbSWMeans(curTrials,1:2);
    curMean = mean(stimTable.dbSWMeans(curTrials,1:2),1,'omitnan');
    statsdbSW = struct();
    %statistics:
    % 
    % [p, tbl, stats] = friedman(curData, 1,'off'); % Here, 1 indicates within-subjects design
    % fprintf('p-value for freidman ANOVA test: %.5f\n',p)
    % % % p-valure is very low, post hoc:
    % data for the four groups
    beforedb = curData(:,1);
    duringdb = curData(:,2);
    % afterdb = curData(:,3);

    % Bonferroni-corrected alpha level
    alpha = 0.05 ;

    % Pairwise Wilcoxon signed-rank tests
    [p_pre_during, ~, stats_before_during] = signrank(beforedb, duringdb);
    % [p_during_post, ~, stats_during_after] = signrank(duringdb, afterdb);
    % [p_pre_post, ~, stats_wake_during] = signrank(beforedb, afterdb);

    % Display results with Bonferroni correction
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_pre_during, alpha);
    % fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_post, alpha);
    % fprintf('Pre vs. Post: p-value = %.4f (Significant if < %.4f)\n', p_pre_post, alpha);

    % %save statistics:
    % statsdbSW.(curName).alpha = alpha;
    % statsdbSW.(curName).pAnova = p;
    % statsdbSW.(curName).p_pre_post = p_pre_post;
    % statsdbSW.(curName).p_pre_during = p_pre_during;
    % statsdbSW.(curName).p_during_post = p_during_post;


    if n>0
        
        [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
        curColorMat = animalsColors(animalIndices, :); 
        hold on
        for j = 1:height(curData)
            plot(x,curData(j,:),'Color',curColorMat(j,:),'Marker', '.')
        end
            plot(x,curMean,'color','k','LineWidth',3,'Marker', '.')
        
            % % plot ExLight:
            % ExTrials = contains(stimTable.Remarks,curType) & ...
            %     contains(stimTable.Remarks,'Ex') &...
            %     all(~isnan(stimTable.dbSWMeans),2);
            % if sum(ExTrials)>0
            %     plot(x,stimTable.dbSWMeans(ExTrials,1:2),'Color',[0.2 0.2 0.2],'LineStyle','--','Marker','.')
            % end
        hold off

        annotation('textbox', [.05 + 0.202*type 0.85, 0.03, 0.1], 'String', ...
            sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
            'right', 'VerticalAlignment', 'middle');
    end
    if type ==1
        ylabel('Time[s]')
    end
        ylim([0 600])
    xticklabels({'Pre','During'})
    xticks(1:2); xlim([0.7 2.3])
end

sgtitle ('D/B decrease during SWS bouts according to wavelangth ')


% savefigure
set(gcf,'PaperPosition',[1 1 3.5 1.5]);
fileName=[analysisFolder filesep 'DBdecreaseAllnigthscolors'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% Fig3d- plot D/B decrease - all nights DIFFS

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

%% AC changes 

%% plot sliding AC sith stimulations: Figure 1I
%set the recording:
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
AC = SA.getDelta2BetaAC;
SA.plotDelta2BetaSlidingAC ('stim',1,'stimCh',stimTable.StimTrighCh(i));

%% AC - plot specific before during and after AC - for a specific rec - Figure 1I
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
ACpre = stimTable.ACpre{i};
ACstim = stimTable.ACstim{i};
ACpost = stimTable.ACpost{i};

PDPcolors = [0.2, 0.6, 0.8; 0.9, 0.4, 0.3; 0.5, 0.8, 0.5];
ACstructs = {ACpre,ACstim,ACpost};
labels = {'preStim', 'Stim', 'postStim'};
%
for j = 1:3
    struct2vars(ACstructs{j});
    fAC = figure;
    %plot:
    lineHandles = plot(autocorrTimes/1000,real(xcf),'Color',PDPcolors(j,:),'LineWidth',4);
    ylim([-0.4 1])
    set(lineHandles(1),'MarkerSize',4);
    grid('on');
    xlabel('Period [s]');
    ylabel('Auto corr.');
    hold('on');

    plot(period/1000,real(xcf(pPeriod)),'o','MarkerSize',5,'color','k');
    text(period/1000,0.05+real(xcf(pPeriod)),num2str(period/1000));

    a = axis;
    % plot([a(1) a(1); a(2) a(2)],[xcf_bounds([1 1]) xcf_bounds([2 2])],'-b');
    plot([a(1) a(2)],[0 0],'-k');
    hold('off');
    title (labels{j})
   % save fig:
   set(fAC,'PaperPositionMode','auto');
   fileName=[analysisFolder filesep 'dbAC_ch' num2str(parDbAutocorr.ch) '_t' num2str(round(parDbAutocorr.tStart)) labels{j}];
   print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

end


%% plot AC - only Red nights: - Figure 1J
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ...
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
title ('Perios Times changes - all red nights')

% savefigure
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'ACperiodReds.pdf']);

%% Figure 3B - plot AC period time changes - According to color and animal

% animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
% plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27], [0.5 0.5 0.5]};
% numAnimal = length(animals);
numType = length(stimType);
stimData = [];
groupNum=[];
colorMat = [];
Ns = [];
ns = [];
x=1:3;
%plot Period Times:
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
    % curCol = plotColors{type};
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

%% plot stim data on 1 line:  not in paper
%run previous
pVals = zeros(numType,1);
SDs = zeros(numType,1);
means =zeros(numType,1);
stimType = ["Blue","Green","Red","LED"];
numType = 4;

newPer = 156; %sec
for i =1:numType
    data = stimData(groupNum==i);
    [p, h] = signrank(data/1000, newPer);
    pVals(i) = p;
    SDs(i) = std(data/1000);
    means(i) = mean(data/1000);

end
disp(SDs)
disp(pVals)

%test for siginificantly changes in variance:
pvartes = vartestn(stimData, groupNum, 'TestType', 'LeveneAbsolute');
%post hoc:
correctedA = 0.05/6;
[h,p12] = vartest2(stimData(groupNum==1), stimData(groupNum==2));
[h,p13] = vartest2(stimData(groupNum==1), stimData(groupNum==3));
[h,p14] = vartest2(stimData(groupNum==1), stimData(groupNum==4));
[h,p23] = vartest2(stimData(groupNum==2), stimData(groupNum==3));
[h,p24] = vartest2(stimData(groupNum==2), stimData(groupNum==4));
[h,p34] = vartest2(stimData(groupNum==3), stimData(groupNum==4));


figure;
swarmchart(groupNum,stimData/1000,10,colorMat,'filled','XJitterWidth',0.5);
hold on
scatter(1:4,means,'k','Marker','+')
ylabel('AC Period time (s)')
yline(newPer,'--','Color', [0.4 0.4 0.4])
xticks(1:4);xticklabels(stimType)



% annotation ('textbox', [0.1 0.81, 0.3, 0.05], 'String', ...
%             sprintf('Pvals: B:%.3f,G:%.3f,R:%.3f,LED:%.3f',pVals(1),pVals(2),pVals(3),pVals(4)), ...
%             'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
% 
% annotation ('textbox', [0.1 0.6, 0.5, 0.5], 'String', ...
%             sprintf('SDs: B:%f,G:%f,R:%f,LED:%f',SDs(1),SDs(2),SDs(3),SDs(4)), ...
%             'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
 
% savefigure
set(gcf,'PaperPosition',[1 1 1.8 1.2]);
fileName=[analysisFolder filesep 'ACto156'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% get the P2V at 156 sec: Figure 3C

%plot
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
x=1:3;

%plot:
f=figure;
hold on
for type = 1:numType
    %plot the data
    subplot(2,4,type)
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
    
    hold on
    for i = 1:height(curData)
        plot(x, curData(i,:),'Color',curColorMat(i,:),'Marker','.','MarkerSize',10,'LineWidth',1)
    end


    if type ==1
        ylabel('P2V in 156s')
    end
         ylim([-0.2 1.5])
    xticks(1:3);xticklabels({'Pre','Stim','Post'})
     xlim([0.7 3.3])
end

subplot(2,numType,[1:numType]+4)
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




%% StimSham analysis:

% stimtime minus pre stim - change in D/B:
% take the last 30 sec of stim and substract the 30 sec before the start of
% stim, ans avrage that. 

%% plot D/B diff stimSham - only Red nights - Figure 1H
% type = 'Red';
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ...
    ~contains(stimTable.Remarks,'Ex') & ...
    all(~isnan(stimTable.dbDiffStimM), 2) &...
    all(~isnan(stimTable.dbDiffShamM), 2); 
    
n = sum(curTrials);
animals = unique(stimTable.Animal(curTrials));
N = length(animals);
% animalsColors
% statistical tests for that figure - diff StimSham only red nights:
% trying different methods: 
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


ylim([-40 120])
xlim([0.9 2.1]);
xticks([1, 2]); % Position of the x-ticks
xticklabels({'Sham', 'Stim'}); % Labels for the x-ticks
annotation('textbox', [0.85, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.15, 0.8, 0.3, 0.2], 'String', ...
    sprintf('Wilcoxon test p =%.3f',pWilcoxon), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

set(gcf,'PaperPosition',[1 4 1.2 1.6])
saveas (gcf, [analysisFolder filesep 'DBdiffStimShamRedNights.pdf']);


%% plot the bar plot - D/B change - stim sham -  all colors - Fig 3A
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


%% Lizard Movement Data analysis:

%% Get LM DATA

% LMData = getLMData(SA, stimTable,analysisFolder);
load([analysisFolder filesep 'LMData.mat'])

%% check the correlation of mevement with th D/B
% i = 19;
% recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
% SA.setCurrentRecording(recName);
% LM_DBt= stimTable.LM_DBt{i};
% DB = SA.getDelta2BetaRatio;
% 
% %%
% figure;
% yyaxis left
% plot(DB.bufferedDelta2BetaRatio)
% hold on;
% yyaxis right
% plot(LM_DBt,'k')
%% plot the full movement data for a night - Figure 2E+F

% for one night:
i = 22;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
DB = SA.getDelta2BetaRatio;

curLMwake = LMData.LMwake(i);
curLMpre = LMData.LMpre(i);
curLMstim = LMData.LMstimbin{i};
curLMpost = LMData.LMpost(i);

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
linkaxes([a1, a2, a3, a4], 'y'); ylim([0 6.5]);

% Add 4 shared title
sgtitle('Mean movement during stimulation, PV161, Night18');

% save fig
set(gcf,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'lizMovWholeNightPV161N18'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

%% plot head movements - Red nights! - figure 2G

wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ...
    ~contains(stimTable.Remarks,'Ex');

n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));

LMpre = LMData.LMpre(curTrials, :);
LMwake = LMData.LMwake(curTrials, :);
% LMpost = LMData.LMpost(curTrials, :);
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
ylim([0 30]);

% savefigure
set(fLMr,'PaperPosition',[1 1 3 2]);
fileName=[analysisFolder filesep 'LMRednights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

% clearvars -except stimTable SA analysisFolder LMdata
%% plot head movements - All nights! - not in paper
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
numType = length(stimType);
Groups = ["Wake","Pre","During Stim","Post"];
statsHeadMovAll = struct();

% Head movements:
fLMall=figure;
sgtitle('head movement according to stim color')
for type = 1:numType
    h = subplot(2,2,type);
    %plot the data
    %curAni = animals{animal};
    curName = stimType(type);
    curType = stimWaveL(type);
    curTrials = contains(stimTable.Remarks,curType) & ...
             ~contains(stimTable.Remarks,'Ex') &...
             all(~isnan(LMData.LMpre),2); 
    n = sum(curTrials);
    N = length(unique(stimTable.Animal(curTrials)));
    % curCol = plotColors{type};

    LMpre = LMData.LMpre(curTrials, :);
    LMwake = LMData.LMwake(curTrials, :);
    LMpost = LMData.LMpost(curTrials, :);
    LMstimbinM = cell2mat(LMData.LMstimbin); % takes out the nan val
    LMstimbintrialM = mean(LMstimbinM(curTrials),2); % mean for each night

    LMplotData = [LMwake; LMpre; LMstimbintrialM; LMpost];
    x = [ones(length(LMwake),1); 2*ones(length(LMwake),1);  3*ones(length(LMwake),1); 4*ones(length(LMwake),1)];
    [~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
    curColorMat = animalsColors(animalIndices, :);
    curColorMatRe = repmat(curColorMat,4,1);
    swarmchart(x,LMplotData,15,curColorMatRe,'filled','XJitterWidth',0.1);
    xticks(1:4), xticklabels(Groups); xlim([0.7 4.3]); ylim([0 10]);
    % check the statistics:
    % Assuming data in columns where each row is a subject and each column is a timepoint
    LmStatDat =  [LMwake, LMpre, LMstimbintrialM, LMpost];
    [p, tbl, stats] = friedman(LmStatDat, 1,'off'); % Here, 1 indicates within-subjects design
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
    % 
   
    statsHeadMovAll.(curName).alpha = alpha;
    statsHeadMovAll.(curName).pAnova = p;
    statsHeadMovAll.(curName).p_wake_pre = p_wake_pre;
    statsHeadMovAll.(curName).p_wake_during = p_wake_during;
    statsHeadMovAll.(curName).p_wake_after = p_wake_after;
    statsHeadMovAll.(curName).p_pre_during = p_pre_during;
    statsHeadMovAll.(curName).p_during_after = p_during_after;
    statsHeadMovAll.(curName).p_pre_after = p_pre_after;


    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   % ylim([-40 160])
    if n==0
        plot(0,0)
    end

%     % add titles. labels...
    ylabel('Head Movements')
    title(stimType(type))

end
fileName = [analysisFolder filesep 'postHocPvalLMAllNights.mat'];
    save(fileName, 'statsHeadMovAll')


% savefigure
set(fLMall,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'LMallnights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);


%% plot polar histogram 
% all nights:
mPhasePre.movs = zeros(height(stimTable),1);
mPhasePre.DBs = zeros(height(stimTable),1);
mPhaseStim.movs = zeros(height(stimTable),1);
mPhaseStim.DBs = zeros(height(stimTable),1);
mPhasePost.movs = zeros(height(stimTable),1);
mPhasePost.DBs = zeros(height(stimTable),1);

%% create data set for all nights:
for i = 15:height(stimTable)
% i = 33;
    if stimTable.LizMov(i) ~= 1
        disp('run Lizard movement on this recording. Moving to next rec.');
        mPhasePre.movs(i) = NaN;
        mPhasePre.DBs(i) = NaN;
        mPhaseStim.movs(i) = NaN;
        mPhaseStim.DBs(i) = NaN;
        % continue; % Skip to the next iteration;
    end

    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    stimStartT = stimTable.stimStartT(i);
    stimEndT = stimTable.stimEndT(i);
    sleepStartT = stimTable.sleepStartT(i);
    sleepEndT = stimTable.sleepEndT(i);

    p = 30*60*1000; %some time diff for the cycle to change.
    preWin = stimStartT - sleepStartT; %all the time before stim start
    stimWin = stimEndT - (stimStartT+p); % stimulation period, not including first 30 min
    postWin = sleepEndT - (stimEndT+p); %post stimulations sleep, not including first 30 mins.

    nBins = 16;

    % get AC for the pre:
    SA.getDelta2BetaAC('tStart',sleepStartT ,'win',preWin , 'overwrite', 1);
    SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);
    try
    hOutPre = SA.plotLizardMovementDB('stim',1 ,'part',1,'tStartStim', ...
        stimStartT,'tEndStim',stimEndT,'nBins',nBins);
    mPhasePre.movs(i) = hOutPre.mPhaseMov;
    mPhasePre.DBs(i) = hOutPre.mPhaseDB;
    catch ME
        disp('No cycles in pre part of the recording')
    end

    % get AC for the stimulation:
    SA.getDelta2BetaAC('tStart',stimStartT+p,'win',stimWin,'overwrite', 1);
    SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);
    try    
        hOutStim = SA.plotLizardMovementDB('stim',1 ,'part',2,'tStartStim', ...
        stimStartT,'tEndStim',stimEndT,'nBins',nBins);
        mPhaseStim.movs(i) = hOutStim.mPhaseMov;
        mPhaseStim.DBs(i) = hOutStim.mPhaseDB;
    catch ME
        disp('No cycles in stim part of the recording')
        disp(ME.message);
    end

    % get AC for the Post:
    SA.getDelta2BetaAC('tStart',stimEndT+p,'win',postWin,'overwrite', 1);
    SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);
    try
    hOutPost = SA.plotLizardMovementDB('stim',1 ,'part',3,'tStartStim', ...
        stimStartT,'tEndStim',stimEndT,'nBins',nBins);
    mPhasePost.movs(i) = hOutPost.mPhaseMov;
    mPhasePost.DBs(i) = hOutPost.mPhaseDB;
    catch ME
        disp('No cycles in post part of the recording')
        disp(ME.message);
    end

    close all
end

%% save histogram data
fileName = [analysisFolder filesep 'polarhistoAllNights.mat'];
save(fileName,"mPhasePre" ,"mPhaseStim","mPhasePost")
%% load data:
load([analysisFolder filesep 'polarhistoAllNights.mat'])
%% plot polar histogram - red nights
uniqueAnimals = unique(stimTable.Animal);
relativePhasePre = mPhasePre.movs -mPhasePre.DBs;
relativePhaseStim = mPhaseStim.movs -mPhaseStim.DBs;
relativePhasePost = mPhasePost.movs -mPhasePost.DBs;

% pVal = cellfun(@(x) ~isempty(x),stimTable.LM_DBt);
% pVal = mPhasePost.movs ~= 0;
wavelength = '635';
pVal = mPhasePost.movs ~= 0 & contains(stimTable.Remarks,wavelength)...
    & ~contains(stimTable.Remarks,'Ex');

fMOVdbred=figure;
h1=subplot(1,3,1,polaraxes);hold on;
title('Pre Stimulations')
Rlim=0.5;

hP={};
for i=1:numel(uniqueAnimals)
    p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
    hP{i}=polarplot([relativePhasePre(p)';relativePhasePre(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
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
relativeStimMean = circ_mean(relativePhaseStim(p));
hP={};
for i=1:numel(uniqueAnimals)
    p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
    hP{i}=polarplot([relativePhaseStim(p)';relativePhaseStim(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
end
hold on;
hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);
hP3=polarplot([relativeStimMean relativeStimMean],[0,Rlim],'color','k','LineWidth',3);

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
    hP{i}=polarplot([relativePhasePost(p)';relativePhasePost(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
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
% set(f,'PaperPositionMode','auto');
fileName=[analysisFolder filesep 'PolarMovDBredNight'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);% clearvars -except stimTable SA analysisFolder LMdata

%% plot polser histogram - red nights - combined - not in paper
relativePhasePre = mPhasePre.movs -mPhasePre.DBs;
relativePhaseStim = mPhaseStim.movs -mPhaseStim.DBs;
relativePhasePost = mPhasePost.movs -mPhasePost.DBs;

% pVal = cellfun(@(x) ~isempty(x),stimTable.LM_DBt);
% pVal = mPhasePost.movs ~= 0;
wavelength = '635';
pVal = mPhaseStim.movs ~= 0 & contains(stimTable.Remarks,wavelength) ...
    & ~contains(stimTable.Remarks,'Ex');
n = sum(pVal);
N = numel(unique(stimTable.Animal(pVal)));

fMOVdbred=figure;
h1=polaraxes;hold on;
title('Red nights, Pre,during and Post')
Rlim=0.5;
plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27]};
a = 0.5;
hP={};
% for i=1:numel(plotColors)
p=find(pVal);
relativePreMean = circ_mean(relativePhasePre(p));
relativeStimMean = circ_mean(relativePhaseStim(p));
relativePostMean = circ_mean(relativePhasePost(p));
% hP{1}=polarplot([relativePhasePre(p)';relativePhasePre(p)'],[zeros(1,numel(p)); ...
%         Rlim*ones(1,numel(p))],'color',[plotColors{1} a],'LineWidth',1);
% hP{1}=polarplot([relativePreMean relativePreMean],[0,Rlim],'color',plotColors{1},'LineWidth',3);
hP{2}=polarplot([relativePhaseStim(p)';relativePhaseStim(p)'],[zeros(1,numel(p)); ...
        Rlim*ones(1,numel(p))],'color',[plotColors{2} a],'LineWidth',1);
hP{2}=polarplot([relativeStimMean relativeStimMean],[0,Rlim],'color',plotColors{2},'LineWidth',3);
% hP{3}=polarplot([relativePhasePost(p)';relativePhasePost(p)'], ...
%         [zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',[plotColors{3} a],'LineWidth',1);
% hP{3}=polarplot([relativePostMean relativePostMean],[0,Rlim],'color',plotColors{3},'LineWidth',3);

hold on;
hP3=polarplot([0 0],[0 Rlim],'color','k','linewidth',3);

hRose=polarhistogram(h1,relativePhaseStim(pVal),12,'Normalization','probability');
hRose.FaceColor=[0.7 0.7 0.7];
hRose.FaceAlpha=0.5;

text(0.2, Rlim/2, '\delta/\beta');
h1.ThetaTick=[0:90:330];
h1.RTick=[0.1:0.1:0.4];
%h2.ThetaTickLabels([2 3 5 6 8 9 11 12])=cell(size([2 3 5 6 8 9 11 12]));

% l1=legend([hP{1}(1),hP{2}(1),hP{3}(1),hP3,hRose],{'pre','during','post','\delta/\beta','Prob.'},'box','off');
% l1.Position=[0.7386    0.8238    0.2125    0.1190];


% statistics - using circular statistics.

% 1/ test if the data is skewed - using Rayleigh test:
%test on all data combined (all data from red night)

% get only red n
curRelativePre = relativePhasePre(pVal);
curRelativeStim = relativePhaseStim(pVal);
curRelativePost = relativePhasePost(pVal);

% Test uniformity of each condition
p_skewPre = circ_rtest(curRelativePre);
p_skewStim = circ_rtest(curRelativeStim);
p_skewPost = circ_rtest(curRelativePost);

disp(['Uniformity p-values: Before=', num2str(p_skewPre), ...
      ', During=', num2str(p_skewStim), ', After=', num2str(p_skewPost)]);

% 2/ Test if the groups are different from one another/
% use Paired Test for Circular Data  the paired-sample Watson-Williams test 

% Pairwise Watson-Williams tests
[p_bd, table_bd] = circ_wwtest(curRelativePre, curRelativeStim); % Before vs. During
[p_ba, table_ba] = circ_wwtest(curRelativePre, curRelativePost);  % Before vs. After
[p_da, table_da] = circ_wwtest(curRelativeStim, curRelativePost);  % During vs. After

% Display p-values
disp(['Watson-Williams p-value (Before vs During): ', num2str(p_bd)]);
disp(['Watson-Williams p-value (Before vs After): ', num2str(p_ba)]);
disp(['Watson-Williams p-value (During vs After): ', num2str(p_da)]);

annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.85, 0.3, 0.1], 'String', ...
    sprintf('Uniformity: Pre=%.4f, Stim=%.4f, Post=%.4f',p_skewPre,p_skewStim,p_skewPost), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');
annotation('textbox', [0.1, 0.65, 0.3, 0.1], 'String', ...
    sprintf('WWtest: Pre-Stim:%.4f, Stim-Post=%.4f, Pre-Post=%.4f',p_bd,p_da,p_ba), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');



% savefigure

set(gcf, 'PaperUnits', 'inches');         % Set paper units to inches
set(gcf, 'PaperSize', [6,4]);            % Set the paper size (width x height in inches)
set(gcf, 'PaperPosition', [0, 0, 6,4]);  % Set position on paper to match size exactly

% Print to PDF
fileName=[analysisFolder filesep 'PolarMovDBredNightsCom'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);
%% plot polser histogram - All nights
uniqueAnimals = unique(stimTable.Animal);
relativePhasePre = mPhasePre.movs -mPhasePre.DBs;
relativePhaseStim = mPhaseStim.movs -mPhaseStim.DBs;
relativePhasePost = mPhasePost.movs -mPhasePost.DBs;

% pVal = cellfun(@(x) ~isempty(x),stimTable.LM_DBt);
pVal = mPhasePost.movs ~= 0;
% wavelength = '635';
% pVal = mPhasePost.movs ~= 0 & contains(stimTable.Remarks,wavelength);


fMOVdb=figure;
h1=subplot(1,3,1,polaraxes);hold on;
title('Pre Stimulations')
Rlim=0.5;

hP={};
for i=1:numel(uniqueAnimals)
    p=find(pVal & strcmp(stimTable.Animal,uniqueAnimals(i)));
    hP{i}=polarplot([relativePhasePre(p)';relativePhasePre(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
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
    hP{i}=polarplot([relativePhaseStim(p)';relativePhaseStim(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
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
    hP{i}=polarplot([relativePhasePost(p)';relativePhasePost(p)'],[zeros(1,numel(p));Rlim*ones(1,numel(p))],'color',animalsColors(i,:),'LineWidth',1);
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
fileName=[analysisFolder filesep 'PolarMovDBAllnights'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)]);

% clearvars -except stimTable SA analysisFolder LMdata


%% Head angle to floor analysis: - supplamenrty 5
% 2: go over all records
%
% one night:
% PV159, night34
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


% plot one liner:
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

%% prepare head angle data for plots:
% load([analysisFolder filesep 'HeadAngleAvg.mat'])
% load([analysisFolder filesep 'HeadAngleSD.mat'])
% load([analysisFolder filesep 'LMData.mat'])
HeadAngleAvg = LMData.HeadAngleAvg;
headAngleSD = LMData.headAngleSD;

%% plot Head Angle - RED NIGHTS ONLY - not in paper
% headAngDiff = diff(HeadAngleAvg,[],2);
% set the zero to 90 Deg, according to accelerometer data ( this is the z
% axis, when it is 90 the accelerometer is penpendicular to the ground)
% HeadAngleAvgP = HeadAngleAvg -90;


wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) &...
    ~contains(stimTable.Remarks,'Ex');
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curHeadAvg = HeadAngleAvg(curTrials,1:3);

% STATISTICAL TEST:

[p, tbl, stats] = friedman(curHeadAvg, 1); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
wakeAng = curHeadAvg(:,1);
beforeAng = curHeadAvg(:,2);
duringAng = curHeadAvg(:,3);
% afterAng = curHeadAvg(:,4);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_wake_before, ~, stats_wake_before] = signrank(wakeAng, beforeAng);
[p_before_during, ~, stats_before_during] = signrank(beforeAng, duringAng);
% [p_during_after, ~, stats_during_after] = signrank(duringAng, afterAng);
[p_wake_during, ~, stats_wake_during] = signrank(wakeAng, duringAng);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Wake vs Before: p-value = %.4f (Significant if < %.4f)\n', p_wake_before, alpha);
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
% fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('Wake vs During: p-value = %.4f (Significant if < %.4f)\n', p_wake_during, alpha);

% plot head angles:
figure;
x1 = 1:width(curHeadAvg);
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
for i= 1:height(curHeadAvg)
    plot(x1,curHeadAvg(i,:),'Color',curColorMat(i,:),'Marker','.'); hold on;
end
plot(x1, mean(curHeadAvg), 'Color','k','Marker','.','LineWidth',1.5)
xlim([0.7,3.2]); xticks(x1); xticklabels(["Wake","Sleep","Stim"])
ylabel('Avg Head Angles (Deg)')
yline(0,'--','Headstage penpendicular to floor')
title('Head Angle avg - red nights')


annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
    sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.55, 0.25, 0.1], 'String', ...
    sprintf('wake-before p = %.4f (Significant if < %.4f)', p_wake_before, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.15, 0.2, 0.25, 0.1], 'String', ...
    sprintf('before-during p-value = %.4f ', p_before_during, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
%     sprintf('during after p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('wake during p-value = %.4f', p_wake_during), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');


% savefigure
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'HeadliftsREDNights.pdf']);
%% plot Head Angle avg - RED NIGHTS ONLY"


wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ~contains(stimTable.Remarks,'Ex') ...
      & headAngleSD(:,1)>0.03; %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curHeadAvg = HeadAngleAvg(curTrials,1:3);
% STATISTICAL TEST:

[p, tbl, stats] = friedman(curHeadAvg, 1,'off'); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
wakeAng = curHeadAvg(:,1);
beforeAng = curHeadAvg(:,2);
duringAng = curHeadAvg(:,3);
% afterAng = curHeadAvg(:,4);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_wake_before, ~, stats_wake_before] = signrank(wakeAng, beforeAng);
[p_before_during, ~, stats_before_during] = signrank(beforeAng, duringAng);
% [p_during_after, ~, stats_during_after] = signrank(duringAng, afterAng);
[p_wake_during, ~, stats_wake_during] = signrank(wakeAng, duringAng);

% Display results with Bonferroni correction
fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
fprintf('Wake vs Before: p-value = %.4f (Significant if < %.4f)\n', p_wake_before, alpha);
fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
% fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
fprintf('Wake vs During: p-value = %.4f (Significant if < %.4f)\n', p_wake_during, alpha);

% plot head angles:
figure;
x1 = 1:width(curHeadAvg);
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
for i= 1:height(curHeadAvg)
    plot(x1,curHeadAvg(i,:),'Color',curColorMat(i,:),'Marker','.'); hold on;
end
plot(x1, mean(curHeadAvg), 'Color','k','Marker','.','LineWidth',1.5)
xlim([0.7,3.2]); xticks(x1); xticklabels(["Wake","Sleep","Stim"])
ylabel('Avg Head Angles (Deg)')
yline(0,'--','Headstage penpendicular to floor')
title('Head Angle avg - red nights')


annotation('textbox', [0.8, 0.85, 0.03, 0.1], 'String', ...
    sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.8, 0.4, 0.1], 'String', ...
    sprintf('p-value for Friedman ANOVA test: %.5f',p), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.1, 0.55, 0.25, 0.1], 'String', ...
    sprintf('wake-before p = %.4f (Significant if < %.4f)', p_wake_before, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.15, 0.2, 0.25, 0.1], 'String', ...
    sprintf('before-during p-value = %.4f ', p_before_during, alpha), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% annotation('textbox', [0.55, 0.65, 0.25, 0.1], 'String', ...
%     sprintf('during after p-value = %.4f', p_during_after), 'EdgeColor', 'none', 'HorizontalAlignment', ...
%     'right', 'VerticalAlignment', 'middle');

annotation('textbox', [0.3, 0.1, 0.25, 0.1], 'String', ...
    sprintf('wake during p-value = %.4f', p_wake_during), 'EdgeColor', 'none', 'HorizontalAlignment', ...
    'right', 'VerticalAlignment', 'middle');

% savefigure
set(gcf,'PaperPositionMode','auto')
saveas (gcf, [analysisFolder filesep 'HeadliftsREDNights.pdf']);

%% plot head angles SD - RED nights - fig 2D
figure;
wavelength = '635';
curTrials = contains(stimTable.Remarks,wavelength) & ~contains(stimTable.Remarks,'Ex') ...
      & headAngleSD(:,1)>0.03; %& contains(stimTable.Animal,curAni);
n = sum(curTrials);
N = length(unique(stimTable.Animal(curTrials)));
groupNames = {'Pre', 'During', 'After'};
curHeadAvg = HeadAngleAvg(curTrials,1:3);
curHeadSD = headAngleSD(curTrials,1:3);

x1 = 1:width(curHeadSD);
[~, animalIndices] = ismember(stimTable.Animal(curTrials), uniqueAnimals);
curColorMat = animalsColors(animalIndices, :); 
for i= 1:height(curHeadSD)
    plot(x1,curHeadSD(i,:),'Color',curColorMat(i,:),'Marker','.'); hold on;
end
plot(x1, mean(curHeadSD), 'Color','k','Marker','.','LineWidth',1.5)
xlim([0.7,3.2]); xticks(x1); xticklabels(["Wake","Pre","Stim"])
ylabel('Avg Head SD')
% yline(0,'--','Headstage penpendicular to floor')
title('Head Angle SD - red nights')

% stats:
[p, tbl, stats] = friedman(curHeadSD, 1,'off'); % Here, 1 indicates within-subjects design
fprintf('p-value for freidman ANOVA test: %.5f\n',p)
% p-valure is very low, post hoc:
% Example data for three groups
wakeSD = curHeadSD(:,1);
beforeSD = curHeadSD(:,2);
duringSD = curHeadSD(:,3);
% afterAng = curHeadAvg(:,4);

% Bonferroni-corrected alpha level
alpha = 0.05 / 3;

% Pairwise Wilcoxon signed-rank tests
[p_wake_before, ~, stats_wake_before] = signrank(wakeSD, beforeSD);
[p_before_during, ~, stats_before_during] = signrank(beforeSD, duringSD);
% [p_during_after, ~, stats_during_after] = signrank(duringAng, afterAng);
[p_wake_during, ~, stats_wake_during] = signrank(wakeSD, duringSD);

    raw_pvals = [p_wake_before,p_before_during,p_wake_during];
    num_comparisons = 3;
    % Display results with Bonferroni correction
    % fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    % fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_before_during, alpha);
    % fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
    % fprintf('After vs Before: p-value = %.4f (Significant if < %.4f)\n', p_after_before, alpha);
    corrected_pvals_bonferroni = min(raw_pvals * num_comparisons, 1);
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Wake vs pre: p-value = %.4f \n', corrected_pvals_bonferroni(1));
    fprintf('pre vs During: p-value = %.4f\n', corrected_pvals_bonferroni(2));
    fprintf('Wake vs During: p-value = %.4f\n ', corrected_pvals_bonferroni(3));
% savefigure
set(gcf,'PaperPosition',[1 4 2.2 1.6])
saveas (gcf, [analysisFolder filesep 'HeadLiftsRedNightsSD.pdf']);


%% plot Head Angle - all nights, 4 subplots:

% headAngDiff = diff(HeadAngleAvg,[],2);
% set the zero to 90 Deg, according to accelerometer data ( this is the z
% HeadAngleAvgP = HeadAngleAvg -90;

animals = unique(stimTable.Animal);
stimType = ["Blue","Green","Red","LED"];
stimWaveL = ["47","532","635","LED"];
plotColors = {[0 0.586 0.9766],[0.05 0.81 0.379],[1 0.27 0.27], [0.5 0.5 0.5]};
numAnimal = length(animals);
numType = length(stimType);
statsHeadAng = struct();

% mean normelized change in D/B.
fHA=figure;
sgtitle('Head Angles over night')
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
    curHeadAvg = HeadAngleAvg(curTrials,:);
    curHeadAvgmean = mean(curHeadAvg,1,'omitnan');

    %statistics:

    [p, tbl, stats] = friedman(curHeadAvg, 1,'off'); % Here, 1 indicates within-subjects design
    fprintf('p-value for freidman ANOVA test: %.5f\n',p)
    % % p-valure is very low, post hoc:
    % data for the four groups
    wakeAng = curHeadAvg(:,1);
    beforeAng = curHeadAvg(:,2);
    duringAng = curHeadAvg(:,3);
    afterAng = curHeadAvg(:,4);

    % Bonferroni-corrected alpha level
    alpha = 0.05 / 4;

    % Pairwise Wilcoxon signed-rank tests
    [p_wake_pre, ~, stats_wake_before] = signrank(wakeAng, beforeAng);
    [p_pre_during, ~, stats_before_during] = signrank(beforeAng, duringAng);
    [p_during_after, ~, stats_during_after] = signrank(duringAng, afterAng);
    [p_wake_during, ~, stats_wake_during] = signrank(wakeAng, duringAng);

    % Display results with Bonferroni correction
    fprintf('Wilcoxon signed-rank test results with Bonferroni correction:\n');
    fprintf('Wake vs Before: p-value = %.4f (Significant if < %.4f)\n', p_wake_pre, alpha);
    fprintf('Before vs During: p-value = %.4f (Significant if < %.4f)\n', p_pre_during, alpha);
    fprintf('During vs After: p-value = %.4f (Significant if < %.4f)\n', p_during_after, alpha);
    fprintf('Wake vs During: p-value = %.4f (Significant if < %.4f)\n', p_wake_during, alpha);

    %save statistics:
    statsHeadAng.(curName).alpha = alpha;
    statsHeadAng.(curName).pAnova = p;
    statsHeadAng.(curName).p_wake_pre = p_wake_pre;
    statsHeadAng.(curName).p_wake_during = p_wake_during;
    statsHeadAng.(curName).p_pre_during = p_pre_during;
    statsHeadAng.(curName).p_during_after = p_during_after;

    %plot the data

    if n>0
        x1 = 1:4;
        for j = 1:height(curHeadAvg)
            plot(x1,curHeadAvg(j,:),'Color',curColorMat(j,:),'Marker','.'); hold on;
        end
        plot(x1, mean(curHeadAvg), 'Color','k','Marker','.','LineWidth',1.5)
        xlim([0.7,4.2]); xticks(x1); xticklabels(["Wake","Sleep","Stim","Sleep after"])
        ylabel('Avg Head Angles (Deg)')
        yline(0,'--','Headstage penpendicular to floor')
    end

    annotation('textbox', [.095 + 0.195*type, 0.85, 0.03, 0.1], 'String', ...
        sprintf('n=%i,N=%i',n,N), 'EdgeColor', 'none', 'HorizontalAlignment', ...
        'right', 'VerticalAlignment', 'middle');
   ylim([-30 50])
    if n==0
        plot(0,0)
    end

    % add titles. labels...
    ylabel('head angles from ground')
    title(stimType(type))
 
end

% savefigure
set(fHA,'Position',[50 50 1000 720]);
fileName=[analysisFolder filesep 'headLiftsAllNights-color'];
print(fileName,'-dpdf',['-r' num2str(SA.figResJPG)],'-bestfit');
save([analysisFolder filesep 'statsheadAllNights-color.mat'], "statsHeadAng")

%% check out one night: look at the angles over night/ 
i = 28;
recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
SA.setCurrentRecording(recName);
LM = SA.getLizardMovements;
[angleF, angleF_t] = getHeadLifts(LM.angles,LM.t_static_ms,100,5);
figure; plot(angleF_t/(1000*60*60),angleF-90)
title(recName)
 

%% sup figure 2 - stimSham 3 nights same animals
animal = 'PV161';
curTrials = contains(stimTable.Animal,animal);
numTrials = sum(curTrials);
for i = 1:height(stimTable)
    if curTrials(i) ==0
        continue
    else
        recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
        SA.setCurrentRecording(recName);
        plotStimSham(SA)
        sgtitle(['Type:' stimTable.Remarks(i) stimTable.recNames(i) ])
    end

end




