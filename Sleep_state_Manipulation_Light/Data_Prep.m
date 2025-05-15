%% Preparation of data table
%

%% creating the StimTable:

%% get data from original table - night +light recs and put in a new table
% this part goes over all the records in SA.
%  for every recoerd that is tagged (1/2/3..) 
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWakeTest.xlsx');
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

    % set tth current rec:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    % run all required analysis:
    SA.getDelta2BetaRatio;
    SA.getDelta2BetaAC('tStart',0,'overwrite',1); %for all the rec time
    SA.getDigitalTriggers;
    
   
    % stimulations timings:
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimTable.stimStartT(i) = T.tTrig{t_ch}(1);
    stimTable.stimEndT(i) = T.tTrig{t_ch}(end);
    
    ACfull =  SA.getDelta2BetaAC;
    stimTable.sleepStartT(i) = ACfull.tStartSleep;
    stimTable.sleepEndT(i) = ACfull.tEndSleep;
         
    % get and save the stim avg
    s = getStimSham(SA,stimTable.StimTrighCh(i),1);
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
firstStimInd = 50000/1000; %time from start of array to first trig, in seconds.
win = 30; %in seconds
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
% stimTable.dbDiffStimM1 = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffStim);
% stimTable.dbDiffShamM1 = cellfun(@(x) mean(x,'omitnan'),stimTable.dbDiffSham);

%% AC - get the Data 

ACcomPer = zeros(height(stimTable),3);
ACcomP2V = zeros(height(stimTable),3);

for i = 1:height(stimTable)
   curACpre= stimTable.ACpre{i};
   curACstim = stimTable.ACstim{i};
   curACpost = stimTable.ACpost{i};
   ACcomPer(i,:) = [curACpre.period, curACstim.period, curACpost.period];
   ACcomP2V(i,:) = [curACpre.peak2VallyDiff, curACstim.peak2VallyDiff, curACpost.peak2VallyDiff];

end

stimTable.ACcomPer = ACcomPer;
stimTable.ACcomP2V = ACcomP2V;

%% get the d/b during sws:
dbSWMeans = zeros([height(stimTable),3]);
dbMeans = zeros([height(stimTable),3]);
deltaSWMeans = zeros([height(stimTable),3]);
betaSWMeans = zeros([height(stimTable),3]);

for i = 1:height(stimTable)
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    DB = SA.getDelta2BetaRatio;
    
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
        pCyc = find(SC.TcycleOnset>= curtimings(1)& SC.TcycleOnset<=curtimings(2)); %find cycles indexs
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

%% save stimTable
clearvars -except SA analysisFolder stimTable
save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');

%%