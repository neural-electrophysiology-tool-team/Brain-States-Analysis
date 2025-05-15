%% Preparation of data table

%% creating the StimTable:

%% get all the stim sham avg from all recs and put in a new table
% this part goes over all the records in SA.
%  for every recoerd that is tagged (1/2/3..) 
SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
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


%% save stimTable
save([analysisFolder filesep 'stimTable.mat'], "stimTable",'-mat');
