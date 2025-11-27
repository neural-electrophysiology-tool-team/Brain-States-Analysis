%% Preparation of data table

%% creating the StimTable:

%% get data from original table - night +light recs and put in a new table
% this part goes over all the records in SA.
%  for every recoerd that is tagged (1/2/3..) 
SA=sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
% maniRecs = SA.recTable.Mani>0; % taking all the rows with manipulation
maniRecs = SA.recTable.Mani>0; % taking all the rows with manipulation
stimTable = SA.recTable(maniRecs,{'Animal','recNames','Remarks','Mani','LizMov','StimTrighCh','Headstage','stimDiodeCh','spikes','eyeVideo'});  % creating new table
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
    stimTrigs = SA.getStimTriggers;
    stimTable.stimStartT(i) = stimTrigs(1);
    stimTable.stimEndT(i) = stimTrigs(end);
    s = getStimSham(SA);
  
       
    ACfull =  SA.getDelta2BetaAC;
    stimTable.sleepStartT(i) = ACfull.tStartSleep;
    stimTable.sleepEndT(i) = ACfull.tEndSleep;
         
    % get and save the stim avg
    
    disp('got Stim for this rec')
    stimTable.StimAvg(i) = {mean(s.StimDB,1,"omitnan")};
    stimTable.StimAvgSham(i) = {mean(s.StimDBSham,1,"omitnan")};
    stimTable.times(i) = {s.ts};
    stimTable.stimDuration(i) = s.stimDuration;
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
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
clearvars -except SA analysisFolder stimTable
save([analysisFolder filesep 'stimTableAll.mat'], "stimTable",'-mat');

%% calculate Lizard movement: getLMData

% Movement during stimulation - All Nights:
% get data:
LMwake = zeros(height(stimTable),1);
LMpre = zeros(height(stimTable),1);
LMstim = cell(height(stimTable),1);
LMstimbin = cell(height(stimTable),1);
LMpost = zeros(height(stimTable),1);
LMallMean = zeros(height(stimTable),1);
LM_DBts = cell(height(stimTable),1);

stimlength = 150;
post =stimlength*1000;
binSize = 10;

% Head Angles Data - stationary periods:
HeadAngleAvg = zeros(height(stimTable),4);
headAngleSD = zeros(height(stimTable),4);

cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;
%

for i = 1:height(stimTable)
    %set the recording:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    DB = SA.getDelta2BetaRatio;
    % get the Lizard movement data
    curHS = stimTable.Headstage{i};
    headstageAmpCali = cali_result.(curHS);
    cur_sensativity = headstageAmpCali.sensetivity;
    cur_zeroGbias = headstageAmpCali.zeroGbais;

    SA.getLizardMovements('sensitivity', cur_sensativity(1,:)', 'zeroGBias', ...
        cur_zeroGbias(1,:)') %calculates the movements according to HS calibration/
    LM = SA.getLizardMovements;

    startSleepT = stimTable.sleepStartT(i);
    endSleepT = stimTable.sleepEndT(i);

    % stimulations timings:
    stimStartT = stimTable.stimStartT(i);
    stimEndT = stimTable.stimEndT(i);
    
    if ~isnan(stimTable.StimTrighCh(i))
        t_ch = stimTable.StimTrighCh(i);
        T=SA.getDigitalTriggers;
        stims = T.tTrig{t_ch};
        
    elseif ~isempty(stimTable.stimDiodeCh)
        trig = SA.getStimDiodeTrig;
        stims = trig.diodeTriggers;
        
    end
    firstTrig=stims(1:8:end-2);
    %% movement calculations:
    % change movement data to DB time scale:
    % calculate the number of movements to each DB bin
    % debug =1; 
    debug =0;
    if debug ==0
        LM_DBt = zeros(size(DB.t_ms));
        % Loop through each bin in DB and count the events in LM that fall within each bin
        for j = 1:length(DB.t_ms)-1
            % Count events from LM that fall within the current bin (DB(i) to DB(i+1))
            LM_DBt(j) = sum(LM.t_mov_ms >= DB.t_ms(j)& LM.t_mov_ms < DB.t_ms(j+1));
        end
        % Count any events at the last bin edge
        LM_DBt(end) = sum(LM.t_mov_ms >= DB.t_ms(end));
        LM_DBts(i) = {LM_DBt};
    else
        LM_DBt = LMdata.LM_DBts{i};
    end

    LMallMean(i) = mean(LM_DBt); %general mean movement for all night- fo rnorm
     
    win = 60*60*1000; %1 hr in ms
    p = 5*60*1000; % buffer time between stages
    if (win+p)<startSleepT %make sure it is in the wake time zomd
        wakind = find(DB.t_ms> p & DB.t_ms< (p+win));
    else
        wakind = find(DB.t_ms> p & DB.t_ms < startSleepT);
    end
    LMwake(i) = mean(LM_DBt(wakind));

    befind = find(DB.t_ms>stimStartT-p-win & DB.t_ms<stimStartT-p); % 1 hour before stimulation
    LMpre(i) = mean(LM_DBt(befind));

    if stimEndT +p + win < endSleepT
        postind = find(DB.t_ms>stimEndT+p & DB.t_ms<stimEndT+p+win);
    else
        postind = find(DB.t_ms>stimEndT+p & DB.t_ms<endSleepT);
    end
    LMpost(i) = mean(LM_DBt(postind));
    
    % stimulations:
    stimLM = zeros(numel(firstTrig),stimlength);
    stimLMbin = zeros(numel(firstTrig),stimlength/binSize);
    
    for j=1:numel(firstTrig)
        pTmp=find(DB.t_ms>(firstTrig(j)) & DB.t_ms<=(firstTrig(j)+post));
        if length(pTmp) ~=stimlength
            pTmp = pTmp(1:stimlength);
        end
        curstimLM = LM_DBt(pTmp);
        stimLM(i,:) = curstimLM;
        % Reshape the array to 15x10
        curReshLM = reshape(curstimLM, binSize,stimlength/binSize);

        % Sum along the rows (dim=1) to get a 1x15 array
        curStimLMBin = sum(curReshLM, 1);
        stimLMbin(j,:) = curStimLMBin/(binSize); %normelized to sec - to bin size.

    end
    % mean across trials.
    LMstim(i) = {mean(stimLM,1)}; 
    LMstimbin(i) = {mean(stimLMbin,1)};



    %% calculate lizards head angles:
    if isfield(LM, 'angles') == 0
            disp('No Angle data after re-running getLizardMovement')
            continue
    end
    
    % calculate the angles from the LM output:
    pitchAngles = LM.angles(2,:);% the angle of the headstage from floor
    if i>9
        pitchAngles = -pitchAngles; %headstage was fliped towards the tail in those recs
    end
    [angleF, angleF_t] = getHeadLifts(pitchAngles,LM.t_static_ms,100,5);

    % get the angles avgs for each part:
    % check wake times: if it started the sleep before 1 hr after the
    % start of the recoring, make the wake from start of the recording
    % until the start of the sleep:

    wakeEnd = 1000*60*60*1;
    if wakeEnd>stimTable.sleepStartT(i)
        wakeEnd = stimTable.sleepStartT(i);
    end
    p = 1000*60*45;
    parts ={[0, wakeEnd],[stimTable.sleepStartT(i)+p,stimTable.stimStartT(i)],...
        [stimTable.stimStartT(i),stimTable.stimEndT(i)], ...
        [stimTable.stimEndT(i),stimTable.sleepEndT(i)]};
    numParts = numel(parts);
    
    for j = 1:numParts
        curPart = parts{j};
        pTmp = find(angleF_t>curPart(1)&angleF_t<curPart(2));
        curAng = angleF(pTmp);
        meanCurAng = mean(curAng);
        curSD = std(curAng);
        HeadAngleAvg(i,j) = meanCurAng;
        headAngleSD(i,j) = curSD;
    end
       isplot = 0;
    if isplot
        figure;
        plot(angleF_t/(1000*60*60), angleF,'k');
        xlabel('Time (hours)'); ylabel('Head Angle')
        xline(stimTable.sleepStartT(i)/(1000*60*60),'k')
        xline(stimTable.stimStartT(i)/(1000*60*60),'r')
        xline(stimTable.stimEndT(i)/(1000*60*60),'Color','r');
        xline(stimTable.sleepEndT(i)/(1000*60*60),'b')
    end
    

fprintf('Head angles avgrages for this recordings: %f,%f, %f, %f',[HeadAngleAvg(i,:)])

end
%%
fileName = [analysisFolder filesep 'LMDataAll.mat'];
maniRecs = SA.recTable.Mani==5; % taking all the rows with manipulation
LMData = SA.recTable(maniRecs,{'Animal','recNames','Remarks','Mani','LizMov'});  % creating new table
LMData = [LMData table(LM_DBts, LMwake, LMpre, LMstimbin ,LMpost, LMallMean,HeadAngleAvg,headAngleSD,...
               'VariableNames', {'LM_DBts','LMwake', 'LMpre', 'LMstimbin','LMpost', ...
               'LMallMean','HeadAngleAvg','headAngleSD'})];

save(fileName,"LMData",'-mat')

%% Ploar data preperation:
wavelength = 'white';
curTrials = (contains(stimTable.Remarks,wavelength)|contains(stimTable.Remarks,'DayTime')) & ...
    ~contains(stimTable.Remarks,'Ex');
recsInd = find(curTrials);

for i = 1:height(stimTable)
    if stimTable.LizMov(i) ~= 1
        disp('run Lizard movement on this recording. Moving to next rec.');
        mPhasePre.movs(i) = NaN;
        mPhasePre.DBs(i) = NaN;
        mPhaseStim.movs(i) = NaN;
        mPhaseStim.DBs(i) = NaN;
        % continue; % Skip to the next iteration;
    end
    
    if ~ismember(recsInd,i)
        continue
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
    postWin = sleepEndT - (stimEndT); %post stimulations sleep, not including first 30 mins.

    nBins = 16;

    % % get AC for the pre:
    % SA.getDelta2BetaAC('tStart',sleepStartT ,'win',preWin , 'overwrite', 1);
    % SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);
    % try
    % hOutPre = SA.plotLizardMovementDB('stim',1 ,'part',1,'tStartStim', ...
    %     stimStartT,'tEndStim',stimEndT,'nBins',nBins);
    % mPhasePre.movs(i) = hOutPre.mPhaseMov;
    % mPhasePre.DBs(i) = hOutPre.mPhaseDB;
    % catch ME
    %     disp('No cycles in pre part of the recording')
    % end

    % get AC for the stimulation:
    SA.getDelta2BetaAC('tStart',stimStartT,'win',stimWin,'overwrite', 1);
    SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);
    SA.plotSlowCycles;
    SA.plotDelta2BetaSlidingAC(stim=1);
    try    
        hOutStim = SA.plotLizardMovementDB('stim',1 ,'part',2,'tStartStim', ...
        stimStartT,'tEndStim',stimEndT,'nBins',nBins);
        SA.plotLizardMovementRaster(stim=1,part=2,nBins=nBins);
        % mPhaseStim.movs(i) = hOutStim.mPhaseMov;
        % mPhaseStim.DBs(i) = hOutStim.mPhaseDB;
    catch ME
        disp('No cycles in stim part of the recording')
        disp(ME.message);
    end
    % 
    % % get AC for the Post:
    % SA.getDelta2BetaAC('tStart',stimEndT+p,'win',postWin,'overwrite', 1);
    % SA.getSlowCycles('excludeIrregularCycles',1,'overwrite',1);
    % try
    % hOutPost = SA.plotLizardMovementDB('stim',1 ,'part',3,'tStartStim', ...
    %     stimStartT,'tEndStim',stimEndT,'nBins',nBins);
    % mPhasePost.movs(i) = hOutPost.mPhaseMov;
    % mPhasePost.DBs(i) = hOutPost.mPhaseDB;
    % catch ME
    %     disp('No cycles in post part of the recording')
    %     disp(ME.message);
    % end

    close all
end


%% save histogram data
fileName = [analysisFolder filesep 'polarhistoAllNightsNew.mat'];
save(fileName,"mPhasePre" ,"mPhaseStim","mPhasePost")


%% spike sorting - revision
SA=sleepAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolder = '/media/sil1/Data/Nitzan/Light Manipulation paper/NitzanAnalysisFiles';
% load([analysisFolder filesep 'stimTableAll.mat'])

recInd = SA.recTable.spikes ==1 & SA.recTable.Mani==4;
animals =SA.recTable.Animal(recInd);recnames = SA.recTable.recNames(recInd);

recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
for i = 1:length(recList)
    SA.setCurrentRecording(recList{i})
    SA.currentDataObj.generateChannelMapFile('40_16x2_FlexLin'); % 120_32x1_H4_CamNeuro for PV24. all layouts in TSV/electrode layouts
    binaryileName = [SA.currentDataObj.recordingDir filesep 'spikeSorting' filesep 'ch1_32.bin'];
    SA.currentDataObj.convert2Binary(binaryileName);
end

%% spikeSorting
% convert to tIc:
spikesInd = stimTable.spikes ==1;
animals =stimTable.Animal(spikesInd);recnames = stimTable.recNames(spikesInd);
recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
for i = 1:length(recList)
    SA.setCurrentRecording(recList{i});
    tIcPath = [SA.currentDataObj.recordingDir '/spikeSorting/kilosort'];
    tIc = SA.currentDataObj.convertPhySorting2tIc(tIcPath);
    
end
%% run spikeRate:
% One of them is: BuildBurstMatrix
spikesInd = stimTable.spikes ==1& (strcmp(stimTable.Remarks,'white')|contains(stimTable.Remarks,'DayTime'));
animals =stimTable.Animal(spikesInd);recnames = stimTable.recNames(spikesInd);
recList = cellfun(@(x,y) ['Animal=' x ',recNames=' y], animals, recnames, 'UniformOutput', false);
bin=1000; 
pre = 20*1000; %ms
width = 150*1000; %150 sec in ms.
meanAllRecs = [];
labelsAll= [];
tiledlayout('flow')
for i = 1:length(recList)
    SA.setCurrentRecording(recList{i});
    tIcPath = [SA.currentDataObj.recordingDir '/spikeSorting/kilosort'];
    tIc = SA.currentDataObj.convertPhySorting2tIc(tIcPath);
    stims = SA.getStimTriggers;
    trial = reshape(stims,[8,length(stims)/8])';
    firstTrig=stims(1:8:end-2);
    startTimes = round(firstTrig-pre);
    curstims = (trial(1,:) -trial(1,1));

    [M]=BuildBurstMatrix(tIc.ic,round(tIc.t/bin),round(startTimes/bin),round(width/bin));
    meanT = squeeze(mean(M,1,'omitnan'));
    meanAllRecs = [meanAllRecs;meanT];
    labelsAll= [labelsAll;tIc.label];
% 
% nexttile;
% plot(meanT'); hold on
%     meanAll=mean(meanT,1);
%     plot(meanAll,'k',LineWidth=2)
%     xline(curstims/1000+20,'r','LineWidth',1.5);
%     title(recList{i})
end

% save
save([analysisFolder filesep 'meanSpikeRate1000msbin.mat'],'meanAllRecs','labelsAll')
%% get the pre-post stim spike rates.
clearvars -except stimTable SA LMData meanAllRecs analysisFolder 
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

save([analysisFolder filesep 'prePostSpikeData.mat'],"preData","postData");


%% Eye movement data.

% get data only for stimulation timings:
stimTable.eyeMovPh = zeros(height(stimTable),1);
stimTable.eyeMovPhDB= zeros(height(stimTable),1);
for i = 1:height(stimTable)
    if stimTable.eyeMov(i) ~= 1
        disp('no Eye movement, moving to next recording.');
        continue; % Skip to the next iteration;
    end
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName)
    stims = SA.getStimTriggers;
    SA.getDelta2BetaAC('overwrite',1,'tStart', stims(1),'win',(stims(end)-stims(1)+30*60*100))
    eye = SA.getSyncedDBEyeMovements('digitalVideoSyncCh',7,'useRobustFloatingAvg',0,'overwrite',1,'stimCyclesOnly',1);
    stimTable.eyeMovPh(i) = eye.mPhaseMov;
    stimTable.eyeMovPhDB(i) = eye.mPhaseDB;

end
