function LMdata = getLMData(SA, stimTable,analysisFolder)

% Movement during stimulation - All Nights:
% get data:
LMwake = zeros(height(stimTable),1);
LMpre = zeros(height(stimTable),1);
LMstim = cell(height(stimTable),1);
LMstimbin = cell(height(stimTable),1);
LMpost = zeros(height(stimTable),1);
LMallMean = zeros(height(stimTable),1);

stimlength = 150;
post =150*1000;
binSize = 10;

for i = 1:height(stimTable)
    %set the recording:
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);
    DB = SA.getDelta2BetaRatio;
    LM_DBt = stimTable.LM_DBt{i};
    startSleepT = stimTable.sleepStartT(i);
    endSleepT = stimTable.sleepEndT(i);

    %check if LM analysis was already done for this Rec
    if isempty(LM_DBt)
        disp('run Lizard movement on this recording. Moving to next rec.');
        LMpre(i) = NaN;
        LMwake(i) = NaN;
        LMpost(i) = NaN;
        continue; % Skip to the next iteration;
    end

    % stimulations timings:
    stimStartT = stimTable.stimStartT(i);
    stimEndT = stimTable.stimEndT(i);
    
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch};
    firstTrig=stims(1:8:end-2);


    % get and save mov for each part
    LMallMean(i) = mean(LM_DBt(LM_DBt<450)); %general mean movement for all night- fo rnorm

    % deivided byt the general mean of the recording:
    % timeBin = DB.parDBRatio.movWin-DB.parDBRatio.movOLWin;
    win = 60*60*1000; %ms
    p = 5*60*1000; % buffer time between stages
    if win+p<startSleepT %make sure it is in the wake time zomd
        wakind = find(DB.t_ms> p & DB.t_ms< p+win);
    else
        wakind = find(DB.t_ms> p & DB.t_ms < startSleepT);
    end
    LMwake(i) = mean(LM_DBt(wakind))/LMallMean(i);

    befind = find(DB.t_ms>stimStartT-p-win & DB.t_ms<stimStartT-p); % 1 hour before stimulation
    LMpre(i) = mean(LM_DBt(befind))/LMallMean(i);

    if stimEndT +p + win < endSleepT
        postind = find(DB.t_ms>stimEndT+p & DB.t_ms<stimEndT+p+win);
    else
        postind = find(DB.t_ms>stimEndT+p & DB.t_ms<endSleepT);
    end

    LMpost(i) = mean(LM_DBt(postind))/LMallMean(i);

    stimLM = zeros(numel(firstTrig),stimlength);
    stimLMbin = zeros(numel(firstTrig),stimlength/binSize);
    % Add the last bin edge
    for j=1:numel(firstTrig)
        pTmp=find(DB.t_ms>(firstTrig(j)) & DB.t_ms<=(firstTrig(j)+post));
        %StimDB(i,:)=1./DB.bufferedDelta2BetaRatio(pTmp);
        if length(pTmp) ~=stimlength
            pTmp = pTmp(1:stimlength);
        end
        curstimLM = LM_DBt(pTmp);
        stimLM(i,:) = curstimLM;

        % Reshape the array to 15x10
        curReshLM = reshape(curstimLM, binSize,stimlength/binSize);

        % Sum along the rows (dim=1) to get a 1x15 array
        curStimLMBin = sum(curReshLM, 1);
        stimLMbin(j,:) = curStimLMBin/(binSize*LMallMean(i)); %normelized to sec - to bin size.

    end
    LMstim(i) = {mean(stimLM,1)};
    LMstimbin(i) = {mean(stimLMbin,1)};

end

fileName = [analysisFolder filesep 'LMdata.mat'];
save(fileName,"LMwake","LMpre","LMstim","LMpost","LMallMean",'-mat')
LMdata = table(LMwake, LMpre, LMstimbin ,LMpost, LMallMean,...
               'VariableNames', {'LMwake', 'LMpre', 'LMstimbin','LMpost','LMallMean'});


save(fileName,"LMdata",'-mat')

end