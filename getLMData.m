% function LMdata = getLMData(SA, stimTable,analysisFolder)

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
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stims = T.tTrig{t_ch};
    firstTrig=stims(1:8:end-2);

    %% movement calculations:
    % change movement data to DB time scale:
    % calculate the number of movements to each DB bin
    debug =1; 
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

fileName = [analysisFolder filesep 'LMdata.mat'];
LMData = table(LM_DBts, LMwake, LMpre, LMstimbin ,LMpost, LMallMean,HeadAngleAvg,headAngleSD,...
               'VariableNames', {'LM_DBts','LMwake', 'LMpre', 'LMstimbin','LMpost', ...
               'LMallMean','HeadAngleAvg','headAngleSD'});

save(fileName,"LMData",'-mat')

