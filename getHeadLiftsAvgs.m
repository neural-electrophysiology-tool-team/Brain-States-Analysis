%% plot for all nights:
% get the full data for head angles in the wake time (before sleep), to 


stimTable.stimStartT = zeros(height(stimTable),1); %for the start and end of stimulation
stimTable.stimEndT = zeros(height(stimTable),1);
stimTable.sleepStartT = zeros(height(stimTable),1); % for the start and end of sleep.
stimTable.sleepEndT = zeros(height(stimTable),1);
%%
for i = 1:height(stimTable)
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);

    %get the timings parameters (and save them to stim table! finly)
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimTable.stimStartT(i) = T.tTrig{t_ch}(1);
    stimTable.stimEndT(i) = T.tTrig{t_ch}(end);
    AC= SA.getDelta2BetaAC;
    stimTable.sleepStartT(i) = AC.tStartSleep;
    stimTable.sleepEndT(i) = AC.tEndSleep;

end
%%

HeadAngleAvg = zeros(height(stimTable),4);

%%
for i = 48:height(stimTable)
    recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
    SA.setCurrentRecording(recName);

    %get the timings parameters (and save them to stim table! finly)
    t_ch = stimTable.StimTrighCh(i);
    T=SA.getDigitalTriggers;
    stimTable.stimStartT(i) = T.tTrig{t_ch}(1);
    stimTable.stimEndT(i) = T.tTrig{t_ch}(end);
    AC= SA.getDelta2BetaAC;
    stimTable.sleepStartT(i) = AC.tStartSleep;
    stimTable.sleepEndT(i) = AC.tEndSleep;

    % get the angle data
    if ~stimTable.LizMov(i)
        continue
    end
    
    LM = SA.getLizardMovements;
   
    if isfield(LM, 'angles') == 0
        SA.getLizardMovements('overwrite',1)
        LM = SA.getLizardMovements;
        if isfield(LM, 'angles') == 0
            disp('No Angle data after re-running getLizardMovement')
            continue
        end
    end
    
    % calculate the angles from the LM output:
    [angleF, angleF_t] = getHeadLifts(LM.angles,LM.t_static_ms,100,5);

    % get the angles avgs for each part:
    % check wake times: if it started the sleep before 1 hr after the
    % start of the recoring, make the wake from start of the recording
    % until the start of the sleep:

    wakeEnd = 1000*60*60*1;
    if wakeEnd>AC.tStartSleep
        wakeEnd = AC.tStartSleep;
    end

    parts ={[0, wakeEnd],[AC.tStartSleep,stimStartT],[stimStartT,stimEndT], [stimEndT,AC.tEndSleep]};
    numParts = numel(parts);
    
    for j = 1:numParts
        curPart = parts{j};
        pTmp = find(angleF_t>curPart(1)&angleF_t<curPart(2));
        curAng = angleF(pTmp);
        meanCurAng = mean(curAng);
        HeadAngleAvg(i,j) = meanCurAng;

    end
 
    isplot = 1;
    if isplot
        figure;
        plot(angleF_t/(1000*60*60), angleF,'k');
        xlabel('Time (hours)'); ylabel('Head Angle')
        xline(AC.tStartSleep/(1000*60*60),'b')
        xline(stimStartT/(1000*60*60),'r')
        xline(stimEndT/(1000*60*60),'Color','r');
        xline(AC.tEndSleep/(1000*60*60),'b')
    end


fprintf('Head angles avgrages for this recordings: %f,%f, %f, %f',[HeadAngleAvg(i,:)])


end


save([analysisFolder filesep 'HeadAngleAvg.mat'],'HeadAngleAvg','-mat')