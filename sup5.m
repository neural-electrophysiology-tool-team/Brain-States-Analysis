%% sup 5 - head ansles no stim:

recList = {'Animal=PV149,recNames=Night5',...
    'Animal=PV149,recNames=Night8',...
    'Animal=PV159,recNames=Night4',...
    'Animal=PV159,recNames=Night11',...
    'Animal=PV159,recNames=Night31',...
    'Animal=PV161,recNames=Night4',...
    'Animal=PV161,recNames=Night12',...
    'Animal=PV161,recNames=Night22',...
    'Animal=PV157,recNames=Night7',...
    'Animal=PV157,recNames=Night15',...
    'Animal=PV157,recNames=Night30',...
    'Animal=PV162,recNames=Night9',...
    'Animal=PV162,recNames=Night20',...
    'Animal=PV162,recNames=Night30',...
    'Animal=PV126,recNames=Night4',...
    'Animal=PV126,recNames=Night9',...
    };


SA.batchProcessData('getLizardMovements',recList)



%% plot all of them:
    p = 1000*60*45;
    cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;

 
% for i = 13:length(recList)
recName = 'Animal=PV159,recNames=Night4';
    % SA.setCurrentRecording(recList{4});
    SA.setCurrentRecording(recName);
       curHS = 'HS4';%stimTable.Headstage(i,:);
        headstageAmpCalib = cali_result.(curHS);
        cur_sensativity = headstageAmpCalib.sensetivity;
        cur_zeroGbias = headstageAmpCalib.zeroGbais;

        % SA.getLizardMovements('sensitivity', cur_sensativity(1,:)', 'zeroGBias', ...
        %     cur_zeroGbias(1,:)','overwrite',1)
        LM = SA.getLizardMovements;
        pitch = -LM.angles(2,:);
    [angleF, angleF_t] = getHeadLifts(pitch,LM.t_static_ms,100,5);

    %get sleep times:
    SA.getDelta2BetaRatio;
    SA.getDelta2BetaAC;
    AC = SA.getDelta2BetaAC;
    curSleepStartT = AC.tStartSleep;
    curSleepEndT = AC.tEndSleep;

    fHA = figure;
    plot(angleF_t/(1000*60*60), angleF,'k');
    xlabel('Time (hours)'); ylabel('Head Angle')
    xline(curSleepStartT/(1000*60*60),'k')
    xline((curSleepStartT+p)/(1000*60*60),'g')
    % xline(stimTable.stimStartT(i)/(1000*60*60),'r')
    % xline(stimTable.stimEndT(i)/(1000*60*60),'Color','r');
    xline(curSleepEndT/(1000*60*60),'b')
    yline(0,'--','Color',[0.4 0.4 0.4])

    %save Figure
    set(fHA,'PaperPositionMode','auto');
    fileName=[SA.currentPlotFolder filesep 'headAngle'];
    print(fileName,'-dpdf','-r300');
        fileName2=[analysisFolder filesep 'headAngle'];
    print(fileName2,'-dpdf','-r300');

% end