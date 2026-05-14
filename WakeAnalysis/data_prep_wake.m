%% data prep for wake states:
% get the WA:
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
%%
wakeTrials = contains(WA.recTable.recNames,'Wake') & ...
             cellfun(@isempty,WA.recTable.Remarks);
wakeTable = WA.recTable(wakeTrials,:);

%% run get lizard movement on all
cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;    
for i = 56:height(wakeTable)
    recName = ['Animal=' wakeTable.Animal{i} ',recNames=' wakeTable.recNames{i}];
    WA.setCurrentRecording(recName);
    if isnan(WA.recTable.Exclude(WA.currentPRec))
        hs = WA.recTable.Headstage{WA.currentPRec};
        curSens = cali_result.(hs).sensetivity';
        curZeroG = cali_result.(hs).zeroGbais';
        curLM = WA.getLizardMovements(sensitivity=curSens,zeroGBias=curZeroG);
        DB = WA.getDelta2BetaRatio;
    end
end