% batch run delta2beta, lizrd movement and freq for all wake
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
wakeTrials = contains(WA.recTable.recNames,'Wake') & ...
            ~contains(WA.recTable.Animal,'PV8') & ...
             isnan(WA.recTable.Exclude) & ...
             cellfun(@isempty,WA.recTable.Remarks);
wakeTable = WA.recTable(wakeTrials,:);

cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;
%% 
%--------------------------
% Changed prameters, re-run get freqBand again
% maybe also re-run get delta 2 beta. not sure yet. 
% stoped running in the meantime.
% 7.7.26
% --------------------------
for i =228:height(wakeTable)
    recName = ['Animal=' wakeTable.Animal{i} ',recNames=' wakeTable.recNames{i}];
    WA.setCurrentRecording(recName);
   
    WA.getDelta2BetaRatio;
    WA.getFreqBandDetection(...
        'fMax',          80, ...
        'segmentLength', 1000, ...
        'binDuration',   5000, ...
        'WelchOL',       0.75, ...
        'overwrite',     1,...
        'dftPoints',     2^10 ...
    );
    
        
    hs = WA.recTable.Headstage{WA.currentPRec};
    curSens = cali_result.(hs).sensetivity';
    curZeroG = cali_result.(hs).zeroGbais';
    WA.getLizardMovements(sensitivity=curSens,zeroGBias=curZeroG);

    
end