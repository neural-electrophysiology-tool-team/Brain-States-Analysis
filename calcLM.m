
%% calculate the lizard movment for all nights/ 

% Headstage calibration data:
% load / change accelerometer calibration values in the correct format: 
cali_result = load(['/media/sil3/Data/accelerometer_calibrations/' ...
    'headtagse_cali_recs/calibration_results.mat']).cali_result;


    for i = 32:height(stimTable)
        %set the recording:
        recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
        SA.setCurrentRecording(recName);
         disp(i)
        % get the Lizard movement data
        curHS = stimTable.Headstage(i,:);
        headstageAmpCalib = cali_result.(curHS);
        cur_sensativity = headstageAmpCalib.sensetivity;
        cur_zeroGbias = headstageAmpCalib.zeroGbais;

        SA.getLizardMovements('sensitivity', cur_sensativity(1,:)', 'zeroGBias', ...
            cur_zeroGbias(1,:)','overwrite',1)
        curLM = SA.getLizardMovements;
    end