%% calculating the zeroGBais and the sensetivity of each accelerometer.


AcAn = sleepAnalysis('/media/sil3/Data/accelerometer_calibrations/headtagse_cali_recs/AccRecs.xlsx');

Headstages = unique(AcAn.recTable.Animal);
numHS = numel(Headstages);
axes = {'X','Y','Z'};
cali_result = struct();

for i = 1:numHS
    currentHS = Headstages{i};
    HSSensitivity = zeros(1,3);
    HSzeroGbias = zeros(1,3);
    %for each axes:
    for j = 1:3
        %get the positive and the negative
        %get the posotive data:
        curRecPos = ['Animal=' currentHS ',recNames=' axes{j} '_Pos'];
        AcAn.setCurrentRecording(curRecPos);
        AcAn.currentDataObj.extractMetaData;
        curPos = AcAn.currentDataObj.getAnalogData(j,0,6000); 
        curPosMean = mean(curPos,3);
        fprintf('positive mean: %.3f\n',curPosMean)
        
        % get the negative data: 
        curRecNeg = ['Animal=' currentHS ',recNames=' axes{j} '_Neg'];
        AcAn.setCurrentRecording(curRecNeg);
        AcAn.currentDataObj.extractMetaData;
        curNeg = AcAn.currentDataObj.getAnalogData(j,0,6000); 
        curNegMean = mean(curNeg,3);
        fprintf('negative mean: %.3f\n',curNegMean)
        % calculate sensetivity
        curSens = (curNegMean - curPosMean) / 2; % volts/g
        HSSensitivity(j) = curSens;
        fprintf('sensativity in %i axis: %.3f\n',j, curSens)
        
        % calculate zeroGbias
        curzGb = (curNegMean + curPosMean) / 2; %the mean of the 2 values
        HSzeroGbias(j) = curzGb;
        fprintf('zeroGbias in %i axis: %.3f\n',j, curzGb)
    
    end
    %save the parameters for each of the headstages. 
    cali_result.(currentHS).sensetivity = HSSensitivity; 
    cali_result.(currentHS).zeroGbais = HSzeroGbias;  

end