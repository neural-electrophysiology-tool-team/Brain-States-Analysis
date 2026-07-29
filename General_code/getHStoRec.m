% Go over all the recording and match the hs to record. It is either 2 or
% 4.
% % 1: get the right headstage parameters for each recording:
%   for each rec: 
%   run LM with the parameters of one of the headstages. 
%   plot the angles in Z axis.
% if good put the number of headstage in the vector, if not, try another
% and plot again, until you find the right one.
% get the head stage prarmeters
load('/media/sil3/Data/accelerometer_calibrations/headtagse_cali_recs/calibration_results.mat')
Headstages = fieldnames(cali_result);

SA=sleepAnalysis('/media/sil1/Data/Pogona Vitticeps/brainStatesWake.xlsx');
analysisFolder = '/media/sil3/Data/Pogona_Vitticeps/NitzanAnalysisFiles';
load([analysisFolder filesep 'stimTable.mat'])
load([analysisFolder filesep 'LMdata.mat'])

%% go over the recordings and plot it with an headstage. 
headstageIDperRec = zeros(height(stimTable),1);
ZGB = [cali_result.HS2.zeroGbais ; cali_result.HS4.zeroGbais]*1000000;
Sens = [cali_result.HS2.sensetivity ; cali_result.HS4.sensetivity] *1000000;

for j = 3:height(scaleInd)
    i = scaleInd(j);
    if stimTable.LizMov(i)==0
        continue
    else
        recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
        SA.setCurrentRecording(recName);
        %try with headstage2:
        SA.getLizardMovements('zeroGBias',ZGB(1,:)','sensitivity',Sens(1,:)','overwrite',1)
        LM = SA.getLizardMovements;
        figure; plot(LM.angles(3,:))

        % Wait for a key press
        disp('Press "s" or "n".');
        waitforbuttonpress; % Wait for key press
        key = get(gcf, 'CurrentCharacter'); % Get the key pressed

        % Handle key press
        if key == 's'
            disp('You pressed "s". Saving this HD as true');
            headstageIDperRec(i) = 'HS2';
            continue
        elseif key == 'n'
            disp('You pressed "n". ploting with HS4:');
            SA.getLizardMovements('zeroGBias',ZGB(2,:)','sensitivity',Sens(2,:)','overwrite',1)
            LM = SA.getLizardMovements;
            plot(LM.angles(3,:))
            
            disp('Press "s" or "n".');
            waitforbuttonpress; % Wait for key press
            key = get(gcf, 'CurrentCharacter'); % Get the key pressed
            
            if key == 's'
            disp('You pressed "s". Saving this HD as true');
            headstageIDperRec(i) = 'HS4';
            continue
            elseif key == 'n'
                disp('You pressed "n". Neither are right, moving to next recording');
                headstageIDperRec(i) = 'neither';
                continue
            end
        elseif key == 'o'
            break
        else
            disp('Key not recognized. Please press "s" or "n".');
        end
    end
end




%%

%% go over the recordings and plot it with an headstage. 
headstageIDperRec = zeros(height(stimTable),1);
ZGB = [cali_result.HS2.zeroGbais ; cali_result.HS4.zeroGbais]*1000000;
Sens = [cali_result.HS2.sensetivity ; cali_result.HS4.sensetivity] *1000000;

for j = 7:height(scaleInd)
    i = scaleInd(j);
    if stimTable.LizMov(i)==0
        continue
    else
        recName = ['Animal=' stimTable.Animal{i} ',recNames=' stimTable.recNames{i}];
        SA.setCurrentRecording(recName);
        %try with headstage2:
        SA.getLizardMovements('zeroGBias',ZGB(1,:)','sensitivity',Sens(1,:)','overwrite',1)
        LM = SA.getLizardMovements;
        figure; plot(LM.angles(3,:)); title ('HS2')

        SA.getLizardMovements('zeroGBias',ZGB(2,:)','sensitivity',Sens(2,:)','overwrite',1)
        LM2 = SA.getLizardMovements;
        figure;   plot(LM.angles(3,:)); title('HS4')
         
    end
end