% Collect all relevant trials for anlysis.
% Created: 10.08.26 by Ofir.
% Spectral analysis: 'movWin' = 2000, 'movOLWin'= 1900.

condition = "hunter";
Animals =["PV213"]; %["PV159","PV161","PV157","PV162","PV126","PV143","PV106","PV153","PV208","PV213","PV228","PV167"];

Nitzan_experiments = 0;
if (Nitzan_experiments)
    SA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx');
else
    SA = wakeAnalysis('/media/sil3/Data/Ofir/BrainStatesExperimentsOfir.xlsx');
end

switch condition
    case "hunter"
        T = table([], [], [], [], [], [], [], ...
            'VariableNames', {'Animal','Record', 'Trial_num','Triger_ms',...
            'Type', 'End_ms', 'Strike_ms'});

        curHunTrials = contains(SA.recTable.recNames,'Hunter') & contains(SA.recTable.Animal,Animals) & isnan(SA.recTable.Exclude) & SA.recTable.Remarks~="Sim";

    case "sound"
        T = table([], [], [], [],  ...
            'VariableNames', {'Animal','Record', 'Trial_num','Triger_ms' ...
            });

        curHunTrials = contains(SA.recTable.recNames,'Wake') & contains(SA.recTable.Animal,Animals) & isnan(SA.recTable.Exclude) & SA.recTable.Remarks=="NS" & SA.recTable.Remarks~="Sim";
        % curHunTrials = (contains(SA.recTable.recNames,'Wake') | contains(SA.recTable.recNames,'Hunter')) & contains(SA.recTable.Animal,Animals) & isnan(SA.recTable.Exclude) & SA.recTable.Remarks~="Sim";

    case "meadow"
        T = table([], [], [], [], [],  ...
            'VariableNames', {'Animal','Record', 'Trial_num','Triger_ms' ...
            ,'Move_ms'}); % 'Triger_ms' = start meadow static, 'Move_ms' = movement phase after static phase (i.e. the same trial)

        curHunTrials = contains(SA.recTable.recNames,'Wake') & contains(SA.recTable.Animal,Animals) & isnan(SA.recTable.Exclude) & SA.recTable.Remarks=="MT";

end

recordsTable = SA.recTable(curHunTrials,:);
hunter2131 = T;

%%
for n = 28:height(recordsTable)

    SA.setCurrentRecording(['Animal=' recordsTable.Animal{n} ',recNames=' recordsTable.recNames{n}]);

% for row = 1:height(reward020926)
% 
%     animal = reward020926.Animal(row); record = reward020926.Record(row); s = reward020926.Sound(row); %for sound    screenDiode for meadow
%     n = find(recordsTable.Animal == animal & recordsTable.recNames == record);
%     SA.setCurrentRecording(['Animal=' recordsTable.Animal{n} ',recNames=' recordsTable.recNames{n}]);

    ch = SA.recTable.defaulLFPCh(SA.currentPRec); 

    win = 5000; % time before triger (ms)

    switch condition

        case "hunter"
            % Get the type of bug movement _for hunters)
            videoPath = recordsTable.VideoFiles{n};        
            infoPath = fullfile(fileparts(fileparts(videoPath)), 'info.yaml');        
            txt = fileread(infoPath);
            tokens = regexp(txt, 'movement_type:\s*(.*?)notes:', 'tokens', 'once');
            movement_type = strtrim(tokens{1});
               
            arenaCSVFile=[SA.currentAnalysisFolder filesep 'ArenaCSV.mat'];
            if exist(arenaCSVFile,'file')
                load(arenaCSVFile); %load data
            else
                disp('running getarenaCSVs')
                SA.getArenaCSVs;
                load(arenaCSVFile); %load data
            end
         
            if ~isfield(arenaCSVs, 'startTrigSh')
                disp('running getarenaCSVs')
                SA.getArenaCSVs;
                load(arenaCSVFile); %load data
            end

            zeroTimes = arenaCSVs.startTrigSh; 
            endTrialT = arenaCSVs.endTrigSh;
            strikeTimeVec = nan(length(zeroTimes),1);
            try
                % strikeTrialNum =arenaCSVs.strikeTrialNum;
                % strikeT = arenaCSVs.strikeTrigSh;
                % strikeTimeVec(strikeTrialNum) = strikeT';


                strikeT = arenaCSVs.strikeTrigSh;
                strikeT = strikeT(:);

                if isfield(arenaCSVs, 'strikeTrialNum')
                    strikeTrialNum = arenaCSVs.strikeTrialNum;
                else
                    % Find the trial of each strike by timing: zeroTimes(k) <= strikeT <= endTrialT(k)
                    strikeTrialNum = nan(size(strikeT));
                    for s = 1:numel(strikeT)
                        k = find(strikeT(s) >= zeroTimes & strikeT(s) <= endTrialT, 1);
                        if ~isempty(k), strikeTrialNum(s) = k; end
                    end
                end
                
                % Keep only strikes that fall inside a trial
                valid = ~isnan(strikeTrialNum);
                if any(~valid)
                    warning('%d strike(s) did not fall within any trial and were skipped.', nnz(~valid));
                end
                
                strikeTimeVec = nan(numel(zeroTimes), 1);
                strikeTimeVec(strikeTrialNum(valid)) = strikeT(valid)';
            % catch
            %     strikeTimeVec = nan(length(zeroTimes),1);
            end
            

            fullwin = 60000; % overall time for trial (ms)

        case "sound"
            

            template = template_RS_PV167_RL3;
            reward = 0;
            microphoneTriggers=getStimMicrophoneTrig(SA,template-mean(template),reward);
            s = microphoneTriggers.microphoneTriggers.Sound;
            zeroTimes = s; 

            fullwin = 10000; % overall time for trial (ms)

        case "meadow"
            
            diodeCh=SA.recTable.ScreenDiodeCh(SA.currentPRec);
            diodeFile=[SA.currentAnalysisFolder filesep 'diodeTrigger_ch' num2str(diodeCh) '.mat'];       
            load(diodeFile); %load data
            thresh = parTrigger.diodeThreshold;
            [screenDiodeTriggers]=getScreenDiodeTrig(SA,"diodeThreshold",thresh);
            move = screenDiodeTriggers.screenDiodeTriggers.Move;
            stat = screenDiodeTriggers.screenDiodeTriggers.Stat;
            zeroTimes = stat; 

            fullwin = 50000; % overall time for trial (ms)
    end




    numTrials = length(zeroTimes);
    
    trialLabels = strings(numTrials, 1);
    for t = 1:numTrials
        trialLabels(t) = sprintf('%d/%d', t, numTrials);
    end
    
    switch condition
    
        case "hunter"   
            newRows = table( ...
                repmat(string(recordsTable.Animal{n}), numTrials, 1), ...
                repmat(string(recordsTable.recNames{n}), numTrials, 1), ...
                trialLabels, ...
                zeroTimes, ...
                repmat(string(movement_type), numTrials, 1), ...
                endTrialT, ...
                strikeTimeVec, ...
                'VariableNames', hunter2131.Properties.VariableNames);
    
        case "sound"
            newRows = table( ...
                repmat(string(recordsTable.Animal{n}), numTrials, 1), ...
                repmat(string(recordsTable.recNames{n}), numTrials, 1), ...
                trialLabels, ...
                zeroTimes, ...
                'VariableNames', hunter2131.Properties.VariableNames);

        case "meadow"
            newRows = table( ...
                repmat(string(recordsTable.Animal{n}), numTrials, 1), ...
                repmat(string(recordsTable.recNames{n}), numTrials, 1), ...
                trialLabels, ...
                zeroTimes', ...
                move', ...
                'VariableNames', hunter2131.Properties.VariableNames);

        end

    hunter2131 = [hunter2131; newRows];

end


%%
% T = reward020926; 
% % Find groups with the same Animal and Record
% [G, ~] = findgroups(T.Animal, T.Record);
% 
% % Number of rows in each group
% groupCounts = splitapply(@numel, T.Record, G);
% 
% % Row number within each group
% rowNum = zeros(height(T), 1);
% 
% for i = 1:max(G)
%     rowNum(G == i) = 1:sum(G == i);
% end
% 
% % Create Trial_num
% Trial_num = compose('%d/%d', rowNum, groupCounts(G));
% 
% % Create new table
% T_new = T;
% 
% % Remove TypeSound
% T_new.TypeSound = [];
% 
% % Insert Trial_num as the third column
% T_new = addvars(T_new, Trial_num, 'Before', 'Sound', ...
%     'NewVariableNames', 'Trial_num');
% 
% % Rename Sound to Triger_ms
% T_new.Properties.VariableNames{'Sound'} = 'Triger_ms';
% 
% hunter1531 = T_new;
