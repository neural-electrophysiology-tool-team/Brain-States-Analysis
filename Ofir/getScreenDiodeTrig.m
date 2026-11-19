 %% getScreenTrigDiode
 function [screenDiodeTriggers]=getScreenDiodeTrig(obj,varargin)

            obj.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'diodeCh',obj.recTable.ScreenDiodeCh(obj.currentPRec),@isnumeric);
            addParameter(parseObj,'win',obj.currentDataObj.recordingDuration_ms,@isnumeric);
            addParameter(parseObj,'tStart',0,@isnumeric);
            addParameter(parseObj,'lowpss_cutoff',150,@isnumeric);
            addParameter(parseObj,'diodeThreshold',3.75e5,@isnumeric); %470000   %3.8e5
            addParameter(parseObj,'minSampleDiff',1.5,@isnumeric); %minimal time between transitions, in sec %2
            addParameter(parseObj,'plotFig',1,@isnumeric);
            addParameter(parseObj,'overwrite',0,@isnumeric);

            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end

            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end

            %make parameter structure
            parTrigger=parseObj.Results;

            %check if analysis was already done
            obj.files.screenDiodeTriggers=[obj.currentAnalysisFolder filesep 'diodeTrigger_ch' num2str(diodeCh) '.mat'];
            if exist(obj.files.screenDiodeTriggers,'file') & ~overwrite
                disp('StimDiode analysis already exists for this recording');
                screenDiodeTriggers=load(obj.files.screenDiodeTriggers);
                return;
            end

            [d,d_ms] = obj.currentDataObj.getAnalogData(diodeCh,tStart,win);
            d=squeeze(d);
            fs = obj.currentDataObj.samplingFrequencyAnalog(diodeCh);
            d_lowpass = lowpass(d,lowpss_cutoff,fs);
            threshold_signal = d_lowpass>diodeThreshold;
            % threshold_signal = d>diodeThreshold;
            X = minSampleDiff*fs; % how much time in samples of 0 before transition
            n = length(threshold_signal);

            if n <= X
                screenDiodeTriggers = [];
                return;
            end

            % Get the triggers:
            % Find 1->0 transitions
            transitions = diff(threshold_signal); %[true; diff(threshold_signal) == 1]; for 0->1
            % find the indexes of transitions, and the diffs between indexes
            transition_ind = find(transitions==-1);
            % find only the transions that the first after X samples (meaning the
            % diffs between indexs are above X)
            diffs = [diff(transition_ind)]; %[0;diff(transition_ind)];
            % large_trans_ind = transition_ind(diffs >X); % find stat phases
            % screenDiodeTriggers_stat = large_trans_ind(1:2:end-1);
            % for most animals, MT
            %%% idx_stat = find (diffs > 19*20000); idx_move = find (diffs > 0.85*20000 & diffs < 2*20000) ;
            % for PV167, MT
            %%% idx_move = (find (diffs > 19*20000))+1; idx_stat = find ((diffs > 17.5*20000 & diffs < 19*20000)) +1;

            %%%screenDiodeTriggers_stat = transition_ind(idx_stat);
            %%% screenDiodeTriggers_move = transition_ind(idx_move); 

            % for MT_W
            idx_stat = find (diffs > 19*20000); idx_move = (find (diffs > 19*20000))+1;

            screenDiodeTriggers_stat = transition_ind(idx_stat); % -2000 minus 100 ms to get the ~start of the square blinking, rather than its end
            screenDiodeTriggers_move = transition_ind(idx_move); 

            % try
            %     screenDiodeTriggers_move = transition_ind(idx_move(find(diffs(idx_move+1) > 17*20000)));  %<27000 % find the beginning of the double flashing - move phases
            % catch
            %     idx_move(end)=[];
            %     screenDiodeTriggers_move = transition_ind(idx_move(find(diffs(idx_move+1) > 17*20000)));  %<27000 % find the beginning of the double flashing - move phases
            % end
            % trig2remove = length(idx_move) - length(screenDiodeTriggers_move);
            % screenDiodeTriggers_stat = transition_ind(idx_stat(1:(length(idx_stat)-trig2remove)));
            if length(screenDiodeTriggers_stat)~=length(screenDiodeTriggers_move)
                disp('different number of stat and move triggers');
            else
                valid = all(screenDiodeTriggers_stat < screenDiodeTriggers_move); % Check that every corresponding element of stat is smaller than move
                if ~valid
                    disp('Non-valid triggers');
                end
            end

           

            % screenDiodeTriggers_move = large_trans_ind(2:2:end-2);
            % screenDiodeTriggers_move = transition_ind(find(diffs >17*20000 & diffs <20*20000));  %<27000 % find the beginning of the double flashing - move phases
            % screenDiodeTriggers_move = transition_ind(diffs >17*20000 & diffs <19.8*20000);  %<27000 % find the beginning of the double flashing - move phases
            % % % % screenDiodeTriggers_move = [];
            % % % % for t = 1:length(screenDiodeTriggers_stat)
            % % % %     relevant_trigger = screenDiodeTriggers_stat(t);
            % % % %     time_window_idx = relevant_trigger:relevant_trigger+20*20000;
            % % % %     idxInWindow = transition_ind(transition_ind >= time_window_idx(1) & transition_ind <= time_window_idx(end));
            % % % %     if ~isempty(idxInWindow)
            % % % %         screenDiodeTriggers_move = [screenDiodeTriggers_move,idxInWindow(1)];
            % % % %     else
            % % % %         screenDiodeTriggers_stat(t) = [];
            % % % %     end
            % % % %     % trigger_move = transition_ind(diffs <18000); 
            % % % % end
            % for move_ind = 1: length(screenDiodeTriggers_move)
            %     current_move_ind = transition_ind(transition_ind==screenDiodeTriggers_move(move_ind));
            %     if transition_ind(current_move_ind+1)<X
            %         screenDiodeTriggers_move(move_ind)=[];
            %     end
            % end

            
            % % From them, remove those how do not have another stimulation within 10
            % % seconds. removes the change in light in the end of recording
            % % - for short recs:
            % minimal_diff = 200000;
            % prev = all_valid_trans_ind -[inf;all_valid_trans_ind(1:end-1)];
            % next =  [all_valid_trans_ind(2:end);inf]-all_valid_trans_ind ;
            % true_ind = (prev<=minimal_diff) | (next<=minimal_diff);
            % valid_Ind = all_valid_trans_ind(true_ind);
            % 
            % trigTimes = d_ms(valid_Ind);
            % maxDiff = 125000;
            % t_diff = diff(trigTimes);
            % 5
            % % --- Find large gaps ---
            % bigGaps = find(t_diff > maxDiff);
            % 
            % if isempty(bigGaps)
            %     % no large gaps > keep all
            %     diodeTriggers = trigTimes;
            % else
            %     % --- Handle start ---
            %     % last big gap before the 8 small-diffs pattern begins
            %     firstKeepIdx = bigGaps(find(bigGaps < numel(trigTimes)/2, 1, 'last')) + 1;
            %     if isempty(firstKeepIdx)
            %         firstKeepIdx = 1; % no junk at start
            %     end
            % 
            %     % --- Handle end ---
            %     % first big gap after the last 8 small-diffs pattern
            %     lastKeepIdx = bigGaps(find(bigGaps > numel(trigTimes)/2, 1, 'first'));
            %     if isempty(lastKeepIdx)
            %         lastKeepIdx = numel(trigTimes); % no junk at end
            %     end
            % 
            %     % --- Extract only the valid block ---
            %     diodeTriggers = trigTimes(firstKeepIdx:lastKeepIdx);
            % end



            if plotFig ==1
                figure; plot(d_ms,d_lowpass,'.'); title('original signal,lowpass')
                % figure; plot(d_ms,d); title('original signal');
                hold on; plot(d_ms(threshold_signal),d(threshold_signal)); 
               % if ~isempty(screenDiodeTriggers)
                    xline(d_ms(screenDiodeTriggers_stat),'k');
                    xline(d_ms(screenDiodeTriggers_move),'g')
                % end
            end

        % screenDiodeTriggers.Stat(end) = []; screenDiodeTriggers.Move(end) = [];

        removedIdx_stat=[]; removedIdx_move=[];

        a=1;
        
        % screenDiodeTriggers_stat(30)=[];
        % screenDiodeTriggers_move(30:31)=[];
        % removedIdx_stat=[30];
        % removedIdx_move=[30,31];


        valid = all(screenDiodeTriggers_stat < screenDiodeTriggers_move); % Check that every corresponding element of stat is smaller than move
        if ~valid
            disp('Non-valid triggers');
        end
        

        
        diode_signal_stat = []; diode_signal_move = [];
        for l=1:length(screenDiodeTriggers_stat)
            d_stat = (d(screenDiodeTriggers_stat(l)-0.5*2e4:screenDiodeTriggers_stat(l)+1.5*2e4-1))';
            d_move = (d(screenDiodeTriggers_move(l)-0.5*2e4:screenDiodeTriggers_move(l)+1.5*2e4-1))';
            diode_signal_stat = [diode_signal_stat;d_stat];
            diode_signal_move = [diode_signal_move;d_move];
        end
        figure; plotShifted(diode_signal_stat'); hold on; xline(0.5*2e4+1,'--');
        figure; plotShifted(diode_signal_move'); hold on; xline(0.5*2e4+1,'--');

        b=1;

        close all; 
        screenDiodeTriggers=struct();
        screenDiodeTriggers.Move = d_ms(screenDiodeTriggers_move);
        screenDiodeTriggers.Stat = d_ms(screenDiodeTriggers_stat);

        save(obj.files.screenDiodeTriggers,'screenDiodeTriggers','parTrigger','removedIdx_stat','removedIdx_move');


        end