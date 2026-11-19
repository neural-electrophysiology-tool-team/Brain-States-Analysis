
%% getting clean signal from stim photodiode:
WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
WA.setCurrentRecording('Animal=PV106,recNames=DayStim1');
% diodeCh = WA.recTable.stimDiodeCh(WA.currentPRec);
diodeCh=3;
[d,d_ms] = WA.currentDataObj.getAnalogData(diodeCh,0,WA.currentDataObj.recordingDuration_ms);
d=squeeze(d);
%%
lowpss_cutoff  = 100;
fs = WA.currentDataObj.samplingFrequencyAnalog(diodeCh);

d_lowpass = lowpass(d,lowpss_cutoff,fs);
threshold = 470000;
T = d_lowpass>threshold;
X = 2*fs; % how much time in samples of 0 before transition

valid_trans_ind= find_transitions(T ,X);





function valid_trans_ind = find_transitions(signal ,X)
   
    
    signal = signal(:);
    n = length(signal);
    
    if n <= X
        valid_trans_ind = [];
        return;
    end
    
    % Find 0->1 transitions
    transitions = [true; diff(signal) == 1]; 
    % find the indexes of transitions, and the diffs between indexes
    transition_ind = find(transitions); 
    % find only the transions that the first after X samples (meaning the
    % diffs between indexs are above X)
    diffs = [0;diff(transition_ind)];
    all_valid_trans_ind = transition_ind(diffs >X); 
    % From them, remove those how do not have another stimulation within 10
    % seconds. removes the change in light in the end of recording.
    minimal_diff = 200000;
    prev = all_valid_trans_ind -[inf;all_valid_trans_ind(1:end-1)];
    next =  [all_valid_trans_ind(2:end);inf]-all_valid_trans_ind ;
    true_ind = (prev<=minimal_diff) | (next<=minimal_diff);
    valid_trans_ind = all_valid_trans_ind(true_ind);
end
