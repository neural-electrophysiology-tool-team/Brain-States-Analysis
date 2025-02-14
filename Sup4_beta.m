%% Supplamentary figure 4: 

% plot the over night beta power for regular and stim nights - red and LED
SA.setCurrentRecording('Animal=PV161,recNames=Night18');
SA.plotDelta2BetaRatio('stim',1,'stimCh',15)
SA.setCurrentRecording('Animal=PV161,recNames=Night17');
SA.plotDelta2BetaRatio

SA.setCurrentRecording('Animal=PV157,recNames=Night18');
SA.plotDelta2BetaRatio('stim',1,'stimCh',15)

%% delta/beta during full cycle + beta during full cycle.

