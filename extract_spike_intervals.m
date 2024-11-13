function ISIt = extract_spike_intervals(spikeRateT, xPositions)
    % Extract the spike rate time intervals between stimulations, 
    % starting 1 second after stimulation and ending 500 ms before the next stimulation.
    
    ISIt = [];
    
    for i = 1:length(xPositions)
        % Find the start and end of the interval
        start= xPositions(i) + 20; % 1 second after stimulation
        finish = xPositions(i)+ 40;

    % Extract the relevant spike rate times and append to the result
    t=  find(spikeRateT >= start*100 & spikeRateT <= finish*100);
    ISIt = [ISIt,t];
    end
end