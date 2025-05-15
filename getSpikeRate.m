function [spikeRate, spikeRate_t] = getSpikeRate(data, units, tStart, win, meanWin, OL, fs)
    % This function computes the spike rate for multiple specified units (neurons) in a dataset
    % across time windows with specified duration, averaging within a window size (meanWin)
    % with a specified overlap (OL).
    %
    % Inputs:
    %   - data: A structure containing spikes and timings, output from Phy after conversion to .mat
    %   - units: Array of unit IDs (e.g., neuron IDs) for which spike rates are calculated
    %   - tStart: The start time in ms (default = 0)
    %   - win: Duration of time window in ms (default = 2 hours)
    %   - meanWin: Duration of the averaging window in ms (default = 1000 ms)
    %   - OL: Overlap duration between consecutive windows in ms (default = 500 ms)
    %
    % Outputs:
    %   - spikeRate: Matrix of spike rates (spikes per second) with each row corresponding to a unit
    %   - spikeRate_t: Time array corresponding to each spike rate value

    % Assign parameters if needed
    if nargin < 7
        fs=20000; % record sampling frequency. defualt for my recs. 
    end
    
    if nargin < 6
        OL = 100; % ms
    end
    if nargin < 5
        meanWin = 1000; % ms
    end
    if nargin < 4
        win = 2 * 60 * 60 * 1000; % 2 hours in ms
    end
    if nargin < 3
        tStart = 0;
    end
    if nargin < 2 
        units = 1; % Default to unit 1 if no units are specified
    end

    % Extract spike data
    spikeClusters = data.spike_clusters;
    spikeTimes_ms = data.spike_times/(fs/1000); % Convert to ms if needed (times are in freq samples)/ 

    % Define the end time of the analysis window
    tEnd = tStart + win;

    % Generate time windows
    spikeRate_t = tStart:OL:(tEnd - meanWin);
    numWindows = length(spikeRate_t);
    numUnits = length(units);

    % Pre-allocate spike rate matrix
    spikeRate = zeros(numUnits, numWindows);

    % Loop over each unit
    for u = 1:numUnits
        % Filter spikes for the current unit
        unitSpikes = spikeTimes_ms(spikeClusters == units(u));

        % Calculate spike rate for each window for the current unit
        for i = 1:numWindows
            % Define the start and end of the current window
            winStart = spikeRate_t(i);
            winEnd = winStart + meanWin;

            % Count spikes within the current window
            spikesInWindow = unitSpikes(unitSpikes >= winStart & unitSpikes < winEnd);
            spikeCount = length(spikesInWindow);

            % Calculate the spike rate (spikes per second)
            spikeRate(u, i) = spikeCount / (meanWin / 1000); % Convert ms to s
        end
    end
end

