
function [angleF, angleF_t] = getHeadLifts(angles,angle_t,newFS,lowpasspass)

    % Input data
    anglesX = angles(1,:); % Head angle X-axis
    anglesY = angles(2,:); % Head angle Y-axis
    anglesZ = angles(3,:); % Head angle Z-axis
    % timeParts = parts; % Define start and end times for each part in a cell. number of parts ins number of cells., in m
  
    %interpulate
    if nargin >2
        samplesPerS = newFS; %samples per s - hz
    else
        samplesPerS = 100;
    end
    duration_s = (angle_t(end)-angle_t(1)) / 1000;
    numSamples = round(duration_s * samplesPerS);
    angle_tL = linspace(angle_t(1),angle_t(end),numSamples);

    [x,uniqueIndx] = unique(angle_t);
    v = anglesZ(uniqueIndx);
    anglesZL = interp1(x,v,angle_tL);

    % Parameters - low pass filter
    if nargin>3
         lowPassFreq = lowpasspass;
    else
    lowPassFreq = 3; % Low-pass cutoff frequency in Hz
    end
    filterOrder = 3; % Filter order

    % Design the filter
    d = designfilt('lowpassiir', 'FilterOrder', filterOrder, ...
        'HalfPowerFrequency', lowPassFreq, 'SampleRate', samplesPerS);

    % Apply the filter
    anglesZLF = filtfilt(d, anglesZL);

    % Convert filtered data from radians to degrees
    angleF =  anglesZLF* (180 / pi);
    angleF_t = angle_tL;

end    


% % claculate the degree, low pass and plot.
