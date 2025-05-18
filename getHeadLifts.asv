
function [angleF, angleF_t] = getHeadLifts(angles,angle_t,newFS,lowpasspass)

    % initiating parameters:
    if nargin <2
        samplesPerS = newFS; %samples per s - hz
    else
        samplesPerS = 100;
    end
    if nargin ==3
        lowpasspass =5;
    end
    
    duration_s = (angle_t(end)-angle_t(1)) / 1000;
    numSamples = round(duration_s * samplesPerS);
    angle_tL = linspace(angle_t(1),angle_t(end),numSamples);

    [x,uniqueIndx] = unique(angle_t);
    v = angles(uniqueIndx);
    anglesL = interp1(x,v,angle_tL);

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
    anglesLF = filtfilt(d, anglesL);

    % Convert filtered data from radians to degrees
    angleF =  anglesLF* (180 / pi);
    angleF_t = angle_tL;

end    


% % claculate the degree, low pass and plot.
