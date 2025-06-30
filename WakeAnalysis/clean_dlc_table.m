function [behavtable] = clean_dlc_table(behavtable)
% Milan's function
%% replaces uncertain values with linear fitting
certain_cutOff = 0.99; % change this to accept more/less values
bodyparts = behavtable.Properties.VariableNames;
for i = 2:length(bodyparts) % from 2 because first is the description "bodyparts"
    body_part = bodyparts{i};
    if contains(body_part, '_2') % doesnt correct the precision values of the labeling
        continue
    end
    % add _2
    if contains(body_part, '_1')
        if contains(body_part, 'Pupil_1') && sum(body_part=='_') < 2 % because pupil_12 and pupil_1 have a _1 in the name
            body_prec = append(body_part, '_2');
        else
            body_prec = append(body_part(1:end-2), '_2'); % deletes the _1 if it exists
        end
    else
        body_prec = append(body_part, '_2');
    end
    body_pos = behavtable.(body_part);
    body_pos(behavtable.(body_prec)<certain_cutOff) = NaN; % deletes uncertain values to make the data clearer for plotting
    %% this is to make values NaN if the bodypart wasnt seen in very long (e.g. eye camera fell off)
    nanLocations = isnan(body_pos); % Get logical array of whether element is NaN or not.
    props = regionprops(nanLocations, 'Area', 'PixelIdxList');
    propss = struct2cell(props)';
    indx = find([propss{:, 1}]>10000); % finds periods with more than 10000 NaN values in a row
    indx(indx==1) = []; %does not delete the first period without labeling
    for i=1:length(indx)
        indxx(i) = propss{indx(i),2}(1);
    end
    %%
    if isnan(body_pos(1)) && any(~isnan(body_pos)) % if the the bodypart is not labeled in the start of the recording, assume that it does not move until it is labeled
        body_pos(1)=body_pos(find(~isnan(body_pos),1));
    end
    body_pos = fillmissing(body_pos,'linear'); %interpolates the uncertain values, to not have gaps when you plot the data
    for i=1:length(indx)
        body_pos(indxx(i):indxx(i)+propss{indx(i),1}-1) = NaN; % makes the long NaN periods NaNs again after the fillmissing
    end
    behavtable.(body_part) = body_pos; %writes the new values to the table
end
end