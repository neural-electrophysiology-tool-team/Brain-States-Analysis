
%% CREATE SIMPLIFIED EXCEL FILE
%
% Output columns:
%   Animal
%   System
%   Folder
%
% Filters:
%   1. Select specific animals
%   2. Include/exclude recNames
%   3. Exclude rows where Exlude == 1
%   4. Optionally restrict recordings to unique recordings
%      in reward_sound and/or neutral_sound tables.
%
% Recording identity = unique combination of Animal + Record.

clear;
clc;


%% =========================================================
%  FILES AND PATHS
%  =========================================================

inputFile = '/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx';

dataFile = '/home/ofir/Documents/MATLAB/data4paper.mat';

outputFolder = '/home/ofir/Documents/MATLAB';

outputFile = fullfile(outputFolder, ...
    'brainStatesWakeN_simplified.xlsx');

% Create destination folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% =========================================================
%  OPTIONS TO EDIT
%  =========================================================

%% 1. ANIMAL FILTER

% [] = include all animals
% Example: ["PV159", "PV167", "PV208"]

includeAnimals = ["PV106", "PV126", "PV143","PV153", "PV157", "PV159","PV161", "PV162"];

%% 2. recNames FILTER

% Include records whose recNames contain at least one term.
% [] = no inclusion restriction.

includeRecNames = ["Hunter","MT","Mt"];

% Exclude records whose recNames contain at least one term.
% [] = no exclusion.

excludeRecNames = "RMHunter";

%% 3. REWARD / NEUTRAL RECORDING FILTER

% true  = include recordings appearing in reward_sound
% false = don't use reward_sound

useRewardSound = true;

% true  = include recordings appearing in neutral_sound
% false = don't use neutral_sound

useNeutralSound = true;

% If both are true, the script uses the UNION of the recordings
% appearing in the two tables.
%
% If both are false, no recording restriction is applied.

%% =========================================================
%  READ EXCEL
%  =========================================================

T = readtable(inputFile, ...
    'Sheet', 1, ...
    'VariableNamingRule', 'preserve');

varNames = string(T.Properties.VariableNames);

% Find column indices (case-insensitive)

animalCol = find(strcmpi(varNames, 'Animal'), 1);

systemCol = find(strcmpi(varNames, 'System'), 1);

folderCol = find(strcmpi(varNames, 'folder'), 1);

recNamesCol = find(strcmpi(varNames, 'recNames'), 1);

recordCol = find(strcmpi(varNames, 'Record'), 1);

excludeCol = find(strcmpi(varNames, 'Exclude'), 1);

% Check required columns

if isempty(animalCol) || isempty(systemCol) || isempty(folderCol)
    error('Missing required Excel columns: Animal, System, or folder.');
end

if isempty(excludeCol)
    error('Could not find the Exlude column in the Excel file.');
end

if (useRewardSound || useNeutralSound) && isempty(recNamesCol)
    error('Recording filtering requires a Record column in the Excel file.');
end

if (~isempty(includeRecNames) || ~isempty(excludeRecNames)) ...
        && isempty(recNamesCol)
    error('recNames filtering requested, but recNames column was not found.');
end

fprintf('Original Excel rows: %d\n', height(T));

%% =========================================================
%  FILTER 1: EXCLUDE Exlude == 1
%  =========================================================

excludeValues = T{:, excludeCol};

if iscell(excludeValues) || isstring(excludeValues) || iscategorical(excludeValues)
    excludeValues = str2double(string(excludeValues));
end

keep = ~isnan(double(excludeValues)) | double(excludeValues) ~= 1;
% The expression above is equivalent to keeping rows with Exlude ~= 1,
% including missing values.

keep = double(excludeValues) ~= 1 | isnan(double(excludeValues));

fprintf('Rows with Exlude == 1: %d\n', sum(~keep));

%% =========================================================
%  FILTER 2: ANIMAL
%  =========================================================

if ~isempty(includeAnimals)

    animalValues = string(T{:, animalCol});

    animalOptions = string(includeAnimals);

    animalMatch = ismember( ...
        lower(animalValues), ...
        lower(animalOptions));

    keep = keep & animalMatch;

end

%% =========================================================
%  FILTER 3: recNames
%  =========================================================

if ~isempty(includeRecNames) || ~isempty(excludeRecNames)

    recNamesValues = string(T{:, recNamesCol});

    recNamesValues(ismissing(recNamesValues)) = "";

    % Inclusion condition

    if ~isempty(includeRecNames)

        includeTerms = string(includeRecNames);

        includeMatch = false(height(T), 1);

        for k = 1:numel(includeTerms)

            includeMatch = includeMatch | ...
                contains(recNamesValues, includeTerms(k), ...
                'IgnoreCase', true);

        end

        keep = keep & includeMatch;

    end

    % Exclusion condition

    if ~isempty(excludeRecNames)

        excludeTerms = string(excludeRecNames);

        excludeMatch = false(height(T), 1);

        for k = 1:numel(excludeTerms)

            excludeMatch = excludeMatch | ...
                contains(recNamesValues, excludeTerms(k), ...
                'IgnoreCase', true);

        end

        keep = keep & ~excludeMatch;

    end

end

%% =========================================================
%  FILTER 4: REWARD / NEUTRAL RECORDINGS
%  =========================================================

if useRewardSound || useNeutralSound

    fprintf('\nLoading reward_sound and neutral_sound...\n');

    % Load only the required tables
    D = load(dataFile, 'reward_sound', 'neutral_sound');

    % Initialize recording list
    uniqueRecordings = table( ...
        strings(0,1), strings(0,1), ...
        'VariableNames', {'Animal', 'Record'});

    %% REWARD SOUND

    if useRewardSound

        if ~isfield(D, 'reward_sound')
            error('reward_sound was not found in data4paper.mat.');
        end

        reward_sound = D.reward_sound;

        requiredCols = {'Animal', 'Record'};

        if ~all(ismember(requiredCols, reward_sound.Properties.VariableNames))
            error('reward_sound must contain Animal and Record columns.');
        end

        rewardRecordings = table( ...
            string(reward_sound.Animal), ...
            string(reward_sound.Record), ...
            'VariableNames', {'Animal', 'Record'});

        uniqueRecordings = [uniqueRecordings; rewardRecordings];

        fprintf('Reward sound recordings: %d\n', ...
            height(unique(rewardRecordings)));

    end

    %% NEUTRAL SOUND

    if useNeutralSound

        if ~isfield(D, 'neutral_sound')
            error('neutral_sound was not found in data4paper.mat.');
        end

        neutral_sound = D.neutral_sound;

        requiredCols = {'Animal', 'Record'};

        if ~all(ismember(requiredCols, neutral_sound.Properties.VariableNames))
            error('neutral_sound must contain Animal and Record columns.');
        end

        neutralRecordings = table( ...
            string(neutral_sound.Animal), ...
            string(neutral_sound.Record), ...
            'VariableNames', {'Animal', 'Record'});

        uniqueRecordings = [uniqueRecordings; neutralRecordings];

        fprintf('Neutral sound recordings: %d\n', ...
            height(unique(neutralRecordings)));

    end

    %% FIND UNIQUE ANIMAL + RECORD COMBINATIONS

    % Remove missing values
    validRows = ~ismissing(uniqueRecordings.Animal) & ...
                ~ismissing(uniqueRecordings.Record);

    uniqueRecordings = uniqueRecordings(validRows, :);

    % Unique combinations of Animal and Record
    uniqueRecordings = unique(uniqueRecordings, 'rows');

    fprintf('Total unique recordings (reward + neutral): %d\n', ...
        height(uniqueRecordings));

    %% MATCH EXCEL RECORDINGS

    excelAnimals = lower(strtrim(string(T{:, animalCol})));
    excelRecords = lower(strtrim(string(T{:, recordCol})));

    targetAnimals = lower(strtrim(uniqueRecordings.Animal));
    targetRecords = lower(strtrim(uniqueRecordings.Record));

    % Create matching keys
    excelKeys = excelAnimals + "|" + excelRecords;

    targetKeys = targetAnimals + "|" + targetRecords;

    recordingMatch = ismember(excelKeys, targetKeys);

    keep = keep & recordingMatch;

end

%% =========================================================
%  APPLY ALL FILTERS
%  =========================================================

T = T(keep, :);

fprintf('\n========================================\n');
fprintf('Filtering completed.\n');
fprintf('Remaining Excel rows: %d\n', height(T));
fprintf('========================================\n');

%% =========================================================
%  EXTRACT OUTPUT COLUMNS
%  =========================================================

Animal = T{:, animalCol};

System = T{:, systemCol};

Folder = string(T{:, folderCol});

%% =========================================================
%  REMOVE RECORD NODE FROM FOLDER PATHS
%  =========================================================

% Remove trailing "/Record Node 101"
% Also supports Windows-style backslashes.

Folder = regexprep(Folder, ...
    '[\\/]Record Node 101[\\/]*$', '');

%% =========================================================
%  CREATE OUTPUT TABLE
%  =========================================================

Tout = table(Animal, System, Folder, ...
    'VariableNames', {'Animal', 'System', 'Folder'});

%% =========================================================
%  SAVE EXCEL
%  =========================================================

writetable(Tout, outputFile, 'Sheet', 'Data');

fprintf('\n========================================\n');
fprintf('Finished successfully!\n');
fprintf('Output file: %s\n', outputFile);
fprintf('Total exported rows: %d\n', height(Tout));
fprintf('========================================\n');
