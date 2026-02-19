%% running the DLC model: preparations of the data:

WA = wakeAnalysis('/media/sil1/Data/Nitzan/Experiments/brainStatesWake.xlsx');
analysisFolderWake = '/media/sil1/Data/Nitzan/WakeStatesPaper/plots';
%% getting or creating the calibration table
%% all wake trials 
wakeTrials = contains(WA.recTable.recNames,'Wake') & ...
             cellfun(@isempty,WA.recTable.Remarks);
wakeTable = WA.recTable(wakeTrials,:);

%% subset to analyse for the effect:
% rng('shuffle')   % optional: for randomness
% 
% animals = unique(wakeTable.Animal);
% rowsToKeep = [];
% 
% for i = 1:numel(animals)
%     idx = find(strcmp(wakeTable.Animal,animals{i}));   % rows of this animal
%     nPick = min(3, numel(idx));   % in case some animals have < 3 rows
%     rowsToKeep = [rowsToKeep; idx(randperm(numel(idx), nPick))];
% end
% 
% wakeSubset = wakeTable(rowsToKeep, :);
%% for each unique animal - the reptilearns optional:
% % animals: cell array or string array of animals you care about
% subT = WA.recTable(ismember(WA.recTable.Animal, animals), {'Animal','system'});

% % keep unique Animal-system pairs
% caliTable = unique(subT, 'rows');
load(["Brain-States-Analysis/WakeAnalysis/calibrationTable.mat"])
% load([analysisFolderWake filesep "calibrationTable.mat"])
wakeSubset = readtable([analysisFolderWake filesep 'wakeSubset.xlsx']);
% save([analysisFolderWake filesep 'wakeSubset.xlsx'],"wakeSubset")
%%
A = wakeSubset(:, {'Animal','recNames'});
B = WA.recTable(:, {'Animal','recNames'});

% Find matching rows
[isMatch, idx] = ismember(A, B, 'rows');
wakeSubsetAll = WA.recTable(idx,:);
% writetable(wakeSubsetAll,[analysisFolderWake filesep 'wakeSubsetAll.xlsx'])
%% Build job list for python loop: one row per video folder with its calib_dir
% Inputs:
%   wakeSubset : table with variables {'Animal','system','folder'}
%   caliTable  : table with variables {'Animal','system','caliPath'}
%
% Output:
%   jobsTable  : table with variables {'Animal','system','folder','calib_dir'}
%   plus a CSV saved to analysisFolderWake

% % --- 0) Make sure types match (string is easiest) ---
% wakeSubset.Animal = string(wakeSubset.Animal);
% wakeSubset.system = string(wakeSubset.system);
% wakeSubset.folder = string(wakeSubset.folder);
% 
% caliTable.Animal  = string(caliTable.Animal);
% caliTable.system  = string(caliTable.system);
% caliTable.caliPath = string(caliTable.caliPath);

% --- 1) Keep only the needed columns (avoid accidental extra vars) ---
leftT  = wakeSubsetAll(:, {'Animal','system','folder'});
rightT = caliTable(:, {'Animal','system','caliPath'});

% --- 2) Join on Animal+system to attach caliPath to each row in wakeSubset ---
jobsTable = outerjoin(leftT, rightT, ...
    'Keys', {'Animal','system'}, ...
    'MergeKeys', true, ...
    'Type', 'left');   % keep all wakeSubset rows

jobsTable.folder = string(fileparts(jobsTable.folder));

% Rename to what python expects
jobsTable.Properties.VariableNames{'caliPath'} = 'calib_dir';

jobsTable.folder    = string(jobsTable.folder);
jobsTable.calib_dir = string(jobsTable.calib_dir);
% add folder paths for hpc: 
jobsTable.folder_hpc    = arrayfun(@fixHPCpath, jobsTable.folder);
jobsTable.calib_dir_hpc = arrayfun(@fixHPCpath, jobsTable.calib_dir);

% --- 5) Save as CSV for python ---
outCsv = fullfile(analysisFolderWake, 'deeplabcut_jobs.csv');
writetable(jobsTable, outCsv);

% --- 5.5) Save as CSV for python in sil3 ---
outCsv = fullfile("/media/sil3/Data/Nitzan/hpc_dlc/deeplabcut_jobs_hpc.csv");
writetable(jobsTable, outCsv);

% --- 6) Also save the MATLAB table (optional) ---
filename = fullfile(analysisFolderWake, 'deeplabcut_jobs.mat');
save(filename, 'jobsTable');

fprintf("Wrote %d jobs to %s\n", height(jobsTableClean), outCsv);
fprintf("saved in: %s\n",outCsv)



%% check prediction file:
%first clean:
%%

function out = fixHPCpath(in)

    parts = split(in, filesep);

    % Remove 3rd folder
    parts(4) = [];

    % Rebuild path
    out = "/a/home/cc/lifesci/nitzanalbeck" + strjoin(parts, filesep);

end
