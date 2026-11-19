%% ============================================================
%  ONE SVM ON ALL FREQUENCY BANDS - SEPARATELY FOR EACH ANIMAL
%  + OPTIONAL ABLATION (5 extra SVMs, each with ONE band removed)
%  averaged over random train/test splits, with an optional permutation null
%
%  Analysis modes (analysisModes):
%   BUG = hunter tables (hunterXXX), split into 'circle' / 'other' trials (no-strike trials only)
%   pre_vs_onset          : bug PRE window vs bug ONSET window                                (PAIRED)
%   <mode>_strike         : for pre_vs_onset, pre_vs_offset, onset_vs_offset and circle_vs_other the SAME analysis on
%                           STRIKE trials only (instead of the no-strike trials used by the plain modes), e.g.
%                           'pre_vs_onset_strike'. Strike trials = Strike_ms is not NaN and
%                           Strike_ms > strikeMinTime_s*1000; Triger_ms is not used.
%   pre_vs_offset         : bug PRE window vs bug OFFSET window                               (PAIRED)
%   onset_vs_offset       : bug ONSET window (class 0) vs bug OFFSET window (class 1)         (PAIRED)
%   circle_vs_other       : circle trials vs other trials, onset windows                      (UNPAIRED)
%   bug_vs_meadow_onset   : bug onset    vs meadow onset   (separately for circle / other bugs)  (UNPAIRED)
%   bug_vs_meadow_offset  : bug offset   vs meadow offset  (separately for circle / other bugs)  (UNPAIRED)
%   meadow_move_vs_stat   : meadow STAT onset (class 0) vs meadow MOVE onset (class 1), same row (PAIRED)
%   bug_vs_sound_onset    : bug onset    vs sound onset    (separately for circle / other bugs)  (UNPAIRED)
%   bug_vs_sound_offset   : bug offset   vs sound offset   (separately for circle / other bugs)  (UNPAIRED)
%   category_pairs        : EVERY couple of categories that appear together in one figure (see pairDefs):
%        PAIRED   (two windows of the same trials): circle, circle strike, other, other strike, meadow, sound
%        UNPAIRED (different events): onset (circle/other/meadow/sound), offset (circle/other/sound),
%                 onset strike, offset strike (circle vs other), and
%                 strike vs no-strike separately for circle / other and separately for ONSET / OFFSET.
%        Paired comparisons: trials whose two windows overlap in time are removed from THAT comparison only.
%        Each comparison = its own SVM (+ permutation null if runPermutations).
%
%  Time windows (all windows are 5 s = 5 bins of 1 s, the end point is NOT included):
%   PRE            : [-6, -1) s from stimulus ONSET
%   ONSET window   : [-1,  4) s from stimulus ONSET
%   OFFSET window  : [ 1,  6) s from stimulus OFFSET
%   PRE-OFFSET     : the PRE window taken relative to the OFFSET (category_pairs; types in preOffsetTypes)
%
%   PAIRED   : row i of class 0 and row i of class 1 are the same trial; a pair stays together in train or test.
%   UNPAIRED : the two groups are different events. Per animal, the SAME number of events is drawn at random
%              from each group (n = size of the smaller group after removing events with missing data).
%
%  Preprocessing: band-averaged RAW power -> (optional) each band divided by the SUM OF THE 5 BANDS' mean powers
%  at every time point -> (optional) z-scoring, set by zMode:
%     'recording' : z-score each band over the WHOLE recording
%     'window'    : z-score each band inside the event windows, before binning
%                     PAIRED modes -> both windows of the SAME trial z-scored TOGETHER
%                     UNPAIRED modes -> every event window z-scored by itself
%     'trial'     : ALL windows of the same trial (hunter: pre, onset, [pre-offset], offset; meadow row: pre-move,
%                   onset-move, pre-stat, onset-stat; sound: pre, onset, offset) are z-scored TOGETHER (one
%                   mean/std per band over the union of their time points, overlapping samples counted once).
%                   Identical normalization in every comparison. A trial with a window without data is skipped
%                   everywhere. Events without trial info (meadow offsets) use their own window.
%     'none'      : no z-scoring
%
%  Model 1 = SVM on all bands together (features = nBands x nBins, band-fastest).
%  If runAblation = true, models 2..6 = the same SVM with ONE band removed.
%  All models use IDENTICAL train/test splits (rng(1) is reset before each fit).
%  Weights are multiplied by the SVM's accuracy (weightScaleFcn).
%  Permutation null is OPTIONAL (runPermutations), computed for the all-bands model only.
%
%  Figures: per group of comparisons one accuracy figure (one panel per comparison: star per animal + mean,
%  boxplot of the permutations, p-value). Weight figures (plotWeights) and the summary figure
%  (plotSummaryFig) are optional.
%
%  Saved: svm_allbands_rawnorm_<mode>.mat with batchResults (acc, pVal, nullAcc, beta, betaScaled, nPerClass,
%  binCenters, xLab, design, class0, class1 for every animal / comparison / model).
% =============================================================

clear; clc;

%% ------------------------------------------------------------
%  User Parameters
% ------------------------------------------------------------
analysisModes = {'category_pairs'};   % can be combined with the old modes, e.g. {'category_pairs','pre_vs_onset'}
% analysisModes = {'pre_vs_onset', 'pre_vs_offset', 'onset_vs_offset', 'circle_vs_other', ...
%                  'pre_vs_onset_strike', 'pre_vs_offset_strike', 'onset_vs_offset_strike', 'circle_vs_other_strike'};
% analysisModes = {'bug_vs_meadow_onset', 'bug_vs_meadow_offset', 'meadow_move_vs_stat'};
% analysisModes = {'bug_vs_sound_onset', 'bug_vs_sound_offset'};

hunterTables = {'hunter106', 'hunter126', 'hunter143', 'hunter153', 'hunter157', 'hunter159', 'hunter161','hunter162', 'hunter167', 'hunter208', 'hunter213', 'hunter228'};

minTrialsRequired = 15;
strikeMinTime_s   = 4;   % strike trials: keep only trials with Strike_ms > this many s
minTrialsStrike   = 15;  % minimum number of valid strike trials; lower it if there are few strike trials
binSize           = 1;   % s
winLen            = 5;   % s (length of every window; the window end point is NOT included)
preStart          = -6;  % PRE window    = [preStart,    preStart + winLen)    s from stim ONSET   -> [-6, -1)
onsetStart        = -1;  % ONSET window  = [onsetStart,  onsetStart + winLen)  s from stim ONSET   -> [-1,  4)
offsetStart       = 1;   % OFFSET window = [offsetStart, offsetStart + winLen) s from stim OFFSET  -> [ 1,  6)
trainFraction     = 0.8;
nSplits           = 1;     % random 80/20 splits per dataset (real or permuted)
nPerm             = 10;    % label permutations (used only if runPermutations = true)
runPermutations   = true;  % false = skip the permutation null
runAblation       = false; % true = also run 5 models, each with ONE band removed (4 of 5 bands)
plotWeights       = false; % weight figures (false: only accuracy + permutation boxplots)
plotSummaryFig    = false; % summary figure across analyses (pre_vs_onset / pre_vs_offset / onset_vs_offset modes)

preOffsetTypes = {'circle'};   % trial types that have a pre-offset window; {'circle','other'} adds it to 'other' too

zMode          = 'trial';   % 'recording' | 'window' | 'trial' | 'none'  (see header)
normToSum      = false;     % BEFORE z-score: divide each band's raw mean power by the sum of the 5 bands' powers
weightScaleFcn = @(a) a;    % weight * accuracy  (or @(a) max(a-0.5,0))
barUseAbs      = true;      % bar plots: false = signed weights, true = |weights|
animalColor    = [0.20 0.45 0.80];   % color of the weight plots / the all-bands model

% colors of the 6 models in the accuracy figure (row 1 = All bands, then w/o Delta, Theta+Alpha, Beta, GammaLow, GammaHigh)
modelColors = [0.20 0.45 0.80;    % All bands   - blue
               0.90 0.45 0.10;    % w/o Delta       - orange
               0.20 0.65 0.30;    % w/o Theta+Alpha - green
               0.60 0.30 0.70;    % w/o Beta        - purple
               0.85 0.15 0.15;    % w/o GammaLow    - red
               0.55 0.35 0.15];   % w/o GammaHigh   - brown

% ---- meadow / sound data (EDIT these names if your variables differ) ----
extraMatFiles  = {};                  % optional extra .mat files merged into loadData
meadowVarFmt   = 'reg_meadow%d';      % meadow table of animal 106 -> loadData.reg_meadow106
                                      %   needed columns: Record, Triger_ms, Move_ms (Animal optional)
meadowState    = 'move';              % 'stat' or 'move': which meadow state is "meadow onset/offset"
                                      %   in bug_vs_meadow_onset / bug_vs_meadow_offset
soundVarName   = 'reward_sound';      % sound table (all animals); needed columns: Animal, Record, Triger_ms
soundDur_ms    = 1000;                % sound offset = onset + soundDur_ms
soundExclWin_s = 5;                   % exclude sounds with another sound onset within +- this many s
soundRecordFilter = @(rec) contains(lower(string(rec)), 'wake') & ~contains(lower(string(rec)), 'hunter');
                                      % which sound records are used (wake only, no hunter records)
minTrialsSound = 2;                   % sound: no 15-trial minimum; only a technical floor for an 80/20 split

bandRanges = [0 4; 5 12; 13 30; 31 40; 60 80];
bandNames  = ["Delta (0-4 Hz)", "Theta+Alpha (5-12 Hz)", "Beta (13-30 Hz)", ...
              "GammaLow (31-40 Hz)", "GammaHigh (60-80 Hz)"];
bandShort  = ["Delta", "Theta+Alpha", "Beta", "GammaLow", "GammaHigh"];   % used in ablation labels
nBands   = size(bandRanges, 1);
nBinsWin = round(winLen / binSize);
nPermEff = nPerm * runPermutations;   % 0 -> no permutations are run

% window definitions shared by all trial sets (used by zMode 'trial' and category_pairs)
T = struct('winLen', winLen, 'preStart', preStart, 'onsetStart', onsetStart, 'offsetStart', offsetStart, ...
           'hunterStarts', [preStart onsetStart preStart offsetStart], ...   % windows: 1 pre, 2 onset, 3 pre-offset, 4 offset
           'meadowStarts', [preStart onsetStart preStart onsetStart], ...    % 1 pre-move, 2 onset-move, 3 pre-stat, 4 onset-stat
           'soundStarts',  [preStart onsetStart offsetStart], ...            % 1 pre, 2 onset, 3 offset
           'preOffsetTypes', {preOffsetTypes}, 'strikeMin_s', strikeMinTime_s, ...
           'soundDur_ms', soundDur_ms, 'soundRecordFilter', soundRecordFilter);

loadData = load('/home/ofir/Documents/MATLAB/data4paper.mat');
for k = 1:numel(extraMatFiles)
    ex = load(extraMatFiles{k});
    for fn = fieldnames(ex)', loadData.(fn{1}) = ex.(fn{1}); end
end


% Caches (speed): each recording is loaded + processed ONCE for the whole script,
% wakeAnalysis is built once per excel file.
recCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
saCache  = containers.Map('KeyType', 'char', 'ValueType', 'any');

% accuracy of the all-bands SVM of every (mode, trial type, animal), for the summary figure
summaryAll = struct('Base', {}, 'Strike', {}, 'TrialType', {}, 'Animal', {}, 'acc', {});

savePPT = true;                         % put all figures into one pptx (one slide per figure, title = figure name)
pptFile = 'svm_figures.pptx';

%% ============================================================
%  LOOP OVER ANALYSIS MODES
% =============================================================
for mIdx = 1:numel(analysisModes)
    analysisMode = analysisModes{mIdx};
    fprintf('\n===== MODE: %s =====\n', analysisMode);
    isStrike = endsWith(analysisMode, '_strike');          % strike-trial variant of a bug analysis
    baseMode = erase(analysisMode, '_strike');
    isPairs  = strcmp(baseMode, 'category_pairs');

    [needBug, needMeadow, needSound, splitBug] = modeNeeds(analysisMode);
    groups = [];
    if isPairs
        pdefs  = pairDefs(preOffsetTypes);
        types  = {pdefs.name};
        groups = pairGroups(pdefs);
    elseif splitBug, types = {'circle', 'other'}; else, types = {analysisMode}; end
    soundTbl = [];
    if needSound
        if isfield(loadData, soundVarName), soundTbl = loadData.(soundVarName);
        elseif ~isPairs, error('Sound table "%s" not found in loadData (see soundVarName / extraMatFiles).', soundVarName);
        end
    end

    batchResults = struct('Animal', {}, 'RunLabel', {}, 'TrialType', {}, 'ModelIdx', {}, 'ModelLabel', {}, ...
                          'acc', {}, 'beta', {}, 'betaScaled', {}, 'nullAcc', {}, 'nPerClass', {}, 'pVal', {}, ...
                          'binCenters', {}, 'xLab', {}, 'design', {}, 'class0', {}, 'class1', {});
    runCounter = 0;
    binCenters = []; xLab = '';

    % ---- each animal is analysed completely separately ----
    for tIdx = 1:numel(hunterTables)
        hTableName = hunterTables{tIdx};
        animalID = str2double(regexp(hTableName, '\d+', 'match', 'once'));

        currentHunter = [];
        if needBug
            if ~isfield(loadData, hTableName), continue; end
            currentHunter = loadData.(hTableName);
        end
        meadowTbl = [];
        if needMeadow
            mName = sprintf(meadowVarFmt, animalID);
            if isfield(loadData, mName)
                meadowTbl = loadData.(mName);
            elseif ~isPairs
                fprintf('Skipping %s (%s): no meadow table %s\n', hTableName, analysisMode, mName);
                continue;
            end
        end

        if ismember(animalID, [208, 213, 228, 167])
            excelPath = '/media/sil3/Data/Ofir/BrainStatesExperimentsOfir.xlsx';
        else
            excelPath = '/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx';
        end
        if ~isKey(saCache, excelPath), saCache(excelPath) = wakeAnalysis(excelPath); end
        SA = saCache(excelPath);

        X = struct('SA', SA, 'bandRanges', bandRanges, ...
                   'binSize', binSize, 'recCache', recCache, 'normToSum', normToSum, 'zMode', zMode);
        if isPairs
            trialSets = buildTrialSets(currentHunter, meadowTbl, soundTbl, animalID, T);
        end

        for cType = types
            cType = cType{1};
            desc = sprintf('%s | %s | %s', hTableName, analysisMode, cType);

            minN   = minTrialsRequired;
            if needSound, minN = minTrialsSound; end
            if isStrike, minN = minTrialsStrike; end
            nBins  = nBinsWin;
            class0Name = '';  class1Name = '';

            %% ---- build the two groups of events (E0 = class 0, E1 = class 1) ----
            switch baseMode
                case 'pre_vs_onset'
                    E0 = hunterEvents(bugTrials(currentHunter, cType, isStrike, strikeMinTime_s), 'onset', T);
                    E1 = E0;
                    start0 = preStart; start1 = onsetStart; design = 'paired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim ONSET; pre window = [%g, %g) s, onset window = [%g, %g) s)', ...
                                   preStart, preStart + winLen, onsetStart, onsetStart + winLen);

                case 'pre_vs_offset'
                    tblB = bugTrials(currentHunter, cType, isStrike, strikeMinTime_s);
                    E0 = hunterEvents(tblB, 'onset', T);
                    E1 = hunterEvents(tblB, 'offset', T);
                    start0 = preStart; start1 = offsetStart; design = 'paired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim OFFSET; offset window = [%g, %g) s; pre window = [%g, %g) s from ONSET)', ...
                                   offsetStart, offsetStart + winLen, preStart, preStart + winLen);

                case 'onset_vs_offset'
                    tblB = bugTrials(currentHunter, cType, isStrike, strikeMinTime_s);
                    E0 = hunterEvents(tblB, 'onset', T);
                    E1 = hunterEvents(tblB, 'offset', T);
                    start0 = onsetStart; start1 = offsetStart; design = 'paired'; binStart = start0;
                    xLab = sprintf('Time bin (axis = s from stim ONSET [%g, %g); class 1 = same bin of [%g, %g) s from stim OFFSET)', ...
                                   onsetStart, onsetStart + winLen, offsetStart, offsetStart + winLen);

                case 'circle_vs_other'
                    E0 = hunterEvents(bugTrials(currentHunter, 'circle', isStrike, strikeMinTime_s), 'onset', T);
                    E1 = hunterEvents(bugTrials(currentHunter, 'other',  isStrike, strikeMinTime_s), 'onset', T);
                    start0 = onsetStart; start1 = onsetStart; design = 'unpaired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim ONSET; window = [%g, %g) s)', onsetStart, onsetStart + winLen);

                case {'bug_vs_meadow_onset', 'bug_vs_meadow_offset'}
                    align = 'onset'; wStart = onsetStart;
                    if endsWith(analysisMode, 'offset'), align = 'offset'; wStart = offsetStart; end
                    E0 = hunterEvents(bugSubset(currentHunter, cType), align, T);
                    E1 = meadowEvents(meadowTbl, animalID, [meadowState ternary(strcmp(align,'onset'),'On','Off')], T);
                    start0 = wStart; start1 = wStart; design = 'unpaired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim %s; meadow %s; window = [%g, %g) s)', upper(align), meadowState, wStart, wStart + winLen);

                case 'meadow_move_vs_stat'
                    E0 = meadowEvents(meadowTbl, animalID, 'statOn', T);
                    E1 = meadowEvents(meadowTbl, animalID, 'moveOn', T);
                    start0 = onsetStart; start1 = onsetStart; design = 'paired'; binStart = start1;
                    xLab = sprintf('Time bin (s from meadow ONSET, window = [%g, %g) s; class 0 = stat, class 1 = move)', onsetStart, onsetStart + winLen);

                case {'bug_vs_sound_onset', 'bug_vs_sound_offset'}
                    align = 'onset'; wStart = onsetStart;
                    if endsWith(analysisMode, 'offset'), align = 'offset'; wStart = offsetStart; end
                    E0 = hunterEvents(bugSubset(currentHunter, cType), align, T);
                    E1 = soundEvents(soundTbl, animalID, soundRecordFilter, align, soundDur_ms, soundExclWin_s, T);
                    start0 = wStart; start1 = wStart; design = 'unpaired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim %s; window = [%g, %g) s)', upper(align), wStart, wStart + winLen);

                case 'category_pairs'
                    d = pdefs(strcmp({pdefs.name}, cType));
                    [E0, E1, start0, start1, design, nOverlap] = buildPair(d, trialSets, T);
                    mt = struct('circle', minTrialsRequired, 'other', minTrialsRequired, 'circleS', minTrialsStrike, ...
                                'otherS', minTrialsStrike, 'meadow', minTrialsRequired, 'sound', minTrialsSound);
                    minN = min(mt.(d.setA), mt.(d.setB));
                    binStart = start1;
                    xLab = sprintf('Time bin (s from window start; class 0 = %s, class 1 = %s)', d.labelA, d.labelB);
                    class0Name = sprintf('%s: %s', d.setA, d.labelA);
                    class1Name = sprintf('%s: %s', d.setB, d.labelB);
                    if nOverlap > 0
                        fprintf('   %s: %d trials removed (the two windows overlap)\n', desc, nOverlap);
                    end

                otherwise
                    error('Unknown analysisMode: %s', analysisMode);
            end

            if isStrike, xLab = ['STRIKE trials: ' xLab]; end

            nE0 = height(E0); nE1 = height(E1);
            if nE0 < minN || nE1 < minN
                fprintf('Skipping %s: groups have %d / %d events (minimum %d)\n', desc, nE0, nE1, minN);
                continue;
            end

            spec0 = struct('start', start0, 'nBins', nBins);
            spec1 = struct('start', start1, 'nBins', nBins);
            binCenters = binStart + ((1:nBins) - 0.5) * binSize;

            %% ---- features [nEvents x (nBands*nBins)], band-fastest ----
            if strcmp(design, 'paired')
                [F0, F1] = featuresPaired(E0, spec0, E1, spec1, X);
            else
                F0 = featuresOf(E0, spec0, X);
                F1 = featuresOf(E1, spec1, X);
            end

            %% ---- class matrices ----
            if strcmp(design, 'unpaired')
                ok0 = find(~any(isnan(F0), 2));
                ok1 = find(~any(isnan(F1), 2));
                nUse = min(numel(ok0), numel(ok1));
                if nUse < minN
                    fprintf('Skipping %s: valid events %d / %d (minimum %d)\n', desc, numel(ok0), numel(ok1), minN);
                    continue;
                end
                fprintf('Processing: %s (valid %d / %d -> %d random events from each)\n', desc, numel(ok0), numel(ok1), nUse);
                rng(1);
                c0All = F0(ok0(randperm(numel(ok0), nUse)), :);
                c1All = F1(ok1(randperm(numel(ok1), nUse)), :);
                nPerClass = nUse;
            else
                c0All = F0; c1All = F1;
                nPerClass = sum(~any(isnan(F0), 2) & ~any(isnan(F1), 2));
                fprintf('Processing: %s (N=%d pairs, %d complete)\n', desc, size(F0, 1), nPerClass);
            end

            %% ---- SVM on all bands (+ optionally 5 leave-one-band-out SVMs) ----
            if runAblation, nModels = 1 + nBands; else, nModels = 1; end
            for m = 1:nModels
                if m == 1
                    keepBands = 1:nBands;  modelLabel = "All bands";
                else
                    keepBands = setdiff(1:nBands, m-1);  modelLabel = "w/o " + bandShort(m-1);
                end
                % features are band-fastest: column = (bin-1)*nBands + band
                colMask = repmat(ismember(1:nBands, keepBands), 1, nBins);

                pVal = NaN;
                rng(1);    % identical splits for every model
                tic;
                [accAll, bKeep, nullAcc] = trainTestSVM(c0All(:, colMask), c1All(:, colMask), ...
                                                        trainFraction, nSplits, nPermEff * (m == 1), minN);
                if m == 1 && runPermutations
                    pVal = (sum(nullAcc >= accAll) + 1) / (nPermEff + 1);
                    fprintf('   [%s] acc = %.3f | null mean = %.3f | p = %.3f | %.0f s\n', ...
                            modelLabel, accAll, mean(nullAcc, 'omitnan'), pVal, toc);
                else
                    fprintf('   [%s] acc = %.3f | %.0f s\n', modelLabel, accAll, toc);
                end

                % put weights back into a full [nBands x nBins] matrix (removed band = NaN)
                bFull = nan(nBands, nBins);
                bFull(keepBands, :) = reshape(bKeep, numel(keepBands), nBins);

                runCounter = runCounter + 1;
                batchResults(runCounter).Animal     = hTableName;
                batchResults(runCounter).RunLabel   = sprintf('%s (%s) %s', hTableName, upper(cType), modelLabel);
                batchResults(runCounter).TrialType  = cType;
                batchResults(runCounter).ModelIdx   = m;
                batchResults(runCounter).ModelLabel = char(modelLabel);
                batchResults(runCounter).acc        = accAll;
                batchResults(runCounter).beta       = bFull;
                batchResults(runCounter).nullAcc    = nullAcc;
                batchResults(runCounter).nPerClass  = nPerClass;
                batchResults(runCounter).betaScaled = bFull * weightScaleFcn(accAll);   % weight x accuracy
                batchResults(runCounter).pVal       = pVal;
                batchResults(runCounter).binCenters = binCenters;
                batchResults(runCounter).xLab       = xLab;
                batchResults(runCounter).design     = design;
                batchResults(runCounter).class0     = class0Name;
                batchResults(runCounter).class1     = class1Name;
            end
        end
    end

    %% Verification of the number of SVMs
    if ~isempty(batchResults)
        if runAblation
            assert(max([batchResults.ModelIdx]) == 1 + nBands, 'Verification failed: expected %d models per dataset.', 1 + nBands);
            fprintf('Verification: ablation ON -> %d SVMs per dataset (all bands + %d leave-one-band-out), %d SVMs in total.\n', ...
                    1 + nBands, nBands, numel(batchResults));
        else
            assert(all([batchResults.ModelIdx] == 1), 'Verification failed: extra SVMs were run although runAblation = false.');
            fprintf('Verification: ablation OFF -> only the all-bands SVM was run (%d SVMs in total).\n', numel(batchResults));
        end
    end

    %% Plot + save
    if ~isempty(batchResults)
        P = struct('mode', analysisMode, 'bandNames', bandNames, 'binCenters', binCenters, ...
                   'xLab', xLab, 'weightScaleFcn', weightScaleFcn, 'barUseAbs', barUseAbs, ...
                   'nPerm', nPermEff, 'runPermutations', runPermutations, 'nSplits', nSplits, ...
                   'animalColor', animalColor, 'modelColors', modelColors, 'plotWeights', plotWeights);
        P.types  = types;
        P.groups = groups;
        plotModeResults(batchResults, P);
        for r = 1:numel(batchResults)                       % all-bands model only -> summary figure
            if batchResults(r).ModelIdx ~= 1, continue; end
            summaryAll(end+1) = struct('Base', baseMode, 'Strike', isStrike, 'TrialType', batchResults(r).TrialType, ...
                                       'Animal', batchResults(r).Animal, 'acc', batchResults(r).acc); %#ok<SAGROW>
        end
        save(sprintf('svm_allbands_rawnorm_%s.mat', analysisMode), 'batchResults', 'bandNames', ...
             'binCenters', 'analysisMode', 'winLen', 'preStart', 'onsetStart', 'offsetStart', ...
             'nSplits', 'nPermEff', 'runPermutations', 'runAblation', 'normToSum', 'zMode', 'meadowState', ...
             'trainFraction', 'minTrialsRequired', 'minTrialsStrike', 'minTrialsSound', 'strikeMinTime_s', ...
             'preOffsetTypes', 'bandRanges');
    end
end


%% ---- summary figure across analyses ----
if plotSummaryFig && ~isempty(summaryAll)
    plotSummary(summaryAll, animalColor, nSplits);
    save('svm_summary_accuracy.mat', 'summaryAll');
end

%% ---- all figures -> one PowerPoint ----
if savePPT
    saveFiguresToPPT(pptFile);
end

%% ============================================================
%  HELPER FUNCTIONS
% =============================================================

function [needBug, needMeadow, needSound, splitBug] = modeNeeds(mode)
% mode may carry the suffix '_strike' (allowed only for the four bug analyses below)
    isS  = endsWith(mode, '_strike');
    base = erase(mode, '_strike');
    allModes = {'pre_vs_onset', 'pre_vs_offset', 'onset_vs_offset', 'circle_vs_other', ...
                'bug_vs_meadow_onset', 'bug_vs_meadow_offset', 'meadow_move_vs_stat', ...
                'bug_vs_sound_onset', 'bug_vs_sound_offset', 'category_pairs'};
    strikeOK = {'pre_vs_onset', 'pre_vs_offset', 'onset_vs_offset', 'circle_vs_other'};
    if ~ismember(base, allModes), error('Unknown analysisMode: %s', mode); end
    if isS && ~ismember(base, strikeOK), error('No strike variant for analysisMode: %s', mode); end
    needBug    = ~strcmp(base, 'meadow_move_vs_stat');
    needMeadow = ismember(base, {'bug_vs_meadow_onset', 'bug_vs_meadow_offset', 'meadow_move_vs_stat', 'category_pairs'});
    needSound  = ismember(base, {'bug_vs_sound_onset', 'bug_vs_sound_offset', 'category_pairs'});
    splitBug   = needBug && ~ismember(base, {'circle_vs_other', 'category_pairs'});   % circle / other analysed separately
end


function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end


function tbl = bugSubset(H, cType)
% circle (circle, circle_flip) or other bug trials without a strike
    isCircle = (H.Type == "circle" | H.Type == "circle_flip");
    isRest   = isnan(H.Strike_ms);
    if strcmp(cType, 'circle'), tbl = H(isCircle & isRest, :);
    else,                       tbl = H(~isCircle & isRest, :); end
end


function tbl = bugStrikeSubset(H, cType, strikeMin_s)
% circle / other bug trials WITH a strike (Strike_ms not NaN) and Strike_ms > strikeMin_s seconds.
    isCircle  = (H.Type == "circle" | H.Type == "circle_flip");
    enough    = H.Strike_ms > strikeMin_s * 1000;      % NaN strike -> false
    if strcmp(cType, 'circle'), selType = isCircle; else, selType = ~isCircle; end
    tbl = H(selType & enough, :);
end


function tbl = bugTrials(H, cType, isStrike, strikeMin_s)
% no-strike trials (isStrike = false) or strike trials with Strike_ms > strikeMin_s (isStrike = true)
    if ~isStrike, tbl = bugSubset(H, cType); return; end
    tbl = bugStrikeSubset(H, cType, strikeMin_s);
    isCircle = (H.Type == "circle" | H.Type == "circle_flip");
    if strcmp(cType, 'circle'), selType = isCircle; else, selType = ~isCircle; end
    nAll = sum(selType & ~isnan(H.Strike_ms));
    fprintf('   strike trials (%s): %d total, %d kept (Strike_ms > %g s)\n', cType, nAll, height(tbl), strikeMin_s);
end


%% ---------------- event tables (with trial info for zMode 'trial') ----------------

function E = hunterSet(tbl, T)
% one row per hunter trial; refs = [onset onset offset offset], wStarts = T.hunterStarts
% (pre-offset column = NaN / not applicable for trial types not in T.preOffsetTypes)
    n   = height(tbl);
    on  = tbl.Triger_ms(:);  off = tbl.End_ms(:);
    isC = (tbl.Type == "circle" | tbl.Type == "circle_flip");  isC = isC(:);
    hasPO = (isC & ismember('circle', T.preOffsetTypes)) | (~isC & ismember('other', T.preOffsetTypes));
    refs = [on on off off];             refs(~hasPO, 3) = NaN;
    wSt  = repmat(T.hunterStarts, n, 1); wSt(~hasPO, 3) = NaN;
    E = table(string(tbl.Animal), string(tbl.Record), nan(n, 1), refs, wSt, ...
              'VariableNames', {'Animal', 'Record', 't_ms', 'refs', 'wStarts'});
end


function E = hunterEvents(tbl, align, T)
% onset = Triger_ms, offset = End_ms; refs / wStarts = all windows of the trial (for zMode 'trial')
    E = hunterSet(tbl, T);
    if strcmp(align, 'offset'), E.t_ms = tbl.End_ms(:); else, E.t_ms = tbl.Triger_ms(:); end
end


function E = meadowSet(M, animalID, T)
% one row per meadow row; refs = [moveOn moveOn statOn statOn], windows [pre-move, onset-move, pre-stat, onset-stat]
    n = height(M);
    if ismember('Animal', M.Properties.VariableNames), an = string(M.Animal);
    else, an = repmat("PV" + animalID, n, 1); end
    statOn = M.Triger_ms(:);  moveOn = M.Move_ms(:);
    E = table(an(:), string(M.Record), nan(n, 1), [moveOn moveOn statOn statOn], repmat(T.meadowStarts, n, 1), ...
              'VariableNames', {'Animal', 'Record', 't_ms', 'refs', 'wStarts'});
end


function E = meadowEvents(M, animalID, evtName, T)
% evtName: 'statOn' | 'statOff' | 'moveOn' | 'moveOff'
%   statOn = Triger_ms, moveOn = Move_ms, statOff = Move_ms of the same row,
%   moveOff = statOn of the NEXT row (same record, else NaN)
% The offset events are not one of the row's windows -> no trial info (own-window z-score in zMode 'trial').
    E = meadowSet(M, animalID, T);
    statOn = E.refs(:, 4);  moveOn = E.refs(:, 2);
    nextOn = [statOn(2:end); NaN];
    sameRec = [E.Record(1:end-1) == E.Record(2:end) & E.Animal(1:end-1) == E.Animal(2:end); false];
    moveOff = nextOn;  moveOff(~sameRec) = NaN;
    switch evtName
        case 'statOn',  t = statOn;
        case 'statOff', t = moveOn;
        case 'moveOn',  t = moveOn;
        case 'moveOff', t = moveOff;
        otherwise, error('Unknown meadow event: %s', evtName);
    end
    E.t_ms = t(:);
    if ismember(evtName, {'statOff', 'moveOff'}), E = removevars(E, {'refs', 'wStarts'}); end
end


function E = soundEvents(S, animalID, recFilter, align, soundDur_ms, exclWin_s, T)
% wake sound events of one animal, onset = Triger_ms, offset = onset + soundDur_ms.
% Dropped if ANOTHER sound onset (same record) lies within +-exclWin_s s of the event time.
    animalNum = cellfun(@(x) str2double(regexp(char(x), '\d+', 'match', 'once')), cellstr(string(S.Animal)));
    sel = animalNum(:) == animalID & logical(recFilter(S.Record));
    S   = S(sel, :);
    rec = string(S.Record);
    tOn = S.Triger_ms(:);
    if strcmp(align, 'offset'), t = tOn + soundDur_ms; else, t = tOn; end
    keep = true(height(S), 1);
    for i = 1:height(S)
        other = rec == rec(i);  other(i) = false;
        if any(abs(tOn(other) - t(i)) <= exclWin_s * 1000), keep(i) = false; end
    end
    nK   = nnz(keep);
    refs = [tOn(keep) tOn(keep) tOn(keep) + soundDur_ms];
    E = table(string(S.Animal(keep)), rec(keep), t(keep), refs, repmat(T.soundStarts, nK, 1), ...
              'VariableNames', {'Animal', 'Record', 't_ms', 'refs', 'wStarts'});
end


function E = soundSet(S, animalID, T)
% wake sounds of one animal; windows [pre, onset, offset]; a sound is dropped if another sound interval
% [onset, onset + dur] of the same record overlaps ANY of its 3 windows
    num = cellfun(@(x) str2double(regexp(char(x), '\d+', 'match', 'once')), cellstr(string(S.Animal)));
    S   = S(num(:) == animalID & logical(T.soundRecordFilter(S.Record)), :);
    rec = string(S.Record);  tOn = S.Triger_ms(:);
    refs = [tOn tOn tOn + T.soundDur_ms];
    keep = true(height(S), 1);
    for i = 1:height(S)
        oth = find(rec == rec(i));  oth(oth == i) = [];
        for w = 1:3
            wS = refs(i, w) + T.soundStarts(w) * 1000;  wE = wS + T.winLen * 1000;
            if any(tOn(oth) < wE & tOn(oth) + T.soundDur_ms > wS), keep(i) = false; break; end
        end
    end
    nK = nnz(keep);
    E = table(string(S.Animal(keep)), rec(keep), nan(nK, 1), refs(keep, :), repmat(T.soundStarts, nK, 1), ...
              'VariableNames', {'Animal', 'Record', 't_ms', 'refs', 'wStarts'});
end


function E = emptySetTable(nW)
    E = table(strings(0, 1), strings(0, 1), zeros(0, 1), zeros(0, nW), zeros(0, nW), ...
              'VariableNames', {'Animal', 'Record', 't_ms', 'refs', 'wStarts'});
end


function sets = buildTrialSets(H, M, Snd, animalID, T)
% all trial sets of ONE animal: circle / other (no strike), circleS / otherS (strike), meadow, sound
    e4 = struct('E', emptySetTable(4), 'starts', T.hunterStarts);
    sets = struct('circle', e4, 'other', e4, 'circleS', e4, 'otherS', e4, ...
                  'meadow', struct('E', emptySetTable(4), 'starts', T.meadowStarts), ...
                  'sound',  struct('E', emptySetTable(3), 'starts', T.soundStarts));
    if ~isempty(H)
        sets.circle.E  = hunterSet(bugTrials(H, 'circle', false, T.strikeMin_s), T);
        sets.other.E   = hunterSet(bugTrials(H, 'other',  false, T.strikeMin_s), T);
        sets.circleS.E = hunterSet(bugTrials(H, 'circle', true,  T.strikeMin_s), T);
        sets.otherS.E  = hunterSet(bugTrials(H, 'other',  true,  T.strikeMin_s), T);
    end
    if ~isempty(M),   sets.meadow.E = meadowSet(M, animalID, T); end
    if ~isempty(Snd), sets.sound.E  = soundSet(Snd, animalID, T); end
end


%% ---------------- category_pairs: definitions of all comparisons ----------------

function D = pairDefs(preOffsetTypes)
% all comparisons: paired = every couple of windows of one trial set; unpaired = couples of sets
% (onset / offset figures) and strike vs no-strike (separately for circle / other and onset / offset)
    D = struct('group', {}, 'name', {}, 'setA', {}, 'winA', {}, 'labelA', {}, 'setB', {}, 'winB', {}, 'labelB', {});
    labs4 = {'Baseline', 'Onset', 'Pre-offset', 'Offset'};   labs3 = {'Baseline', 'Onset', 'Offset'};
    if ismember('circle', preOffsetTypes), cL = labs4; cW = [1 2 3 4]; else, cL = labs3; cW = [1 2 4]; end
    if ismember('other',  preOffsetTypes), oL = labs4; oW = [1 2 3 4]; else, oL = labs3; oW = [1 2 4]; end
    paired = {'circle',        'circle',  cL, cW
              'circle strike', 'circleS', cL, cW
              'other',         'other',   oL, oW
              'other strike',  'otherS',  oL, oW
              'meadow',        'meadow',  {'Pre move', 'Onset move', 'Pre stat', 'Onset stat'}, [1 2 3 4]
              'sound',         'sound',   labs3, [1 2 3]};
    for i = 1:size(paired, 1)
        L = paired{i, 3};  W = paired{i, 4};
        for a = 1:numel(W) - 1
            for b = a + 1:numel(W)
                D(end+1) = mkDef(paired{i, 1}, paired{i, 2}, W(a), L{a}, paired{i, 2}, W(b), L{b}); %#ok<AGROW>
            end
        end
    end
    unp = {'onset',         {'circle', 'other', 'meadow', 'sound'}, [2 2 2 2], {'Circle bug', 'Other bug', 'Meadow move', 'Sound'}
           'offset',        {'circle', 'other', 'sound'},           [4 4 3],   {'Circle bug', 'Other bug', 'Sound'}
           'onset strike',  {'circleS', 'otherS'},                  [2 2],     {'Circle bug', 'Other bug'}
           'offset strike', {'circleS', 'otherS'},                  [4 4],     {'Circle bug', 'Other bug'}};
    for i = 1:size(unp, 1)
        S = unp{i, 2};  W = unp{i, 3};  L = unp{i, 4};
        for a = 1:numel(S) - 1
            for b = a + 1:numel(S)
                D(end+1) = mkDef(unp{i, 1}, S{a}, W(a), L{a}, S{b}, W(b), L{b}); %#ok<AGROW>
            end
        end
    end

    % strike vs no-strike, separately for circle / other and separately for onset / offset
    % (unpaired: equal random subsets per animal; window 2 = onset, window 4 = offset)
    for tp = {'circle', 'other'}
        tp = tp{1};
        D(end+1) = mkDef([tp ' strike vs no-strike'], tp, 2, 'No-strike Onset',  [tp 'S'], 2, 'Strike Onset');  %#ok<AGROW>
        D(end+1) = mkDef([tp ' strike vs no-strike'], tp, 4, 'No-strike Offset', [tp 'S'], 4, 'Strike Offset'); %#ok<AGROW>
    end
end


function d = mkDef(group, setA, winA, labA, setB, winB, labB)
    d = struct('group', group, 'name', sprintf('%s | %s vs %s', group, labA, labB), ...
               'setA', setA, 'winA', winA, 'labelA', labA, 'setB', setB, 'winB', winB, 'labelB', labB);
end


function G = pairGroups(D)
% one figure per group: G(g).name, G(g).types = names of the comparisons of that group
    [~, ia] = unique({D.group}, 'stable');
    names = {D(ia).group};
    G = struct('name', names, 'types', cell(1, numel(names)));
    for g = 1:numel(names)
        G(g).types = {D(strcmp({D.group}, names{g})).name};
    end
end


function [E0, E1, start0, start1, design, nOverlap] = buildPair(d, sets, T)
% event tables of the two classes of one comparison.
% Same set  -> PAIRED; trials whose two windows overlap in time are removed from this comparison only.
% Two sets  -> UNPAIRED (equal random subsets are drawn later in the main loop).
    A = sets.(d.setA);  B = sets.(d.setB);
    start0 = A.starts(d.winA);  start1 = B.starts(d.winB);
    nOverlap = 0;
    if strcmp(d.setA, d.setB)
        design = 'paired';
        ov = windowsOverlap(A.E.refs(:, d.winA), start0, A.E.refs(:, d.winB), start1, T.winLen);
        nOverlap = sum(ov);
        Ek = A.E(~ov, :);
        E0 = Ek;  E0.t_ms = Ek.refs(:, d.winA);
        E1 = Ek;  E1.t_ms = Ek.refs(:, d.winB);
    else
        design = 'unpaired';
        E0 = A.E;  E0.t_ms = A.E.refs(:, d.winA);
        E1 = B.E;  E1.t_ms = B.E.refs(:, d.winB);
    end
end


function ov = windowsOverlap(r1, s1, r2, s2, winLen)
% windows [r + s*1000, r + s*1000 + winLen*1000) ms; true where the two windows intersect
    a1 = r1 + s1 * 1000;  a2 = r2 + s2 * 1000;
    ov = a1 < a2 + winLen * 1000 & a2 < a1 + winLen * 1000;
end


%% ---------------- zMode 'trial' helpers ----------------

function tf = hasTrialRef(E)
    tf = all(ismember({'refs', 'wStarts'}, E.Properties.VariableNames));
end


function [mu, sd, ok] = trialStats(R, refs, starts, winLen)
% mean / std per band over ALL windows of ONE trial (union of their time points, overlapping samples counted once).
% starts(w) = NaN -> window not applicable for this trial (skipped); refs(w) = NaN or a window without data -> ok = false.
    mu = [];  sd = [];  ok = false;
    use = false(size(R.t_ms));
    for w = 1:numel(refs)
        if isnan(starts(w)), continue; end
        if isnan(refs(w)), return; end
        rel = (R.t_ms - refs(w)) / 1000;
        in  = rel >= starts(w) & rel < starts(w) + winLen;
        if ~any(in), return; end
        use = use | in;
    end
    if ~any(use), return; end
    both = R.bandRaw(:, use);
    mu = mean(both, 2, 'omitnan');
    sd = std(both, 0, 2, 'omitnan');
    sd(sd == 0) = NaN;
    ok = true;
end


%% ---------------- features ----------------

function F = featuresOf(E, spec, X)
% [nEvents x (nBands*nBins)], features band-fastest (-> reshapes to [nBands x nBins])
    D = extractLFPBins(E, X.SA, X.bandRanges, spec, X.binSize, X.recCache, X.normToSum, X.zMode);
    F = reshape(D, size(D, 1) * size(D, 2), [])';
end


function [F0, F1] = featuresPaired(E0, spec0, E1, spec1, X)
% Features of a PAIRED design (row i of E0 and row i of E1 = same trial / same table row).
% zMode 'window': the two windows of the same trial are z-scored TOGETHER (one mean/std per band over both).
% zMode 'trial' : the statistics come from ALL windows of the trial (given in E.refs), so each class is
%                 processed separately by featuresOf. Without trial info -> joint 'window' z-score.
% other zModes  : identical to featuresOf on each class.
    useJoint = strcmp(X.zMode, 'window') || (strcmp(X.zMode, 'trial') && ~hasTrialRef(E0));
    if ~useJoint
        F0 = featuresOf(E0, spec0, X);
        F1 = featuresOf(E1, spec1, X);
        return;
    end

    nBands = size(X.bandRanges, 1);
    nB0 = spec0.nBins;  nB1 = spec1.nBins;
    edges0 = spec0.start + (0:nB0) * X.binSize;
    edges1 = spec1.start + (0:nB1) * X.binSize;
    nTr = height(E0);
    D0 = nan(nBands, nB0, nTr);
    D1 = nan(nBands, nB1, nTr);

    for tr = 1:nTr
        animal = char(E0.Animal(tr));
        record = char(E0.Record(tr));
        if ~(strcmp(animal, char(E1.Animal(tr))) && strcmp(record, char(E1.Record(tr)))), continue; end

        t0 = E0.t_ms(tr);  t1 = E1.t_ms(tr);
        if isnan(t0) || isnan(t1), continue; end

        R = getRec(X.recCache, X.SA, animal, record, X.bandRanges, X.normToSum);
        if isempty(R), continue; end

        rel0 = (R.t_ms - t0) / 1000;
        rel1 = (R.t_ms - t1) / 1000;
        w0 = find(rel0 >= edges0(1) & rel0 < edges0(end));
        w1 = find(rel1 >= edges1(1) & rel1 < edges1(end));
        if isempty(w0) || isempty(w1), continue; end

        seg0 = R.bandRaw(:, w0);
        seg1 = R.bandRaw(:, w1);
        both = [seg0, seg1];                         % [nBands x (n0 + n1)]
        mu = mean(both, 2, 'omitnan');
        sd = std(both, 0, 2, 'omitnan');
        sd(sd == 0) = NaN;                           % flat signal -> NaN instead of Inf

        D0(:, :, tr) = binSeg((seg0 - mu) ./ sd, rel0(w0), edges0);
        D1(:, :, tr) = binSeg((seg1 - mu) ./ sd, rel1(w1), edges1);
    end

    F0 = reshape(D0, nBands * nB0, [])';
    F1 = reshape(D1, nBands * nB1, [])';
end


function out = binSeg(seg, sub, tEdges)
% average seg [nBands x nTime] (time stamps sub) into the bins defined by tEdges -> [nBands x nBins]
    nBins = numel(tEdges) - 1;
    out = nan(size(seg, 1), nBins);
    for tb = 1:nBins
        idx = sub >= tEdges(tb) & sub < tEdges(tb+1);
        if any(idx)
            out(:, tb) = mean(seg(:, idx), 2, 'omitnan');   % all bands at once
        end
    end
end


function R = getRec(recCache, SA, animal, record, bandRanges, normToSum)
% band-averaged traces of a recording (loaded only once, then taken from the cache)
    key = [animal '|' record];
    if ~isKey(recCache, key)
        recCache(key) = loadBandZ(SA, animal, record, bandRanges, normToSum);
    end
    R = recCache(key);
end


%% ---------------- SVM ----------------

function [acc, beta, nullAcc] = trainTestSVM(c0Mat, c1Mat, trainFraction, nSplits, nPerm, minTrials)
% Paired linear SVM (no standardization inside the SVM).
% Row i of c0Mat (class 0) and row i of c1Mat (class 1) form a pair that always stays together
% in train or test.
% acc, beta : mean over nSplits random splits on the REAL labels
% nullAcc   : [nPerm x 1]; each entry = mean over nSplits random splits on one permuted dataset
    acc     = nan;
    beta    = nan(1, size(c0Mat, 2));
    nullAcc = nan(nPerm, 1);

    good = find(~any(isnan(c0Mat), 2) & ~any(isnan(c1Mat), 2));
    if numel(good) < minTrials, return; end
    X0 = c0Mat(good, :);
    X1 = c1Mat(good, :);
    nPairs = numel(good);

    % real data (no swapping)
    [acc, beta] = svmOverSplits(X0, X1, false(nPairs, 1), trainFraction, nSplits);

    % null: swap labels within each pair with prob 0.5 (a fresh swap pattern per
    % permutation, kept fixed across that permutation's nSplits splits)
    for p = 1:nPerm
        swap = rand(nPairs,1) < 0.5;
        nullAcc(p) = svmOverSplits(X0, X1, swap, trainFraction, nSplits);
    end
end


function [acc, beta] = svmOverSplits(X0, X1, swap, trainFraction, nSplits)
% Run the SVM on nSplits random train/test splits of the pairs and average.
% swap(i) = true : pair i has labels flipped (X0 row -> 1, X1 row -> 0)
    nPairs = size(X0, 1);
    nTrain = floor(trainFraction * nPairs);
    y0 = double(swap);          % label of the X0 row of each pair
    y1 = double(~swap);         % label of the X1 row of each pair

    accs  = nan(nSplits, 1);
    betas = nan(nSplits, size(X0, 2));
    for s = 1:nSplits
        order = randperm(nPairs);
        tr = order(1:nTrain); te = order(nTrain+1:end);
        Xtr = [X0(tr, :); X1(tr, :)];   ytr = [y0(tr); y1(tr)];
        Xte = [X0(te, :); X1(te, :)];   yte = [y0(te); y1(te)];
        model = fitcsvm(Xtr, ytr, 'KernelFunction', 'linear', 'ClassNames', [0 1]);
        accs(s)     = mean(predict(model, Xte) == yte);
        betas(s, :) = model.Beta';
    end
    acc  = mean(accs);
    beta = mean(betas, 1);
end


%% ---------------- LFP extraction ----------------

function binnedData = extractLFPBins(evt, SA,bandRanges, spec, binSize, recCache, normToSum, zMode)
% evt        : event table with Animal, Record, t_ms (event time in recording ms)
%              (+ refs, wStarts for zMode 'trial')
% spec.start : window start (s) relative to the event;  spec.nBins : number of bins
%              window = [start, start + nBins*binSize)  (end point not included)
% zMode      : 'recording' -> recording-wide z-scored traces
%              'window'    -> z-score each band within this event's window ALONE, then average into bins
%              'trial'     -> z-score with mean/std over ALL windows of the same trial (evt.refs / evt.wStarts);
%                             events without trial info use their own window
%              'none'      -> raw traces, no z-scoring
    nBands  = size(bandRanges, 1);
    nBins   = spec.nBins;
    tEdges  = spec.start + (0:nBins) * binSize;
    nTrials = height(evt);
    binnedData = nan(nBands, nBins, nTrials);
    hasTrial = hasTrialRef(evt);

    for tr = 1:nTrials
        animal = char(evt.Animal(tr));
        record = char(evt.Record(tr));

        

        alignTime_ms = evt.t_ms(tr);
        if isnan(alignTime_ms), continue; end

        R = getRec(recCache, SA, animal, record, bandRanges, normToSum);
        if isempty(R), continue; end        % recording / LFP file not found

        relTime = (R.t_ms - alignTime_ms) / 1000;
        w = find(relTime >= tEdges(1) & relTime < tEdges(end));   % only this event's window
        if isempty(w), continue; end
        sub = relTime(w);

        switch zMode
            case 'window'
                % [nBands x nTimeInWindow]: mean/std over TIME, separately for each band
                seg = R.bandRaw(:, w);
                sd  = std(seg, 0, 2, 'omitnan');
                sd(sd == 0) = NaN;                   % flat signal -> NaN instead of Inf
                seg = (seg - mean(seg, 2, 'omitnan')) ./ sd;
            case 'trial'
                if hasTrial
                    % mean/std per band over ALL windows of this trial (same normalization in every comparison)
                    [mu, sd, ok] = trialStats(R, evt.refs(tr, :), evt.wStarts(tr, :), nBins * binSize);
                    if ~ok, continue; end           % a window of this trial has no data -> trial skipped everywhere
                    seg = (R.bandRaw(:, w) - mu) ./ sd;
                else                                 % event without trial info -> its own window
                    seg = R.bandRaw(:, w);
                    sd  = std(seg, 0, 2, 'omitnan');
                    sd(sd == 0) = NaN;
                    seg = (seg - mean(seg, 2, 'omitnan')) ./ sd;
                end
            case 'recording'
                seg = R.bandZ(:, w);
            case 'none'
                seg = R.bandRaw(:, w);
            otherwise
                error('Unknown zMode: %s', zMode);
        end

        binnedData(:, :, tr) = binSeg(seg, sub, tEdges);
    end
end


function R = loadBandZ(SA, animal, record, bandRanges, normToSum)
% Load one recording and return t_ms plus band-averaged activity [nBands x nTime]:
%   R.bandRaw : band-averaged raw power (after optional normToSum)
%   R.bandZ   : the same, z-scored per band over the whole recording
% Returns [] if the recording / LFP file cannot be found.
    n = find(strcmp(string(SA.recTable.Animal), string(animal)) & strcmp(string(SA.recTable.recNames), string(record)));
    if isempty(n)
        warning('Recording %s / %s not found in recTable - events skipped.', animal, record);
        R = []; return;
    end
    folder = fullfile(SA.recTable.folder{n}, 'analysis', ['Animal=' animal ',recNames=' record]);
    ch = SA.recTable.defaulLFPCh(n);

    if strcmp(animal, 'PV106') && ismember(record, {'Hunter19', 'Hunter20', 'Hunter21', 'Hunter34'}), ch = 5;
    elseif strcmp(animal, 'PV157') && ismember(record, {'Hunter40', 'Hunter41'}), ch = 17; end

    fileName = [folder filesep 'bgRatio_ch' num2str(ch) '.mat'];
    if ~isfile(fileName)
        warning('File not found: %s - events skipped.', fileName);
        R = []; return;
    end
    S = load(fileName, 'f', 't_ms', 'allFreqProfiles');

    nBands = size(bandRanges, 1);
    B = nan(nBands, numel(S.t_ms));
    for fb = 1:nBands
        freqIdx = S.f >= bandRanges(fb,1) & S.f <= bandRanges(fb,2);
        if any(freqIdx), B(fb, :) = mean(S.allFreqProfiles(freqIdx, :), 1, 'omitnan'); end
    end
    % normalize raw mean power of each band to the sum of the bands' mean powers at every time point
    if normToSum
        B = B ./ sum(B, 1, 'omitnan');
        B(~isfinite(B)) = NaN;        % sum == 0 / all-NaN time points
    end

    R.bandRaw = B;
    R.bandZ   = (B - mean(B, 2, 'omitnan')) ./ std(B, 0, 2, 'omitnan');   % whole-recording z-score (NaN tolerant)
    R.t_ms    = S.t_ms(:)';
end


%% ---------------- plotting ----------------

function plotModeResults(R, P)
    nBands  = numel(P.bandNames);
    nBins   = size(R(1).beta, 2);
    modeStr = upper(strrep(P.mode, '_', ' '));

    if isfield(P, 'groups') && ~isempty(P.groups)
        for g = 1:numel(P.groups)                       % one accuracy figure per group of comparisons
            Rg = R(ismember({R.TrialType}, P.groups(g).types));
            if isempty(Rg), continue; end
            Pg = P;  Pg.mode = ['category pairs: ' P.groups(g).name];  Pg.types = P.groups(g).types;
            plotAccuracy(Rg, Pg);
        end
    else
        plotAccuracy(R, P);
    end
    if ~P.plotWeights, return; end

    %% --- Weights: separate figures for every panel AND every SVM ---
    for tt = P.types
        for m = unique([R.ModelIdx])
            sel = find(strcmp({R.TrialType}, tt{1}) & [R.ModelIdx] == m);
            if isempty(sel), continue; end
            animals = {R(sel).Animal};
            W = nan(nBins, nBands, numel(sel));       % [time bins x bands x animals]
            for i = 1:numel(sel)
                r = sel(i);
                W(:, :, i) = R(r).beta' * P.weightScaleFcn(R(r).acc);
            end
            Pw = P;  Pw.binCenters = R(sel(1)).binCenters;  Pw.xLab = R(sel(1)).xLab;   % axis of THIS comparison
            ttl = sprintf('%s | SVM: %s | per animal (weight x accuracy, mean of %d splits) | %s', ...
                          upper(strrep(tt{1}, '_', ' ')), R(sel(1)).ModelLabel, P.nSplits, modeStr);
            plotWeightsShifted(W, Pw, animals, ttl);   % one curve per animal + mean
            plotWeightsBars(W, Pw, animals, ttl);      % one row per animal + mean
        end
    end
end


function plotAccuracy(R, P)
% x = animals (+ Mean); at each x one star per model (All bands + each leave-one-band-out), each model
% in its own color. Vertical grid lines are drawn BETWEEN the animals, tick labels stay at the animal centers.
% Null boxplot / p-value (if P.runPermutations) is shown for the All-bands model only.
    types   = P.types;
    nT      = numel(types);
    animals = unique({R.Animal}, 'stable');
    nA      = numel(animals);
    showMean = nA > 1;
    nX      = nA + showMean;
    labels  = animals;
    if showMean, labels{end+1} = 'Mean'; end
    modeStr = upper(strrep(P.mode, '_', ' '));

    models = unique([R.ModelIdx]);
    nM     = numel(models);
    if nM > 1, offs = linspace(-0.33, 0.33, nM); else, offs = 0; end
    cols = P.modelColors(1:nM, :);       % distinct color per model
    mLab = cell(1, nM);
    for m = 1:nM, mLab{m} = R(find([R.ModelIdx] == models(m), 1)).ModelLabel; end

    nCol = min(nT, 3);  nRow = ceil(nT / nCol);
    figure('Name', ['SVM accuracy ' P.mode], 'Color', 'w', 'Position', [100 100 300 + 450 * nCol, 80 + 420 * nRow]);
    tl = tiledlayout(nRow, nCol, 'TileSpacing', 'compact');
    hLeg = []; axLeg = [];
    for k = 1:nT
        ax = nexttile; hold(ax, 'on');
        title(ax, upper(strrep(types{k}, '_', ' ')), 'FontSize', 9, 'Interpreter', 'none');
        ylim(ax, [0 1]); xlim(ax, [0.5 nX + 0.5]); box(ax, 'on');
        grid(ax, 'off'); ax.YGrid = 'on';                       % horizontal grid only
        for xb = 1.5:1:(nX - 0.5)                               % vertical lines BETWEEN the animals
            xline(ax, xb, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        end
        xticks(ax, 1:nX); xticklabels(ax, labels); ax.TickLabelInterpreter = 'none';
        yline(ax, 0.5, 'k--', 'Chance (50%)', 'HandleVisibility', 'off');
        if mod(k - 1, nCol) == 0, ylabel(ax, sprintf('Accuracy (mean of %d splits)', P.nSplits)); end

        selT = find(strcmp({R.TrialType}, types{k}));
        if isempty(selT)
            text(ax, 0.5, 0.5, 'no data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
            continue;
        end

        accM  = nan(nM, nX);
        nullM = nan(max(P.nPerm, 1), nX);
        for m = 1:nM
            for a = 1:nA
                s = selT(strcmp({R(selT).Animal}, animals{a}) & [R(selT).ModelIdx] == models(m));
                if isempty(s), continue; end
                accM(m, a) = R(s).acc;
                if m == 1 && P.runPermutations, nullM(:, a) = R(s).nullAcc; end
            end
            if showMean, accM(m, nX) = mean(accM(m, 1:nA), 'omitnan'); end
        end
        if showMean && P.runPermutations, nullM(:, nX) = mean(nullM(:, 1:nA), 2, 'omitnan'); end

        h = [];
        if P.runPermutations
            xg = repelem((1:nX) + offs(1), P.nPerm);
            hB = boxchart(ax, xg(:), nullM(:), 'BoxFaceColor', [.6 .6 .6], 'MarkerColor', [.4 .4 .4], ...
                          'BoxWidth', 0.12, 'DisplayName', sprintf('Null, all bands (%d perm.)', P.nPerm));
            h = hB;
            pVals = (sum(nullM >= accM(1, :), 1) + 1) / (P.nPerm + 1);
            pVals(isnan(accM(1, :))) = NaN;
        end
        for m = 1:nM
            hm = plot(ax, (1:nA) + offs(m), accM(m, 1:nA), 'p', 'Color', cols(m, :), ...
                      'MarkerFaceColor', cols(m, :), 'MarkerSize', 10, 'DisplayName', mLab{m});
            h = [h hm];
            if showMean
                plot(ax, nX + offs(m), accM(m, nX), 'p', 'Color', 'k', 'MarkerFaceColor', cols(m, :), ...
                     'MarkerSize', 14, 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
        end
        if showMean
            yline(ax, accM(1, nX), ':', 'Color', cols(1, :), 'HandleVisibility', 'off');  % full-model mean level
        end
        if P.runPermutations
            for x = 1:nX
                if ~isnan(pVals(x))
                    text(ax, x + offs(1), min(accM(1, x) + 0.07, 0.96), sprintf('p=%.3f', pVals(x)), ...
                         'HorizontalAlignment', 'center', 'FontSize', 8);
                end
            end
        end
        if isempty(hLeg), hLeg = h; axLeg = ax; end
    end
    if nM > 1
        title(tl, sprintf('SVM accuracy: all bands vs. one band removed (%s)', modeStr));
    else
        title(tl, sprintf('SVM accuracy: all bands (%s)', modeStr));
    end
    if ~isempty(hLeg)
        lgd = legend(axLeg, hLeg);
        lgd.Layout.Tile = 'east';
    end
end


function plotWeightsShifted(W, P, labels, ttl) %#ok<INUSD>
% shifted |weight| time courses, tall figure; all animals in ONE color, black = mean across animals
    nRuns = size(W, 3);
    nBands = size(W, 2);
    absW = abs(W);
    vs = 1.1 * max(absW(:), [], 'omitnan');
    if isempty(vs) || ~isfinite(vs) || vs <= 0, vs = 1; end
    shiftOffsets = 0:vs:(vs * (nBands - 1));

    figure('Name', ['Weights shifted: ' ttl], 'Color', 'w', 'Position', [100 50 650 950]);
    hold on;
    for r = 1:nRuns
        plotShifted(P.binCenters, absW(:, :, r), '-o', 'verticalShift', vs, 'LineWidth', 1.5, ...
                    'MarkerSize', 3, 'Color', P.animalColor, 'HandleVisibility', 'off');
    end
    plot(nan, nan, '-o', 'Color', P.animalColor, 'DisplayName', sprintf('Animals (n=%d)', nRuns));
    if nRuns > 1
        plotShifted(P.binCenters, mean(absW, 3, 'omitnan'), '-s', 'verticalShift', vs, 'LineWidth', 3, ...
                    'MarkerSize', 4, 'Color', 'k', 'HandleVisibility', 'off');
        plot(nan, nan, '-s', 'Color', 'k', 'LineWidth', 3, 'DisplayName', 'Mean');
    end
    yticks(shiftOffsets); yticklabels(P.bandNames);
    xlabel(P.xLab, 'FontSize', 8); ylabel('Frequency band (shifted |weight|)');
    title(ttl, 'FontSize', 9, 'Interpreter', 'none');
    pbaspect([1 2 1]);
    legend('Location', 'southoutside'); grid on; hold off;
end


function plotWeightsBars(W, P, labels, ttl)
% bar plots of weights vs time: ONE ROW PER ANIMAL (+ black Mean row), one column per band,
% SAME y axis everywhere
    nBands = size(W, 2); nBins = size(W, 1); nRuns = size(W, 3);
    if P.barUseAbs, W = abs(W); end
    showMean = nRuns > 1;
    if showMean
        W = cat(3, W, mean(W, 3, 'omitnan'));
        labels = [labels, {'Mean'}];
    end
    nRows = size(W, 3);

    lo = min(0, min(W(:), [], 'omitnan'));
    hi = max(W(:), [], 'omitnan');
    pad = 0.05 * (hi - lo);
    if ~isfinite(pad) || pad == 0, pad = 1; end
    yl = [lo - pad * (lo < 0), hi + pad];     % common y limits

    figure('Name', ['Weights bars: ' ttl], 'Color', 'w', 'Position', [50 50 1700 220 * nRows + 120]);
    tl = tiledlayout(nRows, nBands, 'TileSpacing', 'compact');
    for r = 1:nRows
        if showMean && r == nRows, col = 'k'; else, col = P.animalColor; end
        for b = 1:nBands
            nexttile; hold on;
            bar(P.binCenters, reshape(W(:, b, r), nBins, 1), 'FaceColor', col);
            yline(0, 'k-', 'HandleVisibility', 'off');
            ylim(yl); box on;
            if r == 1, title(P.bandNames(b), 'FontSize', 9); end
            if b == 1, ylabel(labels{r}, 'Interpreter', 'none', 'FontWeight', 'bold'); else, yticklabels([]); end
            if r < nRows, xticklabels([]); end
        end
    end
    xlabel(tl, P.xLab, 'FontSize', 8);
    title(tl, ttl, 'Interpreter', 'none', 'FontSize', 10);
end


function varargout = plotShifted(varargin)
    pNumeric = cellfun(@(x) isnumeric(x), varargin);
    pFirst   = find(pNumeric, 1, 'first');
    pData    = pFirst + (nargin > 1 && pNumeric(pFirst+1));

    localPropertyList = {'verticalShift'};
    for i = 1:numel(localPropertyList)
        pProp = find(strcmp(localPropertyList{i}, varargin));
        if ~isempty(pProp)
            Par.(localPropertyList{i}) = varargin{pProp+1};
            varargin(pProp:pProp+1) = [];
        else
            Par.(localPropertyList{i}) = [];
        end
    end

    [~, nCh] = size(varargin{pData});
    if isempty(Par.verticalShift)
        Par.verticalShift = cumsum(nanstd(varargin{pData}) * 3);
    elseif numel(Par.verticalShift) == 1
        Par.verticalShift = 0:Par.verticalShift:(Par.verticalShift*(nCh-1));
    end

    varargin{pData} = bsxfun(@plus, varargin{pData}, Par.verticalShift);
    [varargout] = plot(varargin{:});
    varargout = {varargout(1:nargout)};
end


function plotSummary(S, animalColor, nSplits)
% x = Pre vs Onset | Pre vs Offset | Onset vs Offset; one line (same color) per animal connecting its dots,
% black = mean over animals. One panel per trial type (columns) and per no-strike / strike group (rows).
    xModes  = {'pre_vs_onset', 'pre_vs_offset', 'onset_vs_offset'};
    xLabels = {'Pre vs Onset', 'Pre vs Offset', 'Onset vs Offset'};
    S = S(ismember({S.Base}, xModes));
    if isempty(S), return; end

    strikeRows = unique(double([S.Strike]));          % 0 = no-strike trials, 1 = strike trials
    types   = {'circle', 'other'};
    nR = numel(strikeRows);  nC = numel(types);
    animals = unique({S.Animal}, 'stable');
    nA = numel(animals);

    figure('Name', 'Summary: accuracy across analyses', 'Color', 'w', 'Position', [100 100 450 * nC + 150, 400 * nR + 100]);
    tl = tiledlayout(nR, nC, 'TileSpacing', 'compact');
    hLeg = []; axLeg = [];
    for r = 1:nR
        if strikeRows(r) == 1, grpStr = 'STRIKE trials'; else, grpStr = 'NO-STRIKE trials'; end
        for c = 1:nC
            ax = nexttile; hold(ax, 'on');
            title(ax, sprintf('%s | %s', upper(types{c}), grpStr));
            ylim(ax, [0 1]); xlim(ax, [0.5 numel(xModes) + 0.5]); grid(ax, 'on'); box(ax, 'on');
            xticks(ax, 1:numel(xModes)); xticklabels(ax, xLabels);
            yline(ax, 0.5, 'k--', 'Chance (50%)', 'HandleVisibility', 'off');
            if c == 1, ylabel(ax, sprintf('Accuracy (mean of %d splits)', nSplits)); end

            sel = double([S.Strike]) == strikeRows(r) & strcmp({S.TrialType}, types{c});
            A = nan(nA, numel(xModes));
            for a = 1:nA
                for m = 1:numel(xModes)
                    i = find(sel & strcmp({S.Animal}, animals{a}) & strcmp({S.Base}, xModes{m}), 1);
                    if ~isempty(i), A(a, m) = S(i).acc; end
                end
            end
            nUsed = sum(any(~isnan(A), 2));
            if nUsed == 0
                text(ax, 0.5, 0.5, 'no data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                continue;
            end

            for a = 1:nA
                plot(ax, 1:numel(xModes), A(a, :), '-o', 'Color', animalColor, 'MarkerFaceColor', animalColor, ...
                     'LineWidth', 1.2, 'MarkerSize', 6, 'HandleVisibility', 'off');
            end
            h = plot(ax, nan, nan, '-o', 'Color', animalColor, 'MarkerFaceColor', animalColor, ...
                     'DisplayName', sprintf('Animals (n=%d)', nUsed));
            if nUsed > 1
                hM = plot(ax, 1:numel(xModes), mean(A, 1, 'omitnan'), '-s', 'Color', 'k', 'MarkerFaceColor', 'k', ...
                          'LineWidth', 3, 'MarkerSize', 9, 'DisplayName', 'Mean');
                h = [h hM];
            end
            if isempty(hLeg), hLeg = h; axLeg = ax; end
        end
    end
    title(tl, 'SVM accuracy (all bands) across analyses');
    if ~isempty(hLeg)
        lgd = legend(axLeg, hLeg);
        lgd.Layout.Tile = 'east';
    end
end

function saveFiguresToPPT(pptFile)
% one slide per open figure (in creation order); slide title = the figure's Name.
% Needs the MATLAB Report Generator (mlreportgen.ppt).
    import mlreportgen.ppt.*
    figs = findall(0, 'Type', 'figure');
    if isempty(figs), fprintf('No figures to save.\n'); return; end
    [~, order] = sort([figs.Number]);
    figs = figs(order);

    ppt = Presentation(pptFile);
    open(ppt);
    tmpDir = tempname;  mkdir(tmpDir);
    for k = 1:numel(figs)
        f = figs(k);
        ttl = f.Name;
        if isempty(ttl), ttl = sprintf('Figure %d', f.Number); end
        png = fullfile(tmpDir, sprintf('fig_%03d.png', k));
        exportgraphics(f, png, 'Resolution', 200);

        slide = add(ppt, 'Title and Picture');
        replace(slide, 'Title', ttl);
        replace(slide, 'Picture', Picture(png));
    end
    close(ppt);
    rmdir(tmpDir, 's');
    fprintf('Saved %d figures to %s\n', numel(figs), fullfile(pwd, pptFile));
end