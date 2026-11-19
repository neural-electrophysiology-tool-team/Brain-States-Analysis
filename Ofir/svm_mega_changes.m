%% ============================================================
%  ONE SVM ON ALL FREQUENCY BANDS - SEPARATELY FOR EACH ANIMAL   (computation only; plotting = replot_svm_mega.m)
%  averaged over random train/test splits, with an optional permutation null
%
%  Main mode: 'category_pairs' (see pairDefs). The older modes are still supported (see modeNeeds).
%
%  CHANGES vs. the previous version
%   1. "other" bugs are split into three separate trial categories (typeMask / otherTypes):
%         lowNoise : low horizontal noise
%         lowRand  : low horizontal + random low horizontal
%         random   : random
%      (+ the strike variants lowNoiseS / lowRandS / randomS). Every category gets its own paired comparisons,
%      its own onset / offset comparisons, and its own strike vs no-strike comparisons.
%   2. PRE-OFFSET window of hunters ends preOffsetEnd (= -2) s from the OFFSET  -> [-2-L, -2) s from OFFSET.
%   3. circle vs any "other" category (unpaired onset / offset, no-strike and strike): windows of winLenShort (3 s).
%   4. paired comparisons of circle and meadow (sets in longWinSets): windows of winLenLong (10 s).
%      Everything else: winLen (5 s).
%
%  Windows for a window length L (bins of binSize s, end point NOT included):
%   PRE         : [preEnd - L, preEnd)       from stimulus ONSET   (preEnd = -1  -> [-6,-1) for L = 5)
%   ONSET       : [onsetStart, +L)           from stimulus ONSET   ([-1, 4) for L = 5)
%   OFFSET      : [offsetStart, +L)          from stimulus OFFSET  ([ 1, 6) for L = 5)
%   PRE-OFFSET  : [preOffsetEnd - L, preOffsetEnd) from OFFSET     ([-7,-2) for L = 5, [-12,-2) for L = 10)
%
%  PAIRED   : two windows of the same trials; trials whose two windows overlap in time are removed from THAT
%             comparison only. A pair stays together in train or test.
%  UNPAIRED : different events; the same number of events is drawn at random from each group.
%
%  zMode 'trial': all windows of a trial (at that comparison's window length L) are z-scored together.
%  => normalization is identical for all comparisons WITH THE SAME window length, but differs between
%     L = 3 / 5 / 10 (the trial's windows have different lengths).
%
%  Model 1 = SVM on all bands (features = nBands x nBins, band-fastest); optional ablation = models 2..6.
%  Permutation null optional (all-bands model only).
%
%  Saved: svm_allbands_rawnorm_<mode>.mat with batchResults (acc, pVal, nullAcc, beta, betaScaled, nPerClass,
%  binCenters, winLen, xLab, design, class0, class1 ...) and groupOrder (order of the figure groups).
% =============================================================

clear; clc;

%% ------------------------------------------------------------
%  User Parameters
% ------------------------------------------------------------
analysisModes = {'category_pairs'};   % can be combined with the old modes, e.g. {'category_pairs','pre_vs_onset'}

hunterTables = {'hunter106', 'hunter126', 'hunter143', 'hunter153', 'hunter157', 'hunter159', 'hunter161','hunter162', 'hunter167', 'hunter208', 'hunter213', 'hunter228'};

% ---- the three "other" categories (same mechanism as the time-course script) ----
% 'types' = values of the column Type of the hunter tables that belong to the group (pooled).
% Leave {} to choose them in a dialog at run time (the unique non-circle Type values are read from the data).
% After the first run the chosen lists are printed in the command window: paste them here to skip the dialogs.
% 'key' is the internal name of the category (do not change).
otherGroups = struct('key',  {'lowNoise', 'lowRand', 'random'}, ...
                     'name', {'OTHER: LOW HORIZONTAL NOISE', ...
                              'OTHER: LOW HORIZONTAL + RANDOM LOW HORIZONTAL', ...
                              'OTHER: RANDOM'}, ...
                     'types', {{}, {}, {}});

minTrialsRequired = 15;
strikeMinTime_s   = 4;   % strike trials: keep only trials with Strike_ms > this many s
minTrialsStrike   = 15;
binSize           = 1;   % s
winLen            = 5;   % s, default window length
winLenShort       = 3;   % s, circle vs any "other" category (unpaired)
winLenLong        = 10;  % s, paired comparisons of the sets in longWinSets
longWinSets       = {'circle', 'meadow'};   % add 'circleS' for circle-strike as well
preEnd            = -1;  % PRE window ends here (s from ONSET)            -> PRE = [preEnd - L, preEnd)
onsetStart        = -1;  % ONSET window starts here (s from ONSET)
offsetStart       = 1;   % OFFSET window starts here (s from OFFSET)
preOffsetEnd      = -2;  % PRE-OFFSET window ends here (s from OFFSET, i.e. 2 s before the end of the trial)
preStart          = preEnd - winLen;   % default-length PRE window (used by the old modes only)
trainFraction     = 0.8;
nSplits           = 1;
nPerm             = 10;
runPermutations   = true;
runAblation       = false;

preOffsetTypes = {'circle'};   % trial categories that have a pre-offset window: any of {'circle','lowNoise','lowRand','random'}

zMode          = 'trial';   % 'recording' | 'window' | 'trial' | 'none'
normToSum      = false;
weightScaleFcn = @(a) a;

% ---- meadow / sound data ----
extraMatFiles  = {};
meadowVarFmt   = 'reg_meadow%d';
meadowState    = 'move';
soundVarName   = 'reward_sound';
soundDur_ms    = 1000;
soundExclWin_s = 5;
soundRecordFilter = @(rec) contains(lower(string(rec)), 'wake') & ~contains(lower(string(rec)), 'hunter');
minTrialsSound = 2;

bandRanges = [0 4; 5 12; 13 30; 31 40; 60 80];
bandNames  = ["Delta (0-4 Hz)", "Theta+Alpha (5-12 Hz)", "Beta (13-30 Hz)", ...
              "GammaLow (31-40 Hz)", "GammaHigh (60-80 Hz)"];
bandShort  = ["Delta", "Theta+Alpha", "Beta", "GammaLow", "GammaHigh"];
nBands   = size(bandRanges, 1);
nBinsWin = round(winLen / binSize);
nPermEff = nPerm * runPermutations;

% window setup shared by all trial sets; windowSetup(Tbase, L) adds the window starts for window length L
Tbase = struct('preEnd', preEnd, 'onsetStart', onsetStart, 'offsetStart', offsetStart, 'preOffsetEnd', preOffsetEnd, ...
               'preOffsetTypes', {preOffsetTypes}, 'strikeMin_s', strikeMinTime_s, 'otherTypes', struct(), ...   % otherTypes filled after loading the data
               'soundDur_ms', soundDur_ms, 'soundRecordFilter', soundRecordFilter, 'verbose', true);
T = windowSetup(Tbase, winLen);          % default window length (old modes)
cfgWin = struct('winLen', winLen, 'winLenShort', winLenShort, 'winLenLong', winLenLong, 'longWinSets', {longWinSets});

loadData = load('/home/ofir/Documents/MATLAB/data4paper.mat');
for k = 1:numel(extraMatFiles)
    ex = load(extraMatFiles{k});
    for fn = fieldnames(ex)', loadData.(fn{1}) = ex.(fn{1}); end
end

%% ---- "other" groups: read the Type values from the data and assign them to the groups ----
otherGroups = resolveOtherGroups(loadData, hunterTables, otherGroups);
otherTypes = struct();
for g = 1:numel(otherGroups)
    otherTypes.(otherGroups(g).key) = string(otherGroups(g).types);   % key -> pooled Type values
end
Tbase.otherTypes = otherTypes;
T = windowSetup(Tbase, winLen);          % rebuild with the resolved groups

recCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
saCache  = containers.Map('KeyType', 'char', 'ValueType', 'any');

%% ============================================================
%  LOOP OVER ANALYSIS MODES
% =============================================================
for mIdx = 1:numel(analysisModes)
    analysisMode = analysisModes{mIdx};
    fprintf('\n===== MODE: %s =====\n', analysisMode);
    isStrike = endsWith(analysisMode, '_strike');
    baseMode = erase(analysisMode, '_strike');
    isPairs  = strcmp(baseMode, 'category_pairs');

    [needBug, needMeadow, needSound, splitBug] = modeNeeds(analysisMode);
    groups = [];  groupOrder = {};
    if isPairs
        pdefs  = pairDefs(preOffsetTypes, cfgWin);
        types  = {pdefs.name};
        groups = pairGroups(pdefs);
        groupOrder = {groups.name};
        fprintf('%d comparisons per animal\n', numel(types));
    elseif splitBug, types = {'circle', 'other'}; else, types = {analysisMode}; end
    soundTbl = [];
    if needSound
        if isfield(loadData, soundVarName), soundTbl = loadData.(soundVarName);
        elseif ~isPairs, error('Sound table "%s" not found in loadData (see soundVarName / extraMatFiles).', soundVarName);
        end
    end

    batchResults = struct('Animal', {}, 'RunLabel', {}, 'TrialType', {}, 'ModelIdx', {}, 'ModelLabel', {}, ...
                          'acc', {}, 'beta', {}, 'betaScaled', {}, 'nullAcc', {}, 'nPerClass', {}, 'pVal', {}, ...
                          'binCenters', {}, 'winLen', {}, 'xLab', {}, 'design', {}, 'class0', {}, 'class1', {});
    runCounter = 0;
    binCenters = []; xLab = '';

    for tIdx = 1:numel(hunterTables)
        hTableName = hunterTables{tIdx};
        animalID = str2double(regexp(hTableName, '\d+', 'match', 'once'));

        currentHunter = [];
        if needBug
            if ~isfield(loadData, hTableName), continue; end
            currentHunter = loadData.(hTableName);
            if isPairs, reportTypes(currentHunter, T, hTableName); end
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
        setsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');   % trial sets, one per window length

        for cType = types
            cType = cType{1};
            desc = sprintf('%s | %s | %s', hTableName, analysisMode, cType);

            minN   = minTrialsRequired;
            if needSound, minN = minTrialsSound; end
            if isStrike, minN = minTrialsStrike; end
            Lc = winLen;  TL = T;  nBins = nBinsWin;
            class0Name = '';  class1Name = '';

            %% ---- build the two groups of events (E0 = class 0, E1 = class 1) ----
            switch baseMode
                case 'pre_vs_onset'
                    E0 = hunterEvents(bugTrials(currentHunter, cType, isStrike, T), 'onset', T, cType);
                    E1 = E0;
                    start0 = preStart; start1 = onsetStart; design = 'paired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim ONSET; pre window = [%g, %g) s, onset window = [%g, %g) s)', ...
                                   preStart, preStart + winLen, onsetStart, onsetStart + winLen);

                case 'pre_vs_offset'
                    tblB = bugTrials(currentHunter, cType, isStrike, T);
                    E0 = hunterEvents(tblB, 'onset', T, cType);
                    E1 = hunterEvents(tblB, 'offset', T, cType);
                    start0 = preStart; start1 = offsetStart; design = 'paired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim OFFSET; offset window = [%g, %g) s; pre window = [%g, %g) s from ONSET)', ...
                                   offsetStart, offsetStart + winLen, preStart, preStart + winLen);

                case 'onset_vs_offset'
                    tblB = bugTrials(currentHunter, cType, isStrike, T);
                    E0 = hunterEvents(tblB, 'onset', T, cType);
                    E1 = hunterEvents(tblB, 'offset', T, cType);
                    start0 = onsetStart; start1 = offsetStart; design = 'paired'; binStart = start0;
                    xLab = sprintf('Time bin (axis = s from stim ONSET [%g, %g); class 1 = same bin of [%g, %g) s from stim OFFSET)', ...
                                   onsetStart, onsetStart + winLen, offsetStart, offsetStart + winLen);

                case 'circle_vs_other'
                    E0 = hunterEvents(bugTrials(currentHunter, 'circle', isStrike, T), 'onset', T, 'circle');
                    E1 = hunterEvents(bugTrials(currentHunter, 'other',  isStrike, T), 'onset', T, 'other');
                    start0 = onsetStart; start1 = onsetStart; design = 'unpaired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim ONSET; window = [%g, %g) s)', onsetStart, onsetStart + winLen);

                case {'bug_vs_meadow_onset', 'bug_vs_meadow_offset'}
                    align = 'onset'; wStart = onsetStart;
                    if endsWith(analysisMode, 'offset'), align = 'offset'; wStart = offsetStart; end
                    E0 = hunterEvents(bugSubset(currentHunter, cType, T), align, T, cType);
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
                    E0 = hunterEvents(bugSubset(currentHunter, cType, T), align, T, cType);
                    E1 = soundEvents(soundTbl, animalID, soundRecordFilter, align, soundDur_ms, soundExclWin_s, T);
                    start0 = wStart; start1 = wStart; design = 'unpaired'; binStart = start1;
                    xLab = sprintf('Time bin (s from stim %s; window = [%g, %g) s)', upper(align), wStart, wStart + winLen);

                case 'category_pairs'
                    d  = pdefs(strcmp({pdefs.name}, cType));
                    Lc = d.winLen;                                  % window length of THIS comparison
                    TL = windowSetup(Tbase, Lc);
                    nBins = round(Lc / binSize);
                    skey = sprintf('L%g', Lc);
                    if ~isKey(setsMap, skey)                        % trial sets built once per animal and window length
                        TLv = TL;  TLv.verbose = (setsMap.Count == 0);
                        setsMap(skey) = buildTrialSets(currentHunter, meadowTbl, soundTbl, animalID, TLv);
                    end
                    [E0, E1, start0, start1, design, nOverlap] = buildPair(d, setsMap(skey), TL);
                    minN = min(setMinTrials(d.setA, minTrialsRequired, minTrialsStrike, minTrialsSound), ...
                               setMinTrials(d.setB, minTrialsRequired, minTrialsStrike, minTrialsSound));
                    binStart = start1;
                    xLab = sprintf('Time bin (s from window start; %g s windows; class 0 = %s, class 1 = %s)', Lc, d.labelA, d.labelB);
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
                fprintf('Processing: %s [%g s] (valid %d / %d -> %d random events from each)\n', desc, Lc, numel(ok0), numel(ok1), nUse);
                rng(1);
                c0All = F0(ok0(randperm(numel(ok0), nUse)), :);
                c1All = F1(ok1(randperm(numel(ok1), nUse)), :);
                nPerClass = nUse;
            else
                c0All = F0; c1All = F1;
                nPerClass = sum(~any(isnan(F0), 2) & ~any(isnan(F1), 2));
                fprintf('Processing: %s [%g s] (N=%d pairs, %d complete)\n', desc, Lc, size(F0, 1), nPerClass);
            end

            %% ---- SVM on all bands (+ optionally 5 leave-one-band-out SVMs) ----
            if runAblation, nModels = 1 + nBands; else, nModels = 1; end
            for m = 1:nModels
                if m == 1
                    keepBands = 1:nBands;  modelLabel = "All bands";
                else
                    keepBands = setdiff(1:nBands, m-1);  modelLabel = "w/o " + bandShort(m-1);
                end
                colMask = repmat(ismember(1:nBands, keepBands), 1, nBins);   % column = (bin-1)*nBands + band

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
                batchResults(runCounter).betaScaled = bFull * weightScaleFcn(accAll);
                batchResults(runCounter).pVal       = pVal;
                batchResults(runCounter).binCenters = binCenters;
                batchResults(runCounter).winLen     = Lc;
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
            fprintf('Verification: ablation ON -> %d SVMs per dataset, %d SVMs in total.\n', 1 + nBands, numel(batchResults));
        else
            assert(all([batchResults.ModelIdx] == 1), 'Verification failed: extra SVMs were run although runAblation = false.');
            fprintf('Verification: ablation OFF -> only the all-bands SVM was run (%d SVMs in total).\n', numel(batchResults));
        end
    end

    %% Save (plot with replot_svm_mega.m)
    if ~isempty(batchResults)
        save(sprintf('svm_allbands_rawnorm_%s.mat', analysisMode), 'batchResults', 'bandNames', ...
             'binCenters', 'analysisMode', 'winLen', 'winLenShort', 'winLenLong', 'longWinSets', ...
             'preEnd', 'preStart', 'onsetStart', 'offsetStart', 'preOffsetEnd', 'otherTypes', 'otherGroups', 'groupOrder', ...
             'nSplits', 'nPermEff', 'runPermutations', 'runAblation', 'normToSum', 'zMode', 'meadowState', ...
             'trainFraction', 'minTrialsRequired', 'minTrialsStrike', 'minTrialsSound', 'strikeMinTime_s', ...
             'preOffsetTypes', 'bandRanges');
        fprintf('Saved svm_allbands_rawnorm_%s.mat\n', analysisMode);
    end
end

%% ============================================================
%  HELPER FUNCTIONS
% =============================================================

function T = windowSetup(Tb, L)
% window starts (s relative to the reference event) for window length L
    T = Tb;
    T.winLen   = L;
    T.preStart = Tb.preEnd - L;                                   % PRE        = [preEnd - L, preEnd) from ONSET
    preOffStart = Tb.preOffsetEnd - L;                            % PRE-OFFSET = [preOffsetEnd - L, preOffsetEnd) from OFFSET
    T.hunterStarts = [T.preStart, Tb.onsetStart, preOffStart, Tb.offsetStart];   % 1 pre, 2 onset, 3 pre-offset, 4 offset
    T.meadowStarts = [T.preStart, Tb.onsetStart, T.preStart, Tb.onsetStart];     % 1 pre-move, 2 onset-move, 3 pre-stat, 4 onset-stat
    T.soundStarts  = [T.preStart, Tb.onsetStart, Tb.offsetStart];                % 1 pre, 2 onset, 3 offset
end


function [needBug, needMeadow, needSound, splitBug] = modeNeeds(mode)
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
    splitBug   = needBug && ~ismember(base, {'circle_vs_other', 'category_pairs'});
end


function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end


function n = setMinTrials(s, nStd, nStrike, nSound)
% minimum number of trials of a trial set (strike sets end with 'S')
    if strcmp(s, 'sound'), n = nSound;
    elseif numel(s) > 1 && s(end) == 'S', n = nStrike;
    else, n = nStd; end
end


%% ---------------- trial categories ----------------

function groups = resolveOtherGroups(loadData, hunterTables, groups)
% Reads the Type column of all hunter tables, lists the non-circle types (with trial counts per animal),
% and lets the user assign types to the groups whose 'types' is empty (dialog). Preset types are checked
% against the data. The final assignment is printed so it can be pasted into otherGroups.
    allT = strings(0, 1);  info = {};
    for a = 1:numel(hunterTables)
        if ~isfield(loadData, hunterTables{a}), continue; end
        Ty = string(loadData.(hunterTables{a}).Type);
        Ty = Ty(~(Ty == "circle" | Ty == "circle_flip") & ~ismissing(Ty));
        [u, ~, ic] = unique(Ty);
        cnt = accumarray(ic, 1);
        info{end+1} = sprintf('  %s: %s', hunterTables{a}, ...
                              strjoin(arrayfun(@(i) sprintf('%s (%d)', u(i), cnt(i)), 1:numel(u), 'UniformOutput', false), ', ')); %#ok<AGROW>
        allT = [allT; u(:)]; %#ok<AGROW>
    end
    allT = unique(allT);
    fprintf('\nNon-circle values of Type in the hunter tables (trials per animal):\n');
    fprintf('%s\n', info{:});
    if isempty(allT), error('No non-circle Type values found in the hunter tables.'); end

    for g = 1:numel(groups)
        if isempty(groups(g).types)
            [sel, ok] = listdlg('PromptString', {['Group: ' groups(g).name], 'Select the Type values that belong to it:'}, ...
                                'ListString', cellstr(allT), 'SelectionMode', 'multiple', 'ListSize', [420 260], ...
                                'Name', 'Assign Type values');
            if ~ok || isempty(sel), error('No Type selected for group "%s".', groups(g).name); end
            groups(g).types = cellstr(allT(sel))';
        else
            bad = setdiff(string(groups(g).types), allT);
            if ~isempty(bad)
                error('Group "%s": Type value(s) not found in the data: %s', groups(g).name, strjoin(bad, ', '));
            end
        end
    end
    fprintf('\nGroups used (paste into otherGroups to skip the dialogs):\n');
    for g = 1:numel(groups)
        fprintf('  %-48s types = {%s}\n', groups(g).name, strjoin("'" + string(groups(g).types) + "'", ','));
    end
    fprintf('\n');
end


function m = typeMask(H, cType, T)
% logical column: trials of category cType
%   'circle'  : circle, circle_flip
%   'other'   : ALL non-circle trials (legacy, old modes)
%   'lowNoise' / 'lowRand' / 'random' : H.Type in T.otherTypes.(cType)
    ty = string(H.Type);
    isCircle = (ty == "circle" | ty == "circle_flip");
    switch cType
        case 'circle', m = isCircle;
        case 'other',  m = ~isCircle;
        otherwise
            if ~isfield(T.otherTypes, cType), error('Unknown trial category: %s', cType); end
            m = ismember(ty, string(T.otherTypes.(cType)));
    end
    m = m(:);
end


function reportTypes(H, T, name)
% prints the Type values found for an animal and how they are assigned to the categories
    ty = string(H.Type);  ty = ty(:);
    u = unique(ty(~ismissing(ty)));
    fprintf('%s: Type values = %s\n', name, strjoin(u', ', '));
    assigned = (ty == "circle" | ty == "circle_flip");
    fprintf('   circle: %d trials', nnz(assigned));
    for c = fieldnames(T.otherTypes)'
        mk = ismember(ty, string(T.otherTypes.(c{1})));
        assigned = assigned | mk;
        fprintf(' | %s: %d', c{1}, nnz(mk));
    end
    fprintf('\n');
    un = unique(ty(~assigned & ~ismissing(ty)));
    if ~isempty(un)
        fprintf('   NOT assigned to any category (ignored): %s\n', strjoin(un', ', '));
    end
end


function tbl = bugSubset(H, cType, T)
% trials of category cType WITHOUT a strike
    tbl = H(typeMask(H, cType, T) & isnan(H.Strike_ms), :);
end


function tbl = bugStrikeSubset(H, cType, T)
% trials of category cType WITH a strike and Strike_ms > T.strikeMin_s seconds
    tbl = H(typeMask(H, cType, T) & H.Strike_ms > T.strikeMin_s * 1000, :);
end


function tbl = bugTrials(H, cType, isStrike, T)
    if ~isStrike, tbl = bugSubset(H, cType, T); return; end
    tbl = bugStrikeSubset(H, cType, T);
    if isfield(T, 'verbose') && T.verbose
        nAll = sum(typeMask(H, cType, T) & ~isnan(H.Strike_ms));
        fprintf('   strike trials (%s): %d total, %d kept (Strike_ms > %g s)\n', cType, nAll, height(tbl), T.strikeMin_s);
    end
end


%% ---------------- event tables (with trial info for zMode 'trial') ----------------

function E = hunterSet(tbl, T, cType)
% one row per hunter trial; refs = [onset onset offset offset], wStarts = T.hunterStarts
% pre-offset column = NaN for categories not in T.preOffsetTypes
    n   = height(tbl);
    on  = tbl.Triger_ms(:);  off = tbl.End_ms(:);
    refs = [on on off off];
    wSt  = repmat(T.hunterStarts, n, 1);
    if ~ismember(cType, T.preOffsetTypes)
        refs(:, 3) = NaN;  wSt(:, 3) = NaN;
    end
    E = table(string(tbl.Animal), string(tbl.Record), nan(n, 1), refs, wSt, ...
              'VariableNames', {'Animal', 'Record', 't_ms', 'refs', 'wStarts'});
end


function E = hunterEvents(tbl, align, T, cType)
    E = hunterSet(tbl, T, cType);
    if strcmp(align, 'offset'), E.t_ms = tbl.End_ms(:); else, E.t_ms = tbl.Triger_ms(:); end
end


function E = meadowSet(M, animalID, T)
    n = height(M);
    if ismember('Animal', M.Properties.VariableNames), an = string(M.Animal);
    else, an = repmat("PV" + animalID, n, 1); end
    statOn = M.Triger_ms(:);  moveOn = M.Move_ms(:);
    E = table(an(:), string(M.Record), nan(n, 1), [moveOn moveOn statOn statOn], repmat(T.meadowStarts, n, 1), ...
              'VariableNames', {'Animal', 'Record', 't_ms', 'refs', 'wStarts'});
end


function E = meadowEvents(M, animalID, evtName, T)
% evtName: 'statOn' | 'statOff' | 'moveOn' | 'moveOff'
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
% all trial sets of ONE animal for window setup T:
%   circle / lowNoise / lowRand / random (no strike), the same + 'S' (strike), meadow, sound
    e4 = struct('E', emptySetTable(4), 'starts', T.hunterStarts);
    sets = struct('meadow', struct('E', emptySetTable(4), 'starts', T.meadowStarts), ...
                  'sound',  struct('E', emptySetTable(3), 'starts', T.soundStarts));
    for c = {'circle', 'lowNoise', 'lowRand', 'random'}
        cn = c{1};
        sets.(cn)        = e4;
        sets.([cn 'S'])  = e4;
        if ~isempty(H)
            sets.(cn).E       = hunterSet(bugTrials(H, cn, false, T), T, cn);
            sets.([cn 'S']).E = hunterSet(bugTrials(H, cn, true,  T), T, cn);
        end
    end
    if ~isempty(M),   sets.meadow.E = meadowSet(M, animalID, T); end
    if ~isempty(Snd), sets.sound.E  = soundSet(Snd, animalID, T); end
end


%% ---------------- category_pairs: definitions of all comparisons ----------------

function D = pairDefs(preOffsetTypes, cfg)
% all comparisons; every comparison carries its own window length D.winLen (see pairLen)
    D = struct('group', {}, 'name', {}, 'setA', {}, 'winA', {}, 'labelA', {}, 'setB', {}, 'winB', {}, 'labelB', {}, 'winLen', {});
    labs4 = {'Baseline', 'Onset', 'Pre-offset', 'Offset'};   labs3 = {'Baseline', 'Onset', 'Offset'};

    % bug categories: set key | label of the bug (unpaired figures) | name of the paired group
    bugs = {'circle',   'Circle bug',         'circle'
            'lowNoise', 'Low-H noise bug',    'low-H noise'
            'lowRand',  'Low-H + rand low-H bug', 'low-H + rand low-H'
            'random',   'Random bug',         'random'};

    paired = cell(0, 4);
    for i = 1:size(bugs, 1)
        key = bugs{i, 1};  nm = bugs{i, 3};
        if ismember(key, preOffsetTypes), lb = labs4; wi = [1 2 3 4]; else, lb = labs3; wi = [1 2 4]; end
        paired(end+1, :) = {nm,             key,        lb, wi}; %#ok<AGROW>
        paired(end+1, :) = {[nm ' strike'], [key 'S'],  lb, wi}; %#ok<AGROW>
    end
    paired(end+1, :) = {'meadow', 'meadow', {'Pre move', 'Onset move', 'Pre stat', 'Onset stat'}, [1 2 3 4]};
    paired(end+1, :) = {'sound',  'sound',  labs3, [1 2 3]};
    for i = 1:size(paired, 1)
        lb = paired{i, 3};  wi = paired{i, 4};
        for a = 1:numel(wi) - 1
            for b = a + 1:numel(wi)
                D(end+1) = mkDef(paired{i, 1}, paired{i, 2}, wi(a), lb{a}, paired{i, 2}, wi(b), lb{b}, cfg); %#ok<AGROW>
            end
        end
    end

    bs = bugs(:, 1)';  bl = bugs(:, 2)';  bsS = strcat(bs, 'S');
    unp = {'onset',         [bs {'meadow', 'sound'}], [2 2 2 2 2 2], [bl {'Meadow move', 'Sound'}]
           'offset',        [bs {'sound'}],           [4 4 4 4 3],   [bl {'Sound'}]
           'onset strike',  bsS,                      [2 2 2 2],     bl
           'offset strike', bsS,                      [4 4 4 4],     bl};
    for i = 1:size(unp, 1)
        S = unp{i, 2};  W = unp{i, 3};  L = unp{i, 4};
        for a = 1:numel(S) - 1
            for b = a + 1:numel(S)
                D(end+1) = mkDef(unp{i, 1}, S{a}, W(a), L{a}, S{b}, W(b), L{b}, cfg); %#ok<AGROW>
            end
        end
    end

    % strike vs no-strike, separately for every bug category and for onset / offset (window 2 = onset, 4 = offset)
    for i = 1:size(bugs, 1)
        key = bugs{i, 1};  nm = bugs{i, 3};
        D(end+1) = mkDef([nm ' strike vs no-strike'], key, 2, 'No-strike Onset',  [key 'S'], 2, 'Strike Onset',  cfg); %#ok<AGROW>
        D(end+1) = mkDef([nm ' strike vs no-strike'], key, 4, 'No-strike Offset', [key 'S'], 4, 'Strike Offset', cfg); %#ok<AGROW>
    end
end


function L = pairLen(setA, setB, cfg)
% window length of a comparison
%   same set (paired) in cfg.longWinSets             -> winLenLong  (circle, meadow)
%   circle (or circle strike) vs an "other" category -> winLenShort
%   everything else                                  -> winLen
    circ = {'circle', 'circleS'};
    oth  = {'lowNoise', 'lowRand', 'random', 'lowNoiseS', 'lowRandS', 'randomS'};
    if strcmp(setA, setB)
        if ismember(setA, cfg.longWinSets), L = cfg.winLenLong; else, L = cfg.winLen; end
    elseif (ismember(setA, circ) && ismember(setB, oth)) || (ismember(setB, circ) && ismember(setA, oth))
        L = cfg.winLenShort;
    else
        L = cfg.winLen;
    end
end


function d = mkDef(group, setA, winA, labA, setB, winB, labB, cfg)
    d = struct('group', group, 'name', sprintf('%s | %s vs %s', group, labA, labB), ...
               'setA', setA, 'winA', winA, 'labelA', labA, 'setB', setB, 'winB', winB, 'labelB', labB, ...
               'winLen', pairLen(setA, setB, cfg));
end


function G = pairGroups(D)
    [~, ia] = unique({D.group}, 'stable');
    names = {D(ia).group};
    G = struct('name', names, 'types', cell(1, numel(names)));
    for g = 1:numel(names)
        G(g).types = {D(strcmp({D.group}, names{g})).name};
    end
end


function [E0, E1, start0, start1, design, nOverlap] = buildPair(d, sets, T)
% Same set -> PAIRED (trials whose two windows overlap are removed); two sets -> UNPAIRED.
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
    a1 = r1 + s1 * 1000;  a2 = r2 + s2 * 1000;
    ov = a1 < a2 + winLen * 1000 & a2 < a1 + winLen * 1000;
end


%% ---------------- zMode 'trial' helpers ----------------

function tf = hasTrialRef(E)
    tf = all(ismember({'refs', 'wStarts'}, E.Properties.VariableNames));
end


function [mu, sd, ok] = trialStats(R, refs, starts, winLen)
% mean / std per band over ALL windows of ONE trial (union of their time points).
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
    D = extractLFPBins(E, X.SA, X.bandRanges, spec, X.binSize, X.recCache, X.normToSum, X.zMode);
    F = reshape(D, size(D, 1) * size(D, 2), [])';
end


function [F0, F1] = featuresPaired(E0, spec0, E1, spec1, X)
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
        both = [seg0, seg1];
        mu = mean(both, 2, 'omitnan');
        sd = std(both, 0, 2, 'omitnan');
        sd(sd == 0) = NaN;

        D0(:, :, tr) = binSeg((seg0 - mu) ./ sd, rel0(w0), edges0);
        D1(:, :, tr) = binSeg((seg1 - mu) ./ sd, rel1(w1), edges1);
    end

    F0 = reshape(D0, nBands * nB0, [])';
    F1 = reshape(D1, nBands * nB1, [])';
end


function out = binSeg(seg, sub, tEdges)
    nBins = numel(tEdges) - 1;
    out = nan(size(seg, 1), nBins);
    for tb = 1:nBins
        idx = sub >= tEdges(tb) & sub < tEdges(tb+1);
        if any(idx)
            out(:, tb) = mean(seg(:, idx), 2, 'omitnan');
        end
    end
end


function R = getRec(recCache, SA, animal, record, bandRanges, normToSum)
    key = [animal '|' record];
    if ~isKey(recCache, key)
        recCache(key) = loadBandZ(SA, animal, record, bandRanges, normToSum);
    end
    R = recCache(key);
end


%% ---------------- SVM ----------------

function [acc, beta, nullAcc] = trainTestSVM(c0Mat, c1Mat, trainFraction, nSplits, nPerm, minTrials)
% Paired linear SVM; row i of c0Mat and row i of c1Mat form a pair that always stays together in train or test.
    acc     = nan;
    beta    = nan(1, size(c0Mat, 2));
    nullAcc = nan(nPerm, 1);

    good = find(~any(isnan(c0Mat), 2) & ~any(isnan(c1Mat), 2));
    if numel(good) < minTrials, return; end
    X0 = c0Mat(good, :);
    X1 = c1Mat(good, :);
    nPairs = numel(good);

    [acc, beta] = svmOverSplits(X0, X1, false(nPairs, 1), trainFraction, nSplits);

    for p = 1:nPerm
        swap = rand(nPairs,1) < 0.5;
        nullAcc(p) = svmOverSplits(X0, X1, swap, trainFraction, nSplits);
    end
end


function [acc, beta] = svmOverSplits(X0, X1, swap, trainFraction, nSplits)
    nPairs = size(X0, 1);
    nTrain = floor(trainFraction * nPairs);
    y0 = double(swap);
    y1 = double(~swap);

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

function binnedData = extractLFPBins(evt, SA, bandRanges, spec, binSize, recCache, normToSum, zMode)
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
        if isempty(R), continue; end

        relTime = (R.t_ms - alignTime_ms) / 1000;
        w = find(relTime >= tEdges(1) & relTime < tEdges(end));
        if isempty(w), continue; end
        sub = relTime(w);

        switch zMode
            case 'window'
                seg = R.bandRaw(:, w);
                sd  = std(seg, 0, 2, 'omitnan');
                sd(sd == 0) = NaN;
                seg = (seg - mean(seg, 2, 'omitnan')) ./ sd;
            case 'trial'
                if hasTrial
                    % stats over ALL windows of this trial; window length = nBins*binSize of this comparison
                    [mu, sd, ok] = trialStats(R, evt.refs(tr, :), evt.wStarts(tr, :), nBins * binSize);
                    if ~ok, continue; end
                    seg = (R.bandRaw(:, w) - mu) ./ sd;
                else
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
    if normToSum
        B = B ./ sum(B, 1, 'omitnan');
        B(~isfinite(B)) = NaN;
    end

    R.bandRaw = B;
    R.bandZ   = (B - mean(B, 2, 'omitnan')) ./ std(B, 0, 2, 'omitnan');
    R.t_ms    = S.t_ms(:)';
end