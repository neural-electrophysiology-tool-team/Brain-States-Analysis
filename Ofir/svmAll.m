%% ============================================================
%  MULTI-ANIMAL SVM CLASSIFICATION (FLEXIBLE ANALYSIS MODES)
%  - Bin size: 1.0 second
%  - Figure A: Individual accuracies per animal run
%  - Figure B: plotShifted per-run timecourses across 5 frequency bands
%  - Saves plotted variables to svm_plot_variables.mat
% =============================================================

clear; clc;

%% ------------------------------------------------------------
%  User Parameters & Table List
% ------------------------------------------------------------
% Select Analysis Mode: 
% 'pre_vs_post' | 'circle_vs_other' | 'bug_vs_meadow' | 'move_vs_stat' | 'bug_vs_sound'
analysisMode = 'pre_vs_post';

hunterTables = {'hunter106'};%, 'hunter126', 'hunter143', 'hunter153', 'hunter157', 'hunter159', 'hunter161','hunter162', 'hunter167', 'hunter208', 'hunter213', 'hunter228'};
meadowTables = {'reg_meadow106', 'reg_meadow153', 'reg_meadow167', 'reg_meadow208', 'reg_meadow143', 'reg_meadow228'};
soundTables  = {'wake106', 'wake153', 'wake167', 'wake208', 'wake213', 'wake228'};

minTrialsRequired = 15;   % Minimum 15 trials per group
binSize           = 1.0;  % 1.0-sec time bins
trainFraction     = 0.8;  % 80% Train, 20% Test
runBatchPlots     = true; % Generate multi-run figures A & B

% Define the 5 specific frequency bands
bandRanges = [
     0,  4;   % 1. Delta
     5, 12;   % 2. Theta + Alpha
    13, 30;   % 3. Beta
    31, 40;   % 4. GammaLow
    60, 80    % 5. GammaHigh
];
bandNames = ["Delta (0-4 Hz)", "Theta+Alpha (5-12 Hz)", "Beta (13-30 Hz)", ...
             "GammaLow (31-40 Hz)", "GammaHigh (60-80 Hz)"];
nBands = size(bandRanges, 1);

% Load data container
loadData = load('/home/ofir/Documents/MATLAB/data4paper.mat');
recordsNotSync = loadData.recordsNotSync;

% Storage for individual run outputs
batchResults = struct('RunLabel', {}, 'Accuracies', {}, 'Weights', {});
runCounter = 0;


%% ============================================================
%  MAIN PROCESSING LOOP ACCORDING TO ANALYSIS MODE
% =============================================================
for tIdx = 1:numel(hunterTables)
    hTableName = hunterTables{tIdx};
    if ~isfield(loadData, hTableName), continue; end
    
    currentHunter = loadData.(hTableName);
    animalID = str2double(regexp(hTableName, '\d+', 'match', 'once'));
    
    % Select appropriate wakeAnalysis file
    if ismember(animalID, [208, 213, 228, 167])
        excelPath = '/media/sil3/Data/Ofir/BrainStatesExperimentsOfir.xlsx';
    else
        excelPath = '/media/sil1/Data/Nitzan/Experiments/brainStatesWakeN.xlsx';
    end
    SA = wakeAnalysis(excelPath);
    
    trialTypesToRun = {'circle', 'other'};
    
    for typeIdx = 1:numel(trialTypesToRun)
        cType = trialTypesToRun{typeIdx};
        
        % Set up trial masks and data tables based on mode
        switch lower(analysisMode)
            case 'pre_vs_post'
                % 1s-15s Post vs -16s to -2s Pre
                tEdges0 = -6:binSize:-2;
                tEdges1 = 1:binSize:5;
                
                isCircle = (currentHunter.Type == "circle" | currentHunter.Type == "circle_flip");
                isRest   = isnan(currentHunter.Strike_ms);
                if strcmp(cType, 'circle'), selMask = isCircle & isRest;
                else,                       selMask = ~isCircle & isRest; end
                
                table0 = currentHunter(selMask, :);
                table1 = currentHunter(selMask, :);
                
            case 'circle_vs_other'
                % Circle vs Other (1s to 20s)
                tEdges0 = 1:binSize:20;
                tEdges1 = 1:binSize:20;
                
                isRest = isnan(currentHunter.Strike_ms);
                isCircle = (currentHunter.Type == "circle" | currentHunter.Type == "circle_flip");
                
                table0 = currentHunter(~isCircle & isRest, :); % Class 0: Other
                table1 = currentHunter(isCircle & isRest, :);  % Class 1: Circle
                
            case 'bug_vs_meadow'
                % Bug (Hunter) vs Meadow (1s to 15s)
                tEdges0 = 1:binSize:15;
                tEdges1 = 1:binSize:15;
                
                mTableName = sprintf('meadow%d', animalID);
                if ~isfield(loadData, mTableName), continue; end
                currentMeadow = loadData.(mTableName);
                
                isRestH = isnan(currentHunter.Strike_ms);
                isRestM = isnan(currentMeadow.Strike_ms);
                isCircH = (currentHunter.Type == "circle" | currentHunter.Type == "circle_flip");
                isCircM = (currentMeadow.Type == "circle" | currentMeadow.Type == "circle_flip");
                
                if strcmp(cType, 'circle')
                    table0 = currentMeadow(isCircM & isRestM, :);
                    table1 = currentHunter(isCircH & isRestH, :);
                else
                    table0 = currentMeadow(~isCircM & isRestM, :);
                    table1 = currentHunter(~isCircH & isRestH, :);
                end
                
            case 'move_vs_stat'
                % Moving vs Stationary (1s to 15s)
                tEdges0 = 1:binSize:15;
                tEdges1 = 1:binSize:15;
                
                isRest = isnan(currentHunter.Strike_ms);
                isCirc = (currentHunter.Type == "circle" | currentHunter.Type == "circle_flip");
                if strcmp(cType, 'circle')
                    typeMask = isCirc & isRest;
                else
                    typeMask = ~isCirc & isRest;
                end
                
                % Split by locomotion state
                table0 = currentHunter(typeMask & currentHunter.IsMoving == 0, :); % Stationary
                table1 = currentHunter(typeMask & currentHunter.IsMoving == 1, :); % Moving
                
            case 'bug_vs_sound'
                % Bug (Hunter) vs Sound (Wake recordings only, 1s to 15s)
                tEdges0 = 1:binSize:15;
                tEdges1 = 1:binSize:15;
                
                sTableName = sprintf('wake%d', animalID);
                if ~isfield(loadData, sTableName), continue; end
                currentSound = loadData.(sTableName);
                
                isRestH = isnan(currentHunter.Strike_ms);
                isRestS = isnan(currentSound.Strike_ms);
                isCircH = (currentHunter.Type == "circle" | currentHunter.Type == "circle_flip");
                isCircS = (currentSound.Type == "circle" | currentSound.Type == "circle_flip");
                
                if strcmp(cType, 'circle')
                    table0 = currentSound(isCircS & isRestS, :);  % Wake Sound
                    table1 = currentHunter(isCircH & isRestH, :); % Hunter Bug
                else
                    table0 = currentSound(~isCircS & isRestS, :);
                    table1 = currentHunter(~isCircH & isRestH, :);
                end
        end
        
        nTrials0 = height(table0);
        nTrials1 = height(table1);
        
        if nTrials0 < minTrialsRequired || nTrials1 < minTrialsRequired
            fprintf('Skipping %s (%s, %s): Insufficient trials (0: %d, 1: %d)\n', ...
                    hTableName, analysisMode, cType, nTrials0, nTrials1);
            continue;
        end
        
        fprintf('Processing: %s | Mode: %s | Type: %s (N0=%d, N1=%d)\n', ...
                hTableName, analysisMode, cType, nTrials0, nTrials1);
        
        nBins = length(tEdges1) - 1;
        binCenters = tEdges1(1:end-1) + binSize/2;
        
        % Extract LFP Binned Data for both classes
        class0Data = extractLFPBins(table0, SA, recordsNotSync, bandRanges, tEdges0, binSize, strcmp(analysisMode, 'pre_vs_post'));
        class1Data = extractLFPBins(table1, SA, recordsNotSync, bandRanges, tEdges1, binSize, false);
        
        %% ----------------------------------------------------
        %  Train Linear SVM Models per Band
        % ----------------------------------------------------
        rng(1);
        accList = nan(nBands, 1);
        weightMatrix = nan(nBands, nBins);
        
        for fb = 1:nBands
            c0Mat = squeeze(class0Data(fb, :, :))'; 
            c1Mat = squeeze(class1Data(fb, :, :))'; 
            
            good0 = find(~any(isnan(c0Mat), 2));
            good1 = find(~any(isnan(c1Mat), 2));
            
            if numel(good0) < minTrialsRequired || numel(good1) < minTrialsRequired
                continue;
            end
            
            nTrain0 = floor(trainFraction * numel(good0));
            nTrain1 = floor(trainFraction * numel(good1));
            
            shuff0 = good0(randperm(numel(good0)));
            shuff1 = good1(randperm(numel(good1)));
            
            train0 = shuff0(1:nTrain0); test0 = shuff0(nTrain0+1:end);
            train1 = shuff1(1:nTrain1); test1 = shuff1(nTrain1+1:end);
            
            Xtrain = [c0Mat(train0, :); c1Mat(train1, :)];
            yTrain = [zeros(numel(train0),1); ones(numel(train1),1)];
            
            Xtest  = [c0Mat(test0, :);  c1Mat(test1, :)];
            yTest  = [zeros(numel(test0),1); ones(numel(test1),1)];
            
            model = fitcsvm(Xtrain, yTrain, 'KernelFunction', 'linear');
            
            preds = predict(model, Xtest);
            accList(fb) = mean(preds == yTest);
            weightMatrix(fb, :) = model.Beta';
        end
        
        % Store individual run results
        runCounter = runCounter + 1;
        batchResults(runCounter).RunLabel   = sprintf('%s (%s)', hTableName, upper(cType));
        batchResults(runCounter).Accuracies = accList;
        batchResults(runCounter).Weights    = weightMatrix;
    end
end


%% ============================================================
%  EXTRACT AND SAVE PLOTTED VARIABLES FOR FIGURES A & B
% =============================================================
if runBatchPlots && ~isempty(batchResults)
    nRuns = numel(batchResults);
    runColors = lines(nRuns);
    vertShift = 1.5;
    shiftOffsets = 0:vertShift:(vertShift * (nBands - 1));
    
    % Preallocate output matrices
    figA_Accuracies  = nan(nRuns, nBands);             % [Runs x Bands]
    figB_RawWeights  = nan(nBins, nBands, nRuns);       % [Bins x Bands x Runs]
    figB_AbsWeights  = nan(nBins, nBands, nRuns);       % [Bins x Bands x Runs]
    figB_ShiftedData = nan(nBins, nBands, nRuns);       % Shifted values plotted on Y-axis
    runLabels        = cell(nRuns, 1);
    
    for r = 1:nRuns
        runLabels{r} = batchResults(r).RunLabel;
        figA_Accuracies(r, :) = batchResults(r).Accuracies';
        
        rawW = batchResults(r).Weights';                % [nBins x nBands]
        absW = abs(rawW);
        
        figB_RawWeights(:, :, r)  = rawW;
        figB_AbsWeights(:, :, r)  = absW;
        figB_ShiftedData(:, :, r) = bsxfun(@plus, absW, shiftOffsets);
    end
    
    % Save variables to disk
    save('svm_plot_variables_option2_part1.mat', ...
         'figA_Accuracies', 'figB_RawWeights', 'figB_AbsWeights', ...
         'figB_ShiftedData', 'runLabels', 'bandNames', 'binCenters', ...
         'vertShift', 'shiftOffsets', 'analysisMode');
    
    
    %% --- FIGURE A: Accuracy Connected Lines per Run ---
    figure('Name', sprintf('SVM Accuracy (%s)', analysisMode), 'Color', 'w');
    hold on;
    for r = 1:nRuns
        plot(1:nBands, figA_Accuracies(r, :), '-o', 'LineWidth', 1.8, ...
             'Color', runColors(r, :), 'MarkerFaceColor', runColors(r, :), ...
             'DisplayName', runLabels{r});
    end
    yline(0.5, 'r--', 'Chance Level (50%)', 'LineWidth', 1.5, 'DisplayName', 'Chance');
    xlim([0.5, nBands + 0.5]); ylim([0, 1]);
    xticks(1:nBands); xticklabels(bandNames);
    ylabel('Classification Accuracy'); grid on; legend('Location', 'eastoutside');
    title(sprintf('SVM Classification Accuracy (%s)', upper(analysisMode)));
    hold off;
    
    
    %% --- FIGURE B: Weight Magnitudes via plotShifted ---
    figure('Name', sprintf('SVM Weight Magnitudes (%s)', analysisMode), 'Color', 'w');
    hold on;
    for r = 1:nRuns
        plotShifted(binCenters, figB_AbsWeights(:, :, r), '-o', ...
                    'verticalShift', vertShift, 'LineWidth', 1.5, 'MarkerSize', 3, ...
                    'Color', runColors(r, :), 'DisplayName', runLabels{r});
    end
    yticks(shiftOffsets); yticklabels(bandNames);
    xlabel('Time Bins Relative to Stimulus Onset (s)');
    ylabel('Frequency Band (Shifted | \beta | Magnitude)');
    title(sprintf('SVM Weight Magnitudes over Time (%s)', upper(analysisMode)));
    legend('Location', 'eastoutside'); grid on;
    hold off;
end


%% ============================================================
%  HELPER FUNCTIONS
% =============================================================

function binnedData = extractLFPBins(dataTable, SA, recordsNotSync, bandRanges, tEdges, binSize, isPre)
    nBands = size(bandRanges, 1);
    nBins = length(tEdges) - 1;
    nTrials = height(dataTable);
    binnedData = nan(nBands, nBins, nTrials);
    
    currentAnimal = ''; currentRecord = '';
    
    for tr = 1:nTrials
        animal = char(string(dataTable.Animal{tr}));
        record = char(string(dataTable.Record{tr}));
        
        if any(strcmp(recordsNotSync.Animal, animal) & strcmp(recordsNotSync.Record, record))
            continue;
        end
        
        if ~strcmp(animal, currentAnimal) || ~strcmp(record, currentRecord)

            n = find(strcmp(string(SA.recTable.Animal), string(animal)) & strcmp(string(SA.recTable.recNames), string(record)));
            currentAnalysisFolder = fullfile(SA.recTable.folder{n}, 'analysis', ['Animal=' char(animal) ',recNames=' char(record)]);
            ch = SA.recTable.defaulLFPCh(n); 

            % SA.setCurrentRecording(['Animal=' animal ',recNames=' record]);
            % ch = SA.recTable.defaulLFPCh(SA.currentPRec); 
            if strcmp(animal, 'PV106') && ismember(record, {'Hunter19', 'Hunter20', 'Hunter21', 'Hunter34'}), ch = 5;
            elseif strcmp(animal, 'PV157') && ismember(record, {'Hunter40', 'Hunter41'}), ch = 17; end
            
            bgRatioFile = [currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
            % bgRatioFile = [SA.currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
            S = load(bgRatioFile, 'f', 't_ms', 'allFreqProfiles');
            f = S.f; t_ms = S.t_ms; allFreqProfiles = S.allFreqProfiles;
            currentAnimal = animal; currentRecord = record;
        end
        
        stimOnset_ms = dataTable.Triger_ms(tr) + 5000;
        relTime = (t_ms - stimOnset_ms) / 1000;
        
        for fb = 1:nBands
            freqIdx = f >= bandRanges(fb, 1) & f <= bandRanges(fb, 2);
            if ~any(freqIdx), continue; end
            bandActivity = mean(allFreqProfiles(freqIdx, :), 1, 'omitnan');
            bandActivity_z = zscore(bandActivity);

            for tb = 1:nBins
                idx = relTime >= tEdges(tb) & relTime < tEdges(tb+1);
                if any(idx), binnedData(fb, tb, tr) = mean(bandActivity(idx), 'omitnan'); end
            end
        end
    end
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