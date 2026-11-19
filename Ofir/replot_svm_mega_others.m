%% ============================================================
%  RE-PLOT SAVED SVM RESULTS (category_pairs)   -- run svm_mega.m first
%   - loads svm_allbands_rawnorm_category_pairs.mat
%   - accuracy figure per group of comparisons: star per animal + mean, boxplot of the permutations,
%     NO p-values; every panel title shows the window length of that comparison
%   - summary figure / heatmap are OFF by default (showSummary / showHeatmap)
%   - all figures -> one PowerPoint (one slide per figure, slide title = figure name)
% =============================================================

%% ---------------- parameters ----------------
matFile  = 'svm_allbands_rawnorm_category_pairs.mat';
excludeComparisons = {'sound | Onset vs Offset'};
animalColor = [0.20 0.45 0.80];
showSummary = true;                                  % summary figure(s) of the comparisons
summarySplit = 'group';                               % 'group' = one summary figure per group of comparisons (= one per shaded
                                                      %           sub-section of the old combined figure)
                                                      % 'all'   = one big figure with all comparisons (gray bands between groups)
showHeatmap = true;                                  % animals x comparisons heatmap
savePPT  = true;
pptFile  = 'svm_figures.pptx';

%% ---------------- load ----------------
S = load(matFile);
R = S.batchResults;
R = R([R.ModelIdx] == 1);
R = R(~ismember({R.TrialType}, excludeComparisons));
if isempty(R), error('No results left to plot (check matFile / excludeComparisons).'); end
if isfield(S, 'nSplits'), nSplits = S.nSplits; else, nSplits = 1; end

% settings of the SVM run (saved by svm_mega.m)
if isfield(S, 'winLen')
    fprintf('Window lengths: default %g s | circle vs other %g s | longWinSets (%s) %g s\n', S.winLen, S.winLenShort, ...
            strjoin(S.longWinSets, ','), S.winLenLong);
    fprintf('Pre-offset window ends %g s from the offset\n', S.preOffsetEnd);
end
if isfield(S, 'otherGroups')
    for g = 1:numel(S.otherGroups)
        fprintf('  %-48s types = {%s}\n', S.otherGroups(g).name, strjoin("'" + string(S.otherGroups(g).types) + "'", ','));
    end
end

animals = unique({R.Animal}, 'stable');
nA = numel(animals);

% comparisons, ordered by group (group order as saved by svm_mega.m)
cmpRaw = unique({R.TrialType}, 'stable');
grpRaw = cellfun(@(s) char(extractBefore(string(s), ' | ')), cmpRaw, 'UniformOutput', false);
if isfield(S, 'groupOrder') && ~isempty(S.groupOrder), canon = S.groupOrder; else, canon = {}; end
grpNames = unique(grpRaw, 'stable');
grpRank  = zeros(size(grpNames));
for g = 1:numel(grpNames)
    k = find(strcmp(canon, grpNames{g}), 1);
    if isempty(k), k = numel(canon) + g; end
    grpRank(g) = k;
end
[~, o] = sort(grpRank);  grpNames = grpNames(o);
cmpAll = {};  cmpGroup = {};
for g = 1:numel(grpNames)
    idx = find(strcmp(grpRaw, grpNames{g}));
    cmpAll   = [cmpAll,   cmpRaw(idx)];                          %#ok<AGROW>
    cmpGroup = [cmpGroup, repmat(grpNames(g), 1, numel(idx))];   %#ok<AGROW>
end
nC = numel(cmpAll);

% accuracy [animals x comparisons], null [permutations x animals x comparisons], window length per comparison
acc    = nan(nA, nC);
cmpLen = nan(1, nC);
nPerm  = max([cellfun(@numel, {R.nullAcc}), 0]);
nul    = nan(max(nPerm, 1), nA, nC);
for i = 1:numel(R)
    a = find(strcmp(animals, R(i).Animal));
    c = find(strcmp(cmpAll,  R(i).TrialType));
    acc(a, c) = R(i).acc;
    if isfield(R, 'winLen') && ~isempty(R(i).winLen), cmpLen(c) = R(i).winLen; end
    if ~isempty(R(i).nullAcc), nul(1:numel(R(i).nullAcc), a, c) = R(i).nullAcc(:); end
end
hasNull = nPerm > 0 && any(~isnan(nul(:)));

%% ---------------- figures ----------------
figList = gobjects(0);

for g = 1:numel(grpNames)
    idx = find(strcmp(cmpGroup, grpNames{g}));
    figList(end+1) = plotGroupAccuracy(grpNames{g}, cmpAll(idx), cmpLen(idx), acc(:, idx), nul(:, :, idx), ...
                                       animals, hasNull, animalColor, nSplits); %#ok<SAGROW>
end
if showSummary
    if strcmp(summarySplit, 'group')
        for g = 1:numel(grpNames)
            idx = find(strcmp(cmpGroup, grpNames{g}));
            figList(end+1) = plotSummaryAll(cmpAll(idx), cmpGroup(idx), acc(:, idx), nul(:, :, idx), ...
                                            animals, hasNull, animalColor, nSplits, grpNames{g}); %#ok<SAGROW>
        end
    else
        figList(end+1) = plotSummaryAll(cmpAll, cmpGroup, acc, nul, animals, hasNull, animalColor, nSplits, 'all comparisons');
    end
end
if showHeatmap
    figList(end+1) = plotHeatmap(cmpAll, cmpGroup, acc, animals);
end

%% ---------------- PowerPoint ----------------
if savePPT
    saveFiguresToPPT(figList, pptFile);
end


%% ============================================================
%  HELPERS
% =============================================================

function fh = plotGroupAccuracy(gName, cmpNames, cmpLen, acc, nul, animals, hasNull, animalColor, nSplits)
% one panel per comparison: star per animal (+ Mean), gray boxplot of the permutations, no p-values
    nA = numel(animals);  showMean = nA > 1;  nX = nA + showMean;
    labels = animals;  if showMean, labels{end+1} = 'Mean'; end
    nT = numel(cmpNames);
    nCol = min(nT, 3);  nRow = ceil(nT / nCol);
    nPerm = size(nul, 1);

    fh = figure('Name', ['Accuracy: ' gName], 'Color', 'w', ...
                'Position', [100 100 300 + 450 * nCol, 80 + 420 * nRow]);
    tl = tiledlayout(fh, nRow, nCol, 'TileSpacing', 'compact');
    hLeg = [];  axLeg = [];
    for k = 1:nT
        ax = nexttile(tl);  hold(ax, 'on');
        ttl = char(extractAfter(string(cmpNames{k}), ' | '));
        if ~isnan(cmpLen(k)), ttl = sprintf('%s  [%g s windows]', ttl, cmpLen(k)); end
        title(ax, ttl, 'Interpreter', 'none', 'FontSize', 10);
        ylim(ax, [0 1]);  xlim(ax, [0.5 nX + 0.5]);  box(ax, 'on');
        grid(ax, 'off');  ax.YGrid = 'on';
        for xb = 1.5:1:(nX - 0.5)
            xline(ax, xb, '-', 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
        end
        xticks(ax, 1:nX);  xticklabels(ax, labels);  ax.TickLabelInterpreter = 'none';
        yline(ax, 0.5, 'k--', 'Chance (50%)', 'HandleVisibility', 'off');
        if mod(k - 1, nCol) == 0, ylabel(ax, sprintf('Accuracy (mean of %d splits)', nSplits)); end

        a = acc(:, k)';
        if all(isnan(a))
            text(ax, 0.5, 0.5, 'no data', 'Units', 'normalized', 'HorizontalAlignment', 'center');
            continue;
        end
        accRow = a;
        if showMean, accRow = [a, mean(a, 'omitnan')]; end

        h = [];
        if hasNull
            nl = nul(:, :, k);
            if showMean, nl = [nl, mean(nl, 2, 'omitnan')]; end
            xg = repelem(1:nX, nPerm);
            hB = boxchart(ax, xg(:), nl(:), 'BoxFaceColor', [.6 .6 .6], 'MarkerColor', [.4 .4 .4], ...
                          'BoxWidth', 0.3, 'DisplayName', sprintf('Null (%d perm.)', nPerm));
            h = hB;
        end
        hS = plot(ax, 1:nA, a, 'p', 'Color', animalColor, 'MarkerFaceColor', animalColor, ...
                  'MarkerSize', 10, 'DisplayName', 'Accuracy (all bands)');
        h = [h hS]; %#ok<AGROW>
        if showMean
            plot(ax, nX, accRow(end), 'p', 'Color', 'k', 'MarkerFaceColor', animalColor, ...
                 'MarkerSize', 14, 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        if isempty(hLeg), hLeg = h; axLeg = ax; end
    end
    title(tl, sprintf('SVM accuracy (all bands): %s', gName), 'Interpreter', 'none');
    if ~isempty(hLeg)
        lgd = legend(axLeg, hLeg);
        lgd.Layout.Tile = 'east';
    end
end


function fh = plotSummaryAll(cmpAll, cmpGroup, acc, nul, animals, hasNull, animalColor, nSplits, secName)
    nA = numel(animals);  nC = numel(cmpAll);
    fh = figure('Name', ['Summary: ' secName], 'Color', 'w', ...
                'Position', [50 50 min(2600, 500 + 90 * nC) 850]);
    ax = axes(fh);  hold(ax, 'on');

    gNames = unique(cmpGroup, 'stable');
    for g = 1:numel(gNames)
        idx = find(strcmp(cmpGroup, gNames{g}));
        if mod(g, 2) == 1 && numel(gNames) > 1
            patch(ax, [min(idx) - .5, max(idx) + .5, max(idx) + .5, min(idx) - .5], [0 0 1 1], [.93 .93 .93], ...
                  'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end

    legH = [];
    if hasNull
        nPerm = size(nul, 1);
        nm = reshape(mean(nul, 2, 'omitnan'), nPerm, nC);
        xg = repelem(1:nC, nPerm);
        hB = boxchart(ax, xg(:), nm(:), 'BoxFaceColor', [.6 .6 .6], 'MarkerColor', [.4 .4 .4], ...
                      'BoxWidth', 0.5, 'DisplayName', sprintf('Null of the mean (%d perm.)', nPerm));
        legH(end+1) = hB; %#ok<AGROW>
    end
    if nA > 1, offs = linspace(-0.3, 0.3, nA); else, offs = 0; end
    for a = 1:nA
        plot(ax, (1:nC) + offs(a), acc(a, :), 'o', 'Color', animalColor, 'MarkerFaceColor', animalColor, ...
             'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    legH(end+1) = plot(ax, nan, nan, 'o', 'Color', animalColor, 'MarkerFaceColor', animalColor, ...
                       'DisplayName', sprintf('Animals (n=%d)', nA)); %#ok<AGROW>
    if nA > 1
        legH(end+1) = plot(ax, 1:nC, mean(acc, 1, 'omitnan'), 's', 'Color', 'k', 'MarkerFaceColor', 'k', ...
                           'MarkerSize', 7, 'DisplayName', 'Mean'); %#ok<AGROW>
    end

    yline(ax, 0.5, 'k--', 'HandleVisibility', 'off');
    ylim(ax, [0 1]);  xlim(ax, [0.5 nC + 0.5]);  box(ax, 'on');  ax.YGrid = 'on';
    xticks(ax, 1:nC);  xticklabels(ax, cmpAll);  xtickangle(ax, 60);
    ax.TickLabelInterpreter = 'none';  ax.FontSize = 7;
    ylabel(ax, sprintf('Accuracy (mean of %d splits)', nSplits));
    title(ax, ['SVM accuracy (all bands) - ' secName], 'Interpreter', 'none');
    legend(ax, legH, 'Location', 'northeastoutside');
end


function fh = plotHeatmap(cmpAll, cmpGroup, acc, animals)
    nA = numel(animals);  nC = numel(cmpAll);
    M = acc;  rowLab = animals;
    if nA > 1, M = [acc; mean(acc, 1, 'omitnan')]; rowLab = [animals, {'Mean'}]; end
    fh = figure('Name', 'Summary heatmap: accuracy of all comparisons', 'Color', 'w', ...
                'Position', [50 50 min(2600, 400 + 55 * nC) 400 + 40 * size(M, 1)]);
    ax = axes(fh);
    im = imagesc(ax, M);  set(im, 'AlphaData', ~isnan(M));
    ax.Color = 'w';  colormap(ax, parula);  caxis(ax, [0.4 1]);
    cb = colorbar(ax);  cb.Label.String = 'Accuracy';
    hold(ax, 'on');
    for r = 1:size(M, 1)
        for c = 1:nC
            if isnan(M(r, c)), continue; end
            if M(r, c) > 0.7, tc = 'k'; else, tc = 'w'; end
            text(ax, c, r, sprintf('%.2f', M(r, c)), 'HorizontalAlignment', 'center', 'FontSize', 6, 'Color', tc);
        end
    end
    gNames = unique(cmpGroup, 'stable');
    for g = 1:numel(gNames) - 1
        xline(ax, max(find(strcmp(cmpGroup, gNames{g}))) + 0.5, 'k-', 'LineWidth', 1.5);
    end
    if nA > 1, yline(ax, nA + 0.5, 'k-', 'LineWidth', 1.5); end
    xticks(ax, 1:nC);  xticklabels(ax, cmpAll);  xtickangle(ax, 60);
    yticks(ax, 1:size(M, 1));  yticklabels(ax, rowLab);
    ax.TickLabelInterpreter = 'none';  ax.FontSize = 7;
    title(ax, 'SVM accuracy (all bands) - animals x comparisons', 'Interpreter', 'none');
end


function saveFiguresToPPT(figs, pptFile)
% one slide per figure; slide title = the figure's Name. Needs the MATLAB Report Generator.
    import mlreportgen.ppt.*
    figs = figs(isgraphics(figs, 'figure'));
    if isempty(figs), fprintf('No figures to save.\n'); return; end

    maxW = 9;  maxH = 5.6;
    ppt = Presentation(pptFile);
    open(ppt);
    tmpDir = tempname;  mkdir(tmpDir);
    for k = 1:numel(figs)
        f = figs(k);
        ttl = f.Name;
        if isempty(ttl), ttl = sprintf('Figure %d', k); end
        png = fullfile(tmpDir, sprintf('fig_%03d.png', k));
        drawnow;
        try
            exportgraphics(f, png, 'Resolution', 200);
        catch
            print(f, png, '-dpng', '-r150');
        end

        info = imfinfo(png);
        w = maxW;  h = w * info.Height / info.Width;
        if h > maxH, h = maxH;  w = h * info.Width / info.Height; end
        pic = Picture(png);
        pic.Width  = sprintf('%.2fin', w);
        pic.Height = sprintf('%.2fin', h);

        try
            slide = add(ppt, 'Title and Picture');
            replace(slide, 'Title', ttl);
            replace(slide, 'Picture', pic);
        catch
            slide = add(ppt, 'Title and Content');
            replace(slide, 'Title', ttl);
            replace(slide, 'Content', pic);
        end
    end
    close(ppt);
    rmdir(tmpDir, 's');
    fprintf('Saved %d figures to %s\n', numel(figs), fullfile(pwd, pptFile));
end