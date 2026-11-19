function plotPCAFigure(pcaData, freqHzAll,  clustWake, clustSleep, ...
                       nClustersWake, nClustersSleep, clustersNames, ...
                       clusterColorsWake, clusterColorsSleep)

    for r = 1:3
        coeff_r       = pcaData{r,1};
        scoresWake_r  = pcaData{r,2};
        % scoresHun_r   = pcaData{r,3};
        scoresSleep_r = pcaData{r,3};
        explained_r   = pcaData{r,4};
        rowLabel      = pcaData{r,5};
        dx1_r         = pcaData{r,6};
        dx2_r         = pcaData{r,7};
        xlim_r        = pcaData{r,8};
        ylim_r        = pcaData{r,9};

        %-- Col 1: PC Loadings
        subplot(3, 5, (r-1)*5 + 1);
        hold on;
        for pc = 1:3
            plot(freqHzAll(1:size(coeff_r,1)), coeff_r(:,pc), ...
                'DisplayName', sprintf('PC%d (%.1f%%)', pc, explained_r(pc)));
        end
        xlabel('Frequency (Hz)');
        ylabel('Loading');
        title([rowLabel ' - Loadings']);
        legend;

        % %-- Col 2: hist2 Hunter
        % subplot(3, 6, (r-1)*6 + 2);
        % [dHun_r, hHun_r] = hist2(scoresHun_r(:,1), scoresHun_r(:,2), ...
        %     'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
        %     'limits1',xlim_r,'limits2',ylim_r,'plotProb',1);
        % xlabel('PC1'); ylabel('PC2');
        % title([rowLabel ' - Hunter']);
        % hHun_r.XLim = xlim_r; hHun_r.YLim = ylim_r;

        %-- Col 3: hist2 Wake
        subplot(3, 5, (r-1)*5 + 2);
        [dWake_r, hWake_r] = hist2(scoresWake_r(:,1), scoresWake_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r,'plotProb',1);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Wake']);
        hWake_r.XLim = xlim_r; hWake_r.YLim = ylim_r;

        %-- Col 4: hist2 Sleep
        subplot(3, 5, (r-1)*5 + 3);
        [dSleep_r, hSleep_r] = hist2(scoresSleep_r(:,1), scoresSleep_r(:,2), ...
            'dX1',dx1_r,'dX2',dx2_r,'h',gca, ...
            'limits1',xlim_r,'limits2',ylim_r,'plotProb',1);
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Sleep']);
        hSleep_r.XLim = xlim_r; hSleep_r.YLim = ylim_r;

        %-- Col 5: Overlaid RGB hist2
        subplot(3, 5, (r-1)*5 + 4);
        ax = gca;

        edges1_r = (xlim_r(1) : dx1_r : (xlim_r(2)+dx1_r)) + dx1_r/2;
        edges2_r = (ylim_r(1) : dx2_r : (ylim_r(2)+dx2_r)) + dx2_r/2;

        % nHun   = mat2gray(log10(1 + dHun_r));
        nWake  = mat2gray(log10(1 + dWake_r));
        nSleep = mat2gray(log10(1 + dSleep_r));
        nHun = zeros(size(nSleep));

        rgbCombined = cat(3, max(nWake,0), max(nSleep,0), max(nHun,0));

        imagesc(ax, edges1_r, edges2_r, rgbCombined);
        set(ax, 'YDir','normal', 'Color','white');
        ax.XLim = xlim_r; ax.YLim = ylim_r;
        xlabel('PC1'); ylabel('PC2');
        title([rowLabel ' - Overlap: Wake=red, Hunter=blue, Sleep=green']);

        %-- Col 6: Scatter colored by cluster
        subplot(3, 5, (r-1)*5 + 5);
        hold on;

        for c = 1:nClustersWake
            % scatter(scoresHun_r(clustHun==c, 1), scoresHun_r(clustHun==c, 2), ...
            %     10, clusterColorsWakeHun(c,:), 'diamond', ...
            %     'MarkerEdgeAlpha',0.2,'LineWidth',1, ...
            %     'DisplayName', ['Hun C' num2str(c)]);
            scatter(scoresWake_r(clustWake==c, 1), scoresWake_r(clustWake==c, 2), ...
                10, clusterColorsWake(c,:), '*', ...
                'MarkerEdgeAlpha',0.2, ...
                'DisplayName', ['Wake C' num2str(c)]);
        end

        validSleep = ~cellfun(@isempty, clustSleep);
        for c = 1:nClustersSleep
            maskC_sleep = strcmp(string(clustSleep(validSleep)), clustersNames{c});
            scatter(scoresSleep_r(maskC_sleep,1), scoresSleep_r(maskC_sleep,2), ...
                5, clusterColorsSleep(c,:), 'o', ...
                'MarkerEdgeAlpha',0.25, ...
                'DisplayName', clustersNames{c});
        end

        xlabel('PC1'); ylabel('PC2');
        xlim(xlim_r); ylim(ylim_r);
        legend('Location','best','FontSize',6);
        title([rowLabel ' - Clusters']);
        hold off;
    end
end