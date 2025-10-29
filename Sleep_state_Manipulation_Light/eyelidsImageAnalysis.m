%% image analysis for eyelids
% upload the pics for each animal and expusre.
% for each image in a stock, calculate the gausain of noise and reduce it
% from the image. let it be calculate from a specific area. 
% save the light intensity for each exposure and filter.
% calculate the differences. 
AnimalFolder = '/media/sil1/Data/Nitzan/Experiments/Eyelids_obsorption/PV88';
results = analyzeTiffFolders(AnimalFolder);
function resultsTable = analyzeTiffFolders(mainFolder)
% analyzeTiffFolders - Process TIFF images with noise correction and ROI analysis
%
% Syntax: resultsTable = analyzeTiffFolders(mainFolder)
%
% Input:
%   mainFolder - Path to main folder containing subfolders with structure:
%                ANIMAL_BF_TYPE_EXPOSURE
%
% Output:
%   resultsTable - Table with columns: Animal, Type, Exposure, Filter, 
%                  noiseGaus, ROI_Intensity

    % Initialize results storage
    results = struct('Animal', {}, 'BF', {}, 'Type', {}, 'Exposure', {}, ...
                     'Filter', {}, 'NoiseMean', {}, 'NoiseStd', {}, ...
                     'ROI_Intensity', {});
    resultIdx = 1;
    
    % Storage for ROI masks per animal-type combination
    roiStorage = struct();
    
    % Get all subfolders in main folder
    subfolders = dir(mainFolder);
    subfolders = subfolders([subfolders.isdir] & ~ismember({subfolders.name}, {'.', '..'}));
    
    fprintf('Found %d subfolders to process\n', length(subfolders));
    
    % Process each subfolder
    for i = 1:length(subfolders)
        folderName = subfolders(i).name;
        folderPath = fullfile(mainFolder, folderName);
        
        % Parse folder name: ANIMAL_BF_TYPE_EXPOSURE
        parts = strsplit(folderName, '_');
        if length(parts) < 4
            warning('Folder "%s" does not match expected format. Skipping.', folderName);
            continue;
        end
        
        animal = parts{1};
        bf = parts{2};  % Either 'BF' or 'noBF'
        type = parts{3};
        exposure = parts{4};
        
        % Find .tif files in folder
        tifFiles = dir(fullfile(folderPath, '*.tif'));
        if isempty(tifFiles)
            tifFiles = dir(fullfile(folderPath, '*.tiff'));
        end
        
        if isempty(tifFiles)
            warning('No TIFF files found in folder "%s". Skipping.', folderName);
            continue;
        end
        
        % Process each TIFF file (usually should be just one per folder)
        for j = 1:length(tifFiles)
            tifPath = fullfile(folderPath, tifFiles(j).name);
            fprintf('Processing: %s\n', tifPath);
            
            % Read TIFF stack information
            info = imfinfo(tifPath);
            numImages = length(info);
            
            fprintf('  Found %d images in TIFF stack\n', numImages);
            
            % Read first image to get dimensions and select ROIs
            firstImg = imread(tifPath, 1);
            
            % Convert to double for processing
            firstImg = double(firstImg);
            
            % Create unique key for this animal-BF-type combination
            roiKey = sprintf('%s_%s_%s', animal, bf, type);
            
            % Check if we already have ROIs and noise masks for this combination
            if isfield(roiStorage, roiKey)
                % Reuse existing masks
                noiseMask = roiStorage.(roiKey).noiseMask;
                roiMask = roiStorage.(roiKey).roiMask;
                fprintf('  Reusing ROIs for %s\n', roiKey);
            else
                % Need to select new ROIs
                fprintf('  Selecting ROIs for %s\n', roiKey);
                
                % Display first image and select noise region
                figure('Name', sprintf('Select NOISE region - %s', roiKey));
                imshow(firstImg, []);
                title('Draw rectangle for NOISE region, then double-click inside');
                noiseROI = drawrectangle();
                wait(noiseROI);
                noiseMask = createMask(noiseROI);
                delete(noiseROI);
                
                % Display first image and select ROI
                figure('Name', sprintf('Select ROI region - %s', roiKey));
                imshow(firstImg, []);
                title('Draw rectangle for ROI region, then double-click inside');
                roi = drawrectangle();
                wait(roi);
                roiMask = createMask(roi);
                delete(roi);
                
                close all;
                
                % Store masks for future use
                roiStorage.(roiKey).noiseMask = noiseMask;
                roiStorage.(roiKey).roiMask = roiMask;
            end
            
            % Process each image in the stack
            for k = 1:numImages
                % Read image
                img = imread(tifPath, k);
                img = double(img);
                
                % Extract noise region pixels
                noisePixels = img(noiseMask);
                
                % Calculate Gaussian parameters of noise
                noiseMean = mean(noisePixels);
                noiseStd = std(noisePixels);
                
                % Extract ROI pixels
                roiPixels = img(roiMask);
                
                % Subtract noise mean from ROI (noise correction)
                roiCorrected = roiPixels - noiseMean;
                
                % Calculate mean intensity of corrected ROI
                roiIntensity = mean(roiCorrected);
                
                % Store results
                results(resultIdx).Animal = animal;
                results(resultIdx).BF = bf;
                results(resultIdx).Type = type;
                results(resultIdx).Exposure = exposure;
                results(resultIdx).Filter = k;
                results(resultIdx).NoiseMean = noiseMean;
                results(resultIdx).NoiseStd = noiseStd;
                results(resultIdx).ROI_Intensity = roiIntensity;
                
                resultIdx = resultIdx + 1;
            end
        end
    end
    
    % Convert struct array to table
    if ~isempty(results)
        % % Create noiseGaus column combining mean and std
        % noiseGaus = cell(length(results), 1);
        % for i = 1:length(results)
        %     noiseGaus{i} = sprintf('�=%.2f, s=%.2f', ...
        %                            results(i).NoiseMean, results(i).NoiseStd);
        % end
        % 
        resultsTable = struct2table(results);
        % resultsTable.noiseGaus = noiseGaus;
        
        % Reorder columns to match requested format
        resultsTable = resultsTable(:, {'Animal', 'BF', 'Type', 'Exposure', ...
                                         'Filter', 'ROI_Intensity','NoiseMean','NoiseStd' });
        
        fprintf('\nProcessing complete! %d measurements recorded.\n', height(resultsTable));
        
        % Display first few rows
        disp('First few rows of results:');
        disp(head(resultsTable));
        
        % Save results to CSV
        outputFile = fullfile(mainFolder, 'results.csv');
        writetable(resultsTable, outputFile);
        fprintf('Results saved to: %s\n', outputFile);
        
        % Generate plots
        % plotResults(resultsTable, mainFolder);
    else
        warning('No data was processed.');
        resultsTable = table();
    end
end


%% Get the results together and plot:
Animals = {'PV106','PV88','PV143'};
generalFolder = '/media/sil1/Data/Nitzan/Experiments/Eyelids_obsorption';
results = [];
for i = 1:length(Animals)
    curFolder = [generalFolder '/' Animals{i} '/results.csv'];
    curResults = readtable(curFolder);
    results = [results;curResults];
end

results.Exposure = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), results.Exposure);
results.Type(strcmp(results.Type,'noEyelids'))={'noEyelid'};
results.ExpMulti = 5 ./ results.Exposure;
results.ExpRelaInten = results.ROI_Intensity.*results.ExpMulti;

noBFidx = contains(results.BF,'no');
results.Filter(noBFidx) = results.Filter(noBFidx)+1;


%save all results
generalFolder = '/media/sil1/Data/Nitzan/Experiments/Eyelids_obsorption';
outputFile = fullfile(generalFolder, 'results.csv');

writetable(results, outputFile);


%% after relt
filterLabels = ["BF","DAPI","GFP","mCherry","iRFP","YFP"];

figure;
idx = ~contains(results.Type,'no');
         
subData = results(idx, :);
plot(subData.Filter,subData.ExpRelaInten,'.'); hold on; 
idx = ~contains(results.BF, 'no') & ...
    contains(results.Type,'no');
         
% subData2 = results(idx, :);
% plot(subData2.Filter,subData2.ExpRelaInten,'.',color='red');    

ylabel('reltive intesity (to exposure)'),xticks([1:6]);
xticklabels(filterLabels)
legend("Eyelids","noEyelids"), xlim([0.7,6.3])
%% first, check if the exposure is now linear

% types = unique(results.Type);
% Animals = unique(results.Animal);
generalFolder = '/media/sil1/Data/Nitzan/Experiments/Eyelids_obsorption';
plotResults(results, generalFolder);
function plotResults(resultsTable, mainFolder)
    % plotResults - Create subplots for each animal-BF-type combination
    %
    % Each subplot shows different exposures with Filter on x-axis and 
    % ROI_Intensity on y-axis
    
    % Get unique animal-BF-type combinations
    uniqueCombos = unique(resultsTable(:, {'Animal', 'BF', 'Type'}), 'rows');
    numCombos = height(uniqueCombos);
    
    % Determine subplot layout
    numCols = ceil(sqrt(numCombos));
    numRows = ceil(numCombos / numCols);
    
    % Create figure
    figure('Name', 'ROI Intensity Analysis', 'Position', [100, 100, 400*numCols, 300*numRows]);
    
    % Plot each combination
    for i = 1:numCombos
        subplot(numRows, numCols, i);
        
        animal = uniqueCombos.Animal{i};
        bf = uniqueCombos.BF{i};
        type = uniqueCombos.Type{i};
        
        % Filter data for this combination
        idx = strcmp(resultsTable.Animal, animal) & ...
              strcmp(resultsTable.BF, bf) & ...
              strcmp(resultsTable.Type, type);
        subData = resultsTable(idx, :);
        
        % Get unique exposures
        uniqueExposures = unique(subData.Exposure);
        
        % Plot each exposure as a separate line
        hold on;
        colors = lines(length(uniqueExposures));
        
        for j = 1:length(uniqueExposures)
            exposure = uniqueExposures(j);
            expIdx = subData.Exposure == exposure;
            expData = subData(expIdx, :);
            
            % Sort by filter number
            [~, sortIdx] = sort(expData.Filter);
            expData = expData(sortIdx, :);
            
            % Plot
            plot(expData.Filter, expData.ROI_Intensity, '-o', ...
                 'LineWidth', 2, 'MarkerSize', 8, ...
                 'Color', colors(j,:), 'DisplayName', num2str(exposure));
        end
        
        % Formatting
        xlabel('Filter');
        ylabel('ROI Intensity');
        title(sprintf('%s - %s - %s', animal, bf, type), 'Interpreter', 'none');
        legend('Location', 'best');
        grid on;
        hold off;
    end
    
    % Save figure
    saveas(gcf, fullfile(mainFolder, 'ROI_Intensity_Plots.png'));
    saveas(gcf, fullfile(mainFolder, 'ROI_Intensity_Plots.fig'));
    fprintf('Plots saved to: %s\n', mainFolder);
end
%% plot exposure to intensity for each wavelength:
filterLabels = ["BF","DAPI-460","GFP-525","mCherry-630","iRFP-700","YFP-540"];
uniqueExposures = unique(results.Exposure);
colors = parula(length(uniqueExposures));  % or parula(8), hsv(8), etc.
figure;
for i = 1:length(filterLabels)
    % get the sub set of data:
    idx= results.Filter ==i &...
        ~contains(results.Type,'noEyelid');
    subData=results(idx,:);

    subplot(3,2,i)
    plot(subData.Exposure,subData.ROI_Intensity,'k-*'); hold on;
    %plot regresion
    p = polyfit(subData.Exposure,subData.ROI_Intensity,1);
    yfit = polyval(p,subData.Exposure);
    plot(subData.Exposure,yfit,'r-',LineWidth=1.5)
    xlim([-5 1005])
    title(filterLabels(i))
    ylabel('Intensity');xlabel('Exposure (ms)')

end


    saveas(gcf, fullfile(generalFolder, 'Intensity-Exposure.png'));
    % saveas(gcf, fullfile(mainFolder, 'ROI_Intensity_Plots.fig'));
    fprintf('Plots saved to: %s\n', generalFolder);

%% plot reltive intesity to exposure
filterLabels = ["BF","DAPI-460","GFP-525","mCherry-630","iRFP-700","YFP-540"];


figure;
idx = ~contains(results.BF, 'no') & ...
    ~contains(results.Type,'no');
         
subData = results(idx, :);
plot(subData.Filter,subData.ExpRelaInten,'.'); hold on; 
idx = ~contains(results.BF, 'no') & ...
    contains(results.Type,'no');
% 
% subData2 = results(idx, :);
% plot(subData2.Filter,subData2.ExpRelaInten,'.',color='red');    

ylabel('reltive intesity (to exposure)'),xticks([1:6]);
xticklabels(filterLabels)
legend("Eyelids","noEyelids"), xlim([0.7,6.3])



%% calculate the relative light trasferation. 

% function plotRelativeTransmition(resultsTable, resultsFolder)
    results.transmision = NaN(height(results),1);
    % Get unique animal-BF-type combinations
    uniqueCombos = unique(resultsTable(:, {'Animal', 'BF','Filter'}), 'rows');
    numCombos = height(uniqueCombos);

    % calculate for each combination
    for i = 1:numCombos
        
        animal = uniqueCombos.Animal{i};
        bf = uniqueCombos.BF{i};
        filter = uniqueCombos.Filter(i);

        % Filter data for this combination 
        % indexs of all the eyelids with the same parameters and filter
        idxEye = strcmp(results.Animal, animal) & ...
              strcmp(results.BF, bf) & ...
              results.Filter ==filter & ...
              ~contains(results.Type,'no');
        % finding the one no eyelid relevant
        idxnoEye = strcmp(results.Animal, animal) & ...
              strcmp(results.BF, bf) & ...
              results.Filter ==filter & ...
              contains(results.Type,'no');
        %calculate relative transmition
        curbasline = results.ExpRelaInten(idxnoEye);
        results.transmision(idxEye) = results.ExpRelaInten(idxEye) / curbasline;
  
        
    end

 generalFolder = '/media/sil1/Data/Nitzan/Experiments/Eyelids_obsorption';
outputFile = fullfile(generalFolder, 'results.csv');
writetable(results, outputFile);
%% plot reltive tranmition
filterLabels = ["BF","DAPI-460","GFP-525","mCherry-630","iRFP-700","YFP-540"];
Animals = unique(results.Animal);
markers = {'s','o','*'};
uniqueExposures = unique(results.Exposure);
% colors = lines(length(uniqueExposures));
colors = parula(length(uniqueExposures));  % or parula(8), hsv(8), etc.
% colororder(colors);
figure;
for i =1:length(Animals)
    idx = ~contains(results.Type,'no')&...
        strcmp(results.Animal,Animals(i));
    subData = results(idx,:);        
    
        
        % Plot each exposure as a separate color
        hold on;
        curUnique=unique(subData.Exposure);
        for j = 1:length(curUnique)
            exposure = uniqueExposures(j);
            expIdx = subData.Exposure == exposure;
            expData = subData(expIdx, :);
            expIndices = ismember(uniqueExposures, exposure);
            expColor = colors(expIndices, :);
    
            plot(expData.Filter,expData.transmision,'.',Color=expColor,MarkerSize=15); hold on;
        end


end         

ylabel('% Transmition'),xticks([1:6]);
xticklabels(filterLabels)
legend(num2str(uniqueExposures),Location="northeast"), xlim([0.7,6.3])
ylim([0,0.05])
% save image
generalFolder = '/media/sil1/Data/Nitzan/Experiments/Eyelids_obsorption';
outputFile = fullfile(generalFolder, 'Tranmission.png');
saveas(gcf, outputFile);