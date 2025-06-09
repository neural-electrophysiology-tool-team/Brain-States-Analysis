function overlayPointsOnVideo(videoPath, matPath, xColumnName, yColumnName, outputVideoPath)
    % OVERLAYPOINTSONVIDEO Overlays coordinate points from .mat file on video frames
    %
    % Inputs:
    %   videoPath - string, path to input video file
    %   matPath - string, path to .mat file containing table with coordinates
    %   xColumnName - string, column name for X coordinates
    %   yColumnName - string, column name for Y coordinates  
    %   outputVideoPath - string, path for output video file (optional)
    %
    % Example:
    %   overlayPointsOnVideo('input.mp4', 'coords.mat', 'X', 'Y', 'output.mp4')
    
    % Set default output path if not provided
    if nargin < 5
        [path, name, ~] = fileparts(videoPath);
        outputVideoPath = fullfile(path, [name '_with_points.mp4']);
    end
    
    try
        % Load the .mat file
        fprintf('Loading .mat file: %s\n', matPath);
        matData = load(matPath);
        
        % Find the table in the loaded data
        tableVarName = [];
        fieldNames = fieldnames(matData);
        
        % Look for a table variable
        for i = 1:length(fieldNames)
            if isa(matData.(fieldNames{i}), 'table')
                tableVarName = fieldNames{i};
                break;
            end
        end
        
        if isempty(tableVarName)
            error('No table found in .mat file. Please ensure the .mat file contains a table variable.');
        end
        
        coordData = matData.(tableVarName);
        fprintf('Found table variable: %s\n', tableVarName);
        
        % Extract X and Y coordinates
        if ~ismember(xColumnName, coordData.Properties.VariableNames)
            error('Column "%s" not found in table. Available columns: %s', ...
                  xColumnName, strjoin(coordData.Properties.VariableNames, ', '));
        end
        if ~ismember(yColumnName, coordData.Properties.VariableNames)
            error('Column "%s" not found in table. Available columns: %s', ...
                  yColumnName, strjoin(coordData.Properties.VariableNames, ', '));
        end
        
        xCoords = coordData.(xColumnName);
        yCoords = coordData.(yColumnName);
        
        % Read the input video
        fprintf('Reading video file: %s\n', videoPath);
        videoReader = VideoReader(videoPath);
        
        % Check if number of frames matches CSV rows
        numFrames = videoReader.NumFrames;
        numCoordRows = height(coordData);
        
        if numFrames ~= numCoordRows
            if numCoordRows > numFrames
                fprintf('Table has more rows (%d) than video frames (%d). Using only first %d rows.\n', ...
                        numCoordRows, numFrames, numFrames);
                % Trim table data to match video frames
                coordData = coordData(1:numFrames, :);
                xCoords = xCoords(1:numFrames);
                yCoords = yCoords(1:numFrames);
            else
                error('Video has more frames (%d) than table rows (%d). Please ensure table has coordinates for all frames.', ...
                      numFrames, numCoordRows);
            end
        end
        
        % Create video writer
        fprintf('Creating output video: %s\n', outputVideoPath);
        
        % Choose appropriate profile based on file extension
        [~, ~, ext] = fileparts(outputVideoPath);
        
        switch lower(ext)
            case '.mp4'
                try
                    videoWriter = VideoWriter(outputVideoPath, 'MPEG-4');
                catch
                    % Fallback to Motion JPEG AVI if MPEG-4 not supported
                    fprintf('MPEG-4 not supported, using Motion JPEG AVI format\n');
                    outputVideoPath = strrep(outputVideoPath, '.mp4', '.avi');
                    videoWriter = VideoWriter(outputVideoPath, 'Motion JPEG AVI');
                end
            case '.avi'
                videoWriter = VideoWriter(outputVideoPath, 'Motion JPEG AVI');
            otherwise
                % Default to AVI with Motion JPEG
                fprintf('Using default AVI format with Motion JPEG compression\n');
                outputVideoPath = [outputVideoPath '.avi'];
                videoWriter = VideoWriter(outputVideoPath, 'Motion JPEG AVI');
        end
        
        videoWriter.FrameRate = videoReader.FrameRate;
        open(videoWriter);
        
        % Process each frame
        fprintf('Processing frames...\n');
        for frameIdx = 1:numFrames
            % Read frame
            frame = read(videoReader, frameIdx);
            
            % Get coordinates for this frame
            x = xCoords(frameIdx);
            y = yCoords(frameIdx);
            
            % Skip if coordinates are NaN or invalid
            if ~isnan(x) && ~isnan(y) && x > 0 && y > 0 && ...
               x <= size(frame, 2) && y <= size(frame, 1)
                
                % Draw point on frame
                frame = drawPoint(frame, x, y);
            end
            
            % Write frame to output video
            writeVideo(videoWriter, frame);
            
            % Progress indicator
            if mod(frameIdx, 50) == 0 || frameIdx == numFrames
                fprintf('Processed %d/%d frames (%.1f%%)\n', ...
                        frameIdx, numFrames, 100*frameIdx/numFrames);
            end
        end
        
        % Close video writer
        close(videoWriter);
        
        fprintf('Successfully created video with overlaid points: %s\n', outputVideoPath);
        
    catch ME
        % Clean up in case of error
        if exist('videoWriter', 'var') && isvalid(videoWriter)
            close(videoWriter);
        end
        rethrow(ME);
    end
end

function frame = drawPoint(frame, x, y, varargin)
    % DRAWPOINT Draws a point on the frame
    %
    % Optional parameters can be added for customization
    
    % Default point properties
    pointSize = 8;        % Radius of the point
    pointColor = [255, 0, 0]; % Red color (RGB)
    lineWidth = 2;
    
    % Parse optional arguments
    if nargin > 3
        for i = 1:2:length(varargin)
            switch lower(varargin{i})
                case 'size'
                    pointSize = varargin{i+1};
                case 'color'
                    pointColor = varargin{i+1};
                case 'linewidth'
                    lineWidth = varargin{i+1};
            end
        end
    end
    
    % Convert to integer coordinates
    x = round(x);
    y = round(y);
    
    % Get frame dimensions
    [height, width, channels] = size(frame);
    
    % Create circular mask for the point
    [X, Y] = meshgrid(1:width, 1:height);
    dist = sqrt((X - x).^2 + (Y - y).^2);
    
    % Create filled circle
    circleMask = dist <= pointSize;
    
    % Create ring for border if lineWidth > 1
    if lineWidth > 1
        borderMask = (dist <= pointSize) & (dist > pointSize - lineWidth);
        circleMask = borderMask;
    end
    
    % Apply color to the point
    for ch = 1:channels
        channel = frame(:, :, ch);
        channel(circleMask) = pointColor(ch);
        frame(:, :, ch) = channel;
    end
    
    % Add a small cross at the center for better visibility
    crossSize = max(2, round(pointSize/3));
    
    % Horizontal line of cross
    if y-crossSize >= 1 && y+crossSize <= height && x >= 1 && x <= width
        for ch = 1:channels
            frame(y-crossSize:y+crossSize, x, ch) = 255; % White cross
        end
    end
    
    % Vertical line of cross  
    if x-crossSize >= 1 && x+crossSize <= width && y >= 1 && y <= height
        for ch = 1:channels
            frame(y, x-crossSize:x+crossSize, ch) = 255; % White cross
        end
    end
end

% Example usage function
function runExample()
    % This function demonstrates how to use the main function
    % Create sample data for testing
    
    fprintf('Example usage:\n');
    fprintf('overlayPointsOnVideo(''input.mp4'', ''coordinates.mat'', ''X_pixel'', ''Y_pixel'', ''output.mp4'');\n\n');
    
    fprintf('.mat file should contain a table with coordinate data, for example:\n');
    fprintf('coordTable = table([1;2;3;...], [100;105;110;...], [150;148;152;...], ...\n');
    fprintf('                  ''VariableNames'', {''Frame'', ''X_pixel'', ''Y_pixel''});\n');
    fprintf('save(''coordinates.mat'', ''coordTable'');\n\n');
    
    fprintf('The function will:\n');
    fprintf('1. Read the video and CSV files\n');
    fprintf('2. Validate that dimensions match\n');
    fprintf('3. Overlay red circles with white crosses at specified coordinates\n');
    fprintf('4. Save the result as a new video file\n');
end