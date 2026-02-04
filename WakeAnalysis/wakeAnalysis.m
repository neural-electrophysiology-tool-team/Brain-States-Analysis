

classdef wakeAnalysis < sleepAnalysis

    methods
%% class constructor
        function obj=wakeAnalysis(xlsFile)
           
           obj=obj@sleepAnalysis(xlsFile);
        end
        
        %% Spectrogam - under construction
        % this method is made to show the spectrogrm for the waking time
        
        % what sould it do?
        %  - get the recorded signal
        %  - get the channel from the table
        %  - load the sleep time from the relevant analysis (i think it was D2BAC)
        %  - take the data in the specific timing and do the spectrogram. 
        %  - save the mat file
        %
        
        function obj = getSpectrogram(obj, varargin)
            
            % needs the out of the DBAC analysis
            
            %checking that there is a current rec:
            obj.checkFileRecording;
            
            % used in Mark's functions. create a parser for the inputs and
            % results.
            
            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.recTable.defaulLFPCh(obj.currentPRec),@isnumeric);
            addParameter(parseObj,'win',1000,@isnumeric);
            addParameter(parseObj,'tStart',0,@isnumeric);
            %addParameter(parseObj,'overwrite',0,@isnumeric);       
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            
            %evaluate all input parameters in workspace - don't know what
            %it does.
            
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            %make parameter structure
            parSpectrogram=parseObj.Results;

            
            % get the wake time:
            % load the output of DB autocor, that holds the param for sleep
            % times
            DBAutocorFile=[obj.currentAnalysisFolder filesep 'dbAutocorr_ch' num2str(ch) '.mat'];
            load(DBAutocorFile)
            
            % get the signal
            signal = obj.currentDataObj.getData(ch, tStart,1000000);
            %get the spectro
            [s] = spectrogram (signal, win);
            figure;
            plot (s)
            
        end

        %% getHead position:
        function locTable = getHeadPosition(obj)
            % read the parquet 
            videoPath = obj.recTable.VideoFiles(obj.currentPRec);
            [videosFolderPath,vidName,~] = fileparts(videoPath{1});
            locFilename = strcat(videosFolderPath,'/predictions/front_head_ephys_resnet_101__',vidName,'.parquet');

            % check file exists:
            if ~isfile(locFilename)
                error('File does not exist: %s', locFilename);
            end
            locTable = parquetread(locFilename);
            ArenaCSV = obj.getArenaCSVs;
            locTable.t_ms = ArenaCSV.oeCamTrigs(1:height(locTable));


            obj.files.locTable = [obj.currentAnalysisFolder filesep 'locTable.mat'];
            save(obj.files.locTable,"locTable",'-mat')

        end



        %% get pupil diameter
        
        function [rightData, leftData] = getPupilData(obj)
            blockPath = obj.recTable.blockPath{obj.currentPRec};
            eyePath = [blockPath filesep 'Eye_data'];
            rightDataPath = [eyePath filesep 'right_eye_data.csv'];
            leftDataPath = [eyePath filesep 'left_eye_data.csv'];
            rightData = readtable(rightDataPath);
            leftData = readtable(leftDataPath);

        end

        %% Dendrogram + spectral profiles:
        % re-write the pwelch and dendrogram:
        function [] = plotHunterDendroSpectral(obj,varargin)
          % parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            parseObj.FunctionName='wakeAnalysis\plotHunterDendroSpectral';
            
            addParameter(parseObj,'ch',obj.recTable.defaulLFPCh(obj.currentPRec),@isnumeric);
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'trials',1:20,@isnumeric); %trials to incluse in analysis
            addParameter(parseObj,'win',5000,@isnumeric); % before and after start trial, ms
            addParameter(parseObj,'inputParams',false,@isnumeric);
            addParameter(parseObj,'maxDendroClusters',2,@isnumeric);
            addParameter(parseObj,'fs',obj.currentDataObj.samplingFrequency(1),@isnumeric); 
            addParameter(parseObj,'welchOL',0.5,@isnumeric); 
            addParameter(parseObj,'welchWin',5000,@isnumeric); 
            addParameter(parseObj,'fMax',110,@isnumeric);
          

            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
            
            ArenaCSVs = obj.getArenaCSVs(1);
            startTrigSh = ArenaCSVs.startTrigSh;
            if length(startTrigSh)~=20
                trials = 1:length(startTrigSh);
            end
            
            % get the data:
           
            % Bug apperance: before and after
            postStart = obj.currentDataObj.getData(ch,startTrigSh(trials)+400,win);
            preStart = obj.currentDataObj.getData(ch,startTrigSh(trials)-win, win);
            prePost = [squeeze(preStart);squeeze(postStart)];
            
            % welch parameters:
            welchWin = 1*fs;
            samplesOL = welchOL*welchWin;
            
            %calculate the ppx:
            [pxx,f] = pwelch(prePost',welchWin,samplesOL,[],fs);
            
            % dendrogram:
            p=find(f<fMax);
            pp=find(sum(pxx(p,:))<0.4e6); %reject signals with very high amplitudes (probably noise)
            sPxx=pxx(p,pp);
            freqHz=f(p);
            normsPxx=bsxfun(@rdivide,sPxx,mean(sPxx,2));
            corrMat=corrcoef(normsPxx);

            if maxDendroClusters==2

                [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);

                S1=mean(normsPxx(:,clusters==1),2);
                S2=mean(normsPxx(:,clusters==2),2);
                if mean(S1(1:3))>mean(S2(1:3))
                    crossFreq=freqHz(find(S2-S1>=0,1,'first'));
                else
                    crossFreq=freqHz(find(S1-S2>=0,1,'first'));
                end
            else
                [DC,order,clusters]=DendrogramMatrix(corrMat,'linkMetric','euclidean','linkMethod','ward','maxClusters',maxDendroClusters);

                for i=1:maxDendroClusters
                    S(:,i)=mean(normsPxx(:,clusters==i),2);
                end
                crossFreq=[];
            end

        %not normelized PSD:
        sPxx1 = mean(sPxx(:,clusters==1),2);
        sPxx2 = mean(sPxx(:,clusters==2),2);

        % plot
        k = length(trials); % Threshold
        % Display the correlation matrix
        fig1 = figure;
        set(fig1, 'Position', [100, 100, 900, 800]); % [x, y, width, height]

        imagesc(DC);
        colorbar;
        % axis equal;
        hold on;

          % Loop through each segment number and display it
        for i = 1:length(order)
            if order(i) <= k
                text_color = 'black';
                background_color = [0.7, 0.7, 1]; % Light purple background
            else
                text_color = 'black';
                background_color = [1, 0.85, 0.7]; % Light orange background
            end

            % Display the segment number with background
            text(-0.5, i, num2str(order(i)), 'HorizontalAlignment', 'center', ...
                'BackgroundColor', background_color, 'Color',text_color , ...
                'FontWeight', 'bold');
            text(i, -0.5, num2str(order(i)), 'HorizontalAlignment', 'center', ...
                'BackgroundColor', background_color, 'Color',text_color , ...
                'FontWeight', 'bold');

        end
        % Adjust limits
        % xlim([-1, length(order)]);
        % ylim([0.5, length(order) + 0.5]);
        set(gca, 'XTick', [], 'YTick', []);
        hold off;
        set(fig1,'PaperPosition',[.25 3 8 6])
        saveas(gcf,strcat(obj.currentPlotFolder, '/corelogram.pdf'));

        %plot the freqs:
        figure;
        set(gcf, 'Position', [300, 300, 700, 500]); % [x, y, width, height]
        plot(freqHz,S1,'Color',[0.7, 0.7, 1], 'LineWidth',2);
        hold on
        plot(freqHz,S2,'Color',[1, 0.85, 0.7],'LineWidth',2);
        legend('Quiet','Active')

        %save Figure
        set(gcf,'PaperPosition',[.25 3 8 6])
        saveas(gcf,strcat(obj.currentPlotFolder, '/freqBands.pdf'));
        end

        %% getBeta2GamnaRatio
        function data=getBeta2GammaRatio(obj,varargin)
            obj.checkFileRecording;

            parseObj = inputParser;
            addParameter(parseObj,'ch',obj.recTable.defaulLFPCh(obj.currentPRec),@isnumeric);
            addParameter(parseObj,'avgOnCh',[],@isnumeric); %uses several averaged channels for d/b extraction
            addParameter(parseObj,'movLongWin',1000*60*30,@isnumeric); %max freq. to examine
            addParameter(parseObj,'movWin',2000,@isnumeric);
            addParameter(parseObj,'movOLWin',1000,@isnumeric);
            addParameter(parseObj,'segmentWelch',1000,@isnumeric);
            addParameter(parseObj,'dftPointsWelch',2^10,@isnumeric);
            addParameter(parseObj,'OLWelch',0.5);
            addParameter(parseObj,'tStart',0,@isnumeric);
            addParameter(parseObj,'win',0,@isnumeric); %if 0 uses the whole recording duration
            addParameter(parseObj,'gammaLow',60,@isnumeric);
            addParameter(parseObj,'gammaHigh',80,@isnumeric);
            addParameter(parseObj,'betaLow',10,@isnumeric);
            addParameter(parseObj,'betaHigh',30,@isnumeric);
            addParameter(parseObj,'applyNotch',0,@isnumeric);
            addParameter(parseObj,'saveSpectralProfiles',0,@isnumeric);
            addParameter(parseObj,'maxVoltage',1000,@isnumeric);
            addParameter(parseObj,'overwrite',0,@isnumeric);
            addParameter(parseObj,'inputParams',false,@isnumeric);
            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end

            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end

            %make parameter structure
            par=parseObj.Results;

            if isnan(ch)
                disp('Error: no reference channel for beta2gamma extraction');
                return;
            end
            %check if analysis was already done done
            obj.files.bgRatio=[obj.currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
            if exist(obj.files.bgRatio,'file') & ~overwrite
                if nargout==1
                    data=load(obj.files.bgRatio);
                else
                    disp('beta2gamma analysis already exists for this recording');
                end
                return;
            end

            obj.getFilters;
            movWinSamples=movWin/1000*obj.filt.FFs;%obj.filt.FFs in Hz, movWin in samples
            movOLWinSamples=movOLWin/1000*obj.filt.FFs;
            timeBin=(movWin-movOLWin); %ms

            segmentWelchSamples = round(segmentWelch/1000*obj.filt.FFs);
            samplesOLWelch = round(segmentWelchSamples*OLWelch);

            %run welch once to get frequencies for every bin (f) determine frequency bands
            [~,f] = pwelch(randn(1,movWinSamples),segmentWelchSamples,samplesOLWelch,dftPointsWelch,obj.filt.FFs);
            pGamma=find(f>=gammaLow & f<gammaHigh);
            pBeta=find(f>=betaLow & f<betaHigh);

            %if obj.currentDataObj.recordingDuration_ms<movLongWin
            %    movLongWin=obj.currentDataObj.recordingDuration_ms;
            %end

            if win==0
                win=obj.currentDataObj.recordingDuration_ms-tStart;
                endTime=obj.currentDataObj.recordingDuration_ms;
            else
                endTime=min(win+tStart,obj.currentDataObj.recordingDuration_ms);
            end
            startTimes=tStart:(movLongWin-movOLWin):endTime;
            nChunks=numel(startTimes);
            beta2gammaRatioAll=cell(1,nChunks);
            t_ms=cell(1,nChunks);
            %band2to1RatioAllLow=cell(1,nChunks);;band2to1RatioAllHigh=cell(1,nChunks);

            if saveSpectralProfiles
                FMLongB = buffer(true(1,movLongWin/1000*obj.filt.FFs),movWinSamples,movOLWinSamples,'nodelay');
                fftInBuffer=size(FMLongB,2);
                allFreqProfiles=zeros(ceil(dftPointsWelch/2)+1,nChunks*fftInBuffer);
            else
                allFreqProfiles=[];
            end
            if applyNotch
                obj.filt.FN=filterData(obj.currentDataObj.samplingFrequency(1));
                obj.filt.FN.filterDesign='cheby1';
                obj.filt.FN.padding=true;
                obj.filt.FN=obj.filt.FN.designNotch;
            end

            loadCh=ch;
            if ~isempty(avgOnCh)
                loadCh=avgOnCh; %this can not be moved to other positions
            end

            fprintf('\nBeta2Gamma extraction (%d chunks)-',nChunks);
            for i=1:nChunks
                fprintf('%d,',i);
                MLong=obj.currentDataObj.getData(loadCh,startTimes(i),movLongWin);
                if applyNotch
                    MLong=obj.filt.FN.getFilteredData(MLong); %for 50Hz noise
                end
                FMLong=obj.filt.F.getFilteredData(MLong);
                if ~isempty(avgOnCh)
                    FMLong=mean(FMLong,1);
                end

                FMLong(FMLong<-maxVoltage | FMLong>maxVoltage)=nan; %remove high voltage movement artifacts

                FMLongB = buffer(FMLong,movWinSamples,movOLWinSamples,"nodelay");
                pValid=all(~isnan(FMLongB));

                beta2gammaRatioAll{i}=nan(1,numel(pValid)); %changes from zeros to nan in these 3 lines (Mark)
                gammaAll{i}=nan(1,numel(pValid));
                betaAll{i}=nan(1,numel(pValid));
                if any(pValid)
                    [pxx,f] = pwelch(FMLongB(:,pValid),segmentWelchSamples,samplesOLWelch,dftPointsWelch,obj.filt.FFs);
                    beta2gammaRatioAll{i}(pValid)=(mean(pxx(pBeta,:))./mean(pxx(pGamma,:)))';
                    gammaAll{i}(pValid)=mean(pxx(pGamma,:))';
                    betaAll{i}(pValid)=mean(pxx(pBeta,:))';
                else
                    pxx=zeros(dftPointsWelch/2+1,numel(pValid));
                end

                if saveSpectralProfiles
                    allFreqProfiles(:,(fftInBuffer*(i-1)+find(pValid)))=pxx;
                end

                t_ms{i}=startTimes(i)+((movWin/2):timeBin:(movLongWin-movWin/2));
            end
            fprintf('\n');
            beta2gammaRatioAll{end}(t_ms{end}>(endTime-movWin/2))=NaN;
            gammaAll{end}(t_ms{end}>(endTime-movWin/2))=NaN;
            betaAll{end}(t_ms{end}>(endTime-movWin/2))=NaN;

            beta2gammaRatio=cell2mat(beta2gammaRatioAll);beta2gammaRatio=beta2gammaRatio(:);
            gammaBand=cell2mat(gammaAll);gammaBand=gammaBand(:);
            betaBand=cell2mat(betaAll);betaBand=betaBand(:);

            t_ms=cell2mat(t_ms);

            save(obj.files.bgRatio,'t_ms','beta2gammaRatio','par','betaBand','gammaBand','allFreqProfiles');
            if nargout==1
                data=load(obj.files.bgRatio);
            end

               
        end
        
        % %% plot ACCtoDB
        % function [] = plotDBtoACC(obj)
        %     bin = 1000;
        %     convWin = 60;
        %     th= 100;
        % 
        %     [counts, edges] =  histcounts(LM.t_mov_ms,BinWidth=bin);
        %     h_thresh= counts>th;
        %     % h_times = h.BinEdges;
        %     y = ones(1,convWin);
        %     conv = convn(h_thresh,y,'same');
        %     figure; plot(conv);hold on;plot(DB.bufferedDelta2BetaRatio)
        %     figure; scatter(conv(1:length(DB.bufferedDelta2BetaRatio)),DB.bufferedDelta2BetaRatio,'.');
        %     xlabel('Movement');ylabel("D/B")
        % 
        % end 


        %% plotTrialsGB
        function [] = plotTrailsBG(obj,varargin)
          % parameter and settings
            obj.checkFileRecording;
            
            parseObj = inputParser;
            parseObj.FunctionName='wakeAnalysis\plotTrailsGB';
            
            addParameter(parseObj,'ch',obj.recTable.defaulLFPCh(obj.currentPRec),@isnumeric);
            addParameter(parseObj,'saveFigures',1,@isnumeric);
            addParameter(parseObj,'trials',1:20,@isnumeric); %trials to incluse in analysis
            addParameter(parseObj,'win',5000,@isnumeric); % before and after start trial, ms
            addParameter(parseObj,'inputParams',false,@isnumeric);    
            addParameter(parseObj,'fullwin',20000,@isnumeric);    
            addParameter(parseObj,'bug_t',6000,@isnumeric);    

            parseObj.parse(varargin{:});
            if parseObj.Results.inputParams
                disp(parseObj.Results);
                return;
            end
            
            %evaluate all input parameters in workspace
            for i=1:numel(parseObj.Parameters)
                eval([parseObj.Parameters{i} '=' 'parseObj.Results.(parseObj.Parameters{i});']);
            end
        
            bgRatioFile=[obj.currentAnalysisFolder filesep 'bgRatio_ch' num2str(ch) '.mat'];
            obj.checkFileRecording(bgRatioFile,'beta to gamma file missing, please first run getBGRatio');
            load(bgRatioFile); %load data

            arenaCSVFile=[obj.currentAnalysisFolder filesep 'ArenaCSV.mat'];
            obj.checkFileRecording(arenaCSVFile,'arenaCSV file missing, please first run getArenaCSV');
            load(arenaCSVFile); %load data


            zeroTimes = arenaCSVs.startTrigSh-bug_t; %6 seconds before trial started
            endTrialT = arenaCSVs.endTrigSh - zeroTimes;
            strikeTrialNum =arenaCSVs.strikeTrialNum;
            strikeT = arenaCSVs.strikeTrigSh -zeroTimes(strikeTrialNum);
            str_ind = 1;


            matLen = length(find((t_ms>zeroTimes(1))&(t_ms<(zeroTimes(1)+fullwin))));
            BG_mat = zeros([length(zeroTimes),matLen]);
            % BG_mat = [];
            for i=1:length(zeroTimes)
                BG_ind = find((t_ms>zeroTimes(i))&(t_ms<(zeroTimes(i)+fullwin)));
                BG_mc = beta2gammaRatio(BG_ind);

                %     BG_mat= [BG_mat;BG_mc'];
                BG_mat(i,1:length(BG_mc)) = BG_mc';
            end
            
            % plot
            figure;
            %add subplot for the B/G ratio
            ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
            tdiff = (par.movWin-par.movOLWin);
            timeAxis = ((0:matLen-1) * tdiff)/1000; %in s
            % timeAxis = (0:500:(length(BG_mat)-1)*500)/1000; % 0 ms to 19500 ms
            imagesc(timeAxis,1:length(zeroTimes),BG_mat)
            xlabel("Time (s)")
            hold on
            % Plot start time horizontal line:
            xline(bug_t/1000,'r' ,'LineWidth', 1.5, 'Label', 'Bug Apperance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
            ylabel('Trial #')
            % Plot end times and strike times
            for i = 1:size(BG_mat, 1)
                h1 = plot(endTrialT(i)/1000, i, '.', 'Color', 'k','MarkerSize',20); % 'ko' for black circles
                if ismember (i,strikeTrialNum)
                    h2 = plot(strikeT(str_ind)/1000, i, 'pentagram', 'Color', 'g','MarkerFaceColor','g','MarkerSize',8); % 'ko' for black circles
                    str_ind = str_ind+1;
                end
            end
            %add legend and title:
            if ~isempty(strikeTrialNum)
                legend([h2(1), h1(1)], {'Strike','End Trial'},'Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]
            else
                legend(h1(1), 'End Trial','Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]

            end

            %add colorbar
            c = colorbar;
            c.Position = [0.82 0.3 0.02 0.6];
            % Set vertical label for colorbar
            c.Label.String = 'B/G ratio'; % Set label text
            c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
            c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

            % Plot Average
            ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
            plot(timeAxis, mean(BG_mat,1,"omitmissing"),'LineWidth', 1.5)
            xlabel('Time [s]'); ylabel('avg. B2G')
            % ylim([10 85])
            xline(bug_t/1000,'r' ,'LineWidth', 1.5)
            linkaxes([ax1,ax2],'x');

            % Save plot:
            set(gcf,'PaperPosition',[.25 3 8 6])
            saveas(gcf,strcat(obj.currentPlotFolder, '/TrialsBG.pdf'));


        
        end


        

    end

end