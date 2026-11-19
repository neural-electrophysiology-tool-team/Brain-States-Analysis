function data=getBeta2GammaRatioOfir(obj,varargin)
    obj.checkFileRecording;

    parseObj = inputParser;
    addParameter(parseObj,'ch',obj.recTable.defaulLFPCh(obj.currentPRec),@isnumeric); %1:32   obj.recTable.defaulLFPCh(obj.currentPRec)
    addParameter(parseObj,'avgOnCh',[],@isnumeric); %uses several averaged channels for d/b extraction
    addParameter(parseObj,'movLongWin',1000*60*30,@isnumeric); %max freq. to examine
    addParameter(parseObj,'movWin',2000,@isnumeric); %2000 %500
    addParameter(parseObj,'movOLWin',1900,@isnumeric); %1000 %448
    addParameter(parseObj,'segmentWelch',1000,@isnumeric); %248
    addParameter(parseObj,'dftPointsWelch',2^10,@isnumeric);
    addParameter(parseObj,'OLWelch',0.5);
    % addParameter(parseObj,'tStart',event_oe_time-6000,@isnumeric); %window begins 1 sec prior to the sound event in the blockLog. If 0 begins from the start of the recording
    % addParameter(parseObj,'win',12000,@isnumeric); %if 0 uses the whole recording duration
    addParameter(parseObj,'tStart',0,@isnumeric);
    addParameter(parseObj,'win',0,@isnumeric); %if 0 uses the whole recording duration
    addParameter(parseObj,'gammaLow',60,@isnumeric); %60 %21 %31
    addParameter(parseObj,'gammaHigh',80,@isnumeric); %80 %30 %50
    addParameter(parseObj,'betaLow',10,@isnumeric);
    addParameter(parseObj,'betaHigh',30,@isnumeric); %30 %20
    addParameter(parseObj,'applyNotch',0,@isnumeric);
    addParameter(parseObj,'saveSpectralProfiles',1,@isnumeric);
    addParameter(parseObj,'maxVoltage',1000,@isnumeric);
    addParameter(parseObj,'overwrite',1,@isnumeric);
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
    % if exist(obj.files.bgRatio,'file') & ~overwrite
    %     if nargout==1
    %         data=load(obj.files.bgRatio);
    %     else
    %         disp('beta2gamma analysis already exists for this recording');
    %     end
    %     return;
    % end

    obj.getFilters;
    movWinSamples=movWin/1000*obj.filt.FFs;%obj.filt.FFs in Hz, movWin in samples
    movOLWinSamples=movOLWin/1000*obj.filt.FFs;
    timeBin=(movWin-movOLWin); %ms

    segmentWelchSamples = round(segmentWelch/1000*obj.filt.FFs);
    samplesOLWelch = round(segmentWelchSamples*OLWelch);

    %run welch once to get frequencies for every bin (f) determine frequency bands
    [~,f] = pwelch(randn(1,movWinSamples),segmentWelchSamples,samplesOLWelch,dftPointsWelch,obj.filt.FFs);
    % if exist(obj.files.bgRatio,'file') & ~overwrite
    %     if nargout==1
    %         data=load(obj.files.bgRatio);
    %     else
            % disp('beta2gamma analysis already exists for this recording');
    %     end
    %     return;
    % end    
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
        
        % MLongAll=obj.currentDataObj.getData(loadCh,startTimes(i),movLongWin);
        % MLongAll=squeeze(MLongAll);
        % MLongAll=bsxfun(@minus,MLongAll,median(MLongAll,1));
        % 
        % MLongOne=MLongAll(25,:,:);
        % MLong=reshape(MLongOne,1,size(MLongOne,1),size(MLongOne,2));

        MLong=obj.currentDataObj.getData(loadCh,startTimes(i),movLongWin);
        % MLong=signal_s_allSortRes(:,1);

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

    save(obj.files.bgRatio,'t_ms','beta2gammaRatio','par','betaBand','gammaBand','allFreqProfiles','f');
    if nargout==1
        data=load(obj.files.bgRatio);
    end

    a=1;

       
end