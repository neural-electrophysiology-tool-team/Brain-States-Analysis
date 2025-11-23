

function f = plotStimSham(SA)
    % SA is an instance of sleep analysis class,with a record currently
    % selected
    stimShamFile=[SA.currentAnalysisFolder filesep  'stimSham.mat'];
    SA.checkFileRecording(stimShamFile,'stim Sham file missing, please first run getStimSham');
    load(stimShamFile); %load data

    trialsham = height(StimDBSham)-70+1;

    colorLim=[0 600];
    f=figure;
    subplot(4,2,[1:2:6]);imagesc(StimDBSham(trialsham:end,:),colorLim);ylabel('Trial #');title('Sham');hold on;set(gca,'XTick',[]);
    cb=colorbar('Position',[0.47 0.76 0.013 0.17]);ylabel(cb,'\delta/\beta');
    line([pre/1000 pre/1000],ylim,'color','r');
    line([(pre+stimDuration)/1000 (pre+stimDuration)/1000],ylim,'color','r');
    subplot(4,2,7);plot(ts-pre/1000,mean(StimDBSham,'omitnan'),'k');xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/2);
    line([0 0],ylim,'color','r');xlim([-pre/1000 post/1000])
    line([stimDuration/1000 stimDuration/1000],ylim,'color','r');
    subplot(4,2,[2:2:6]);imagesc(StimDB,colorLim);ylabel('Trial #');title('Stim');set(gca,'XTick',[]);
    cb=colorbar('Position',[ 0.91 0.76 0.013 0.17]);ylabel(cb,'\delta/\beta');
    line([pre/1000 pre/1000],ylim,'color','r');
    line([(pre+stimDuration)/1000 (pre+stimDuration)/1000],ylim,'color','r');
    subplot(4,2,8);plot(ts-pre/1000,mean(StimDB,'omitnan'),'k');xlabel(['Time [s]']);ylabel('Avg.');ylim(colorLim/2);
    xlim([-pre/1000 post/1000])
    line([0 0],ylim,'color','r');
    line([stimDuration/1000 stimDuration/1000],ylim,'color','r');

%   save?
    fileName=[SA.currentPlotFolder filesep 'stim_sham_activation.pdf'];
    saveas (f, fileName);
    % saveas (gcf, [analysisFolder filesep 'lizMovAllNights_norm.pdf']);


end
% 
% function plotStimSham(SA, axArray)
%     % SA is an instance of sleep analysis class,with a record currently
%     % selected
%     stimShamFile = [SA.currentAnalysisFolder filesep 'stimSham.mat'];
%     SA.checkFileRecording(stimShamFile, ...
%         'stim Sham file missing, please first run getStimSham');
%     load(stimShamFile); % loads StimDBSham, StimDB, ts, pre, stimDuration
% 
%     trialsham = height(StimDBSham) - 70 + 1;
%     colorLim = [0 250];
% 
%     % If no axes array was given, use nexttile to create them
%     if nargin < 2 || isempty(axArray)
%         axArray(1) = nexttile;
%         axArray(2) = nexttile;
%         axArray(3) = nexttile;
%         axArray(4) = nexttile;
%     end
% 
%     % --- Sham heatmap
%     axes(axArray(1));
%     imagesc(StimDBSham(trialsham:end,:), colorLim);
%     ylabel('Trial #'); title('Sham'); hold on; set(gca,'XTick',[]);
%     line([pre/1000 pre/1000], ylim, 'color','r');
%     cb = colorbar; ylabel(cb,'\delta/\beta');
% 
%     % --- Sham average trace
%     axes(axArray(2));
%     plot(ts-pre/1000, nanmean(StimDBSham));
%     xlabel('Time [s]'); ylabel('Avg.'); ylim(colorLim/3);
%     line([0 0], ylim,'color','r');
%     line([stimDuration/1000 stimDuration/1000], ylim,'color','r');
% 
%     % --- Stim heatmap
%     axes(axArray(3));
%     imagesc(StimDB, colorLim);
%     ylabel('Trial #'); title('Stim'); set(gca,'XTick',[]);
%     line([pre/1000 pre/1000], ylim,'color','r');
%     cb = colorbar; ylabel(cb,'\delta/\beta');
% 
%     % --- Stim average trace
%     axes(axArray(4));
%     plot(ts-pre/1000, nanmean(StimDB));
%     xlabel('Time [s]'); ylabel('Avg.'); ylim(colorLim/3);
%     line([0 0], ylim,'color','r');
%     line([stimDuration/1000 stimDuration/1000], ylim,'color','r');
% end
