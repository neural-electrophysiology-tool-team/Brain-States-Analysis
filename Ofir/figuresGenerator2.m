%% 1. Plot ratios of frequencies per animal and hunter-type

animal = '228';

hunter = hunter228; %current animal

isCircle = hunter.Type == "circle" | hunter.Type == "circle_flip";

hunter_circle = hunter(isCircle, :);
hunter_other  = hunter(~isCircle, :);

T = hunter_circle;

idx_strikes = find(~isnan(T.Strike_ms));
idx_rests = find(isnan(T.Strike_ms)); % no strikes

idx_2plot = idx_strikes;

%% 2. Hunter - bug on: all
% Bug ON: b/g
rows2remove=[];
BG=cell2mat(T.b2g(idx_2plot,:)); endTrialT=T.End_ms(idx_2plot,:); strikeT=T.Strike_ms(idx_2plot,:);
% rows = all(~isnan(BG(:,51:70)), 2); % ignore trials with nans after bug appearance
% BG_2save=BG(rows,:); endTrialT_2plot=endTrialT(rows,:); strikeT_2plot=strikeT(rows,:);

BG_2save=BG; endTrialT_2plot=endTrialT; strikeT_2plot=strikeT;
% Remove circle rest trials shorter than 40 s or longer than 43 s
 % endTrialT_2plot=endTrialT; strikeT_2plot=strikeT;
% short_long_circle_trials = find(endTrialT_2plot<4e4 | endTrialT_2plot>4.3e4); %4   4.3     1.6    2.0 for PV106 

short_long_circle_trials = find(endTrialT_2plot>5e4 | strikeT_2plot>5e4 | strikeT_2plot<0e4); % Remove larger than 5 s non-circle trials 0.5
rows2remove = short_long_circle_trials;
BG_2save(rows2remove,:)=[]; endTrialT_2plot(rows2remove,:)=[]; strikeT_2plot(rows2remove,:)=[];

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2save,2)-1) * tdiff)-win)/1000; %in s

figure(1);
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2save,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis,1:size(BG_2save),BG_2save)
axis xy % important: row 1 at top
xlabel("Time from Bug Appearance (s)")
hold on
title("PV208 No strikes, low horizontal noise");
% Plot start time horizontal line:
xline(0,'g' ,'LineWidth', 1.5, 'Label', 'Bug Apperance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
% Plot end times and strike times
for i = 1:size(BG_2save, 1)
    h1 = plot(endTrialT_2plot(i)/1000, i, '.', 'Color', 'r','MarkerSize',20); % 'ko' for black circles
    % if ismember (i,strikeTrialNum)
        h2 = plot(strikeT_2plot(i)/1000, i, 'pentagram', 'Color', 'k','MarkerFaceColor','k','MarkerSize',8); % 'ko' for black circles
        % str_ind = str_ind+1;
    % end
end
%add legend and title:
 % if ~isempty(strikeTrialNum)
    legend([h2(1), h1(1)], {'Strike','End Trial'},'Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]
% else
    % legend(h1(1), 'End Trial','Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]

% end

thr = prctile(BG_2save(:),95);
clim([min(BG_2save(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'B/G ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, nanmean(BG_2save,1),'LineWidth', 1.5)
xlabel('Time from Bug Appearance (s)'); ylabel('avg. B2G')
% ylim([10 85])
xline(0,'g' ,'LineWidth', 1.5)
linkaxes([ax1,ax2],'x');

% Hunter - Bug ON: d/b

BG=cell2mat(T.d2b(idx_2plot,:)); endTrialT=T.End_ms(idx_2plot,:); strikeT=T.Strike_ms(idx_2plot,:);
rows = all(~isnan(BG(:,51:70)), 2); % ignore trials with nans after bug appearance
% BG_2plot=BG(rows,:); endTrialT_2plot=endTrialT(rows,:); strikeT_2plot=strikeT(rows,:);
BG_2plot=BG; endTrialT_2plot=endTrialT; strikeT_2plot=strikeT;

% Remove circle rest trials shorter than 40 s or longer than 43 s
% short_long_circle_trials = find(endTrialT_2plot<3e4 | endTrialT_2plot>3.5e4); %4   4.3     1.6    2.0

% short_long_circle_trials = find(endTrialT_2plot>0.5e4); % Remove larger than 5 s non-circle trials
% rows2remove = short_long_circle_trials;
BG_2plot(rows2remove,:)=[]; endTrialT_2plot(rows2remove,:)=[]; strikeT_2plot(rows2remove,:)=[];

figure(2);
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2plot,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis,1:size(BG_2plot),BG_2plot)
axis xy % important: row 1 at top
xlabel("Time from Bug Appearance (s)")
hold on
title("PV208 No strikes, low horizontal noise");
% Plot start time horizontal line:
xline(0,'g' ,'LineWidth', 1.5, 'Label', 'Bug Apperance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
% Plot end times and strike times
for i = 1:size(BG_2plot, 1)
    h1 = plot(endTrialT_2plot(i)/1000, i, '.', 'Color', 'r','MarkerSize',20); % 'ko' for black circles
    % if ismember (i,strikeTrialNum)
        h2 = plot(strikeT_2plot(i)/1000, i, 'pentagram', 'Color', 'k','MarkerFaceColor','k','MarkerSize',8); % 'ko' for black circles
        % str_ind = str_ind+1;
    % end
end
%add legend and title:
 % if ~isempty(strikeTrialNum)
    legend([h2(1), h1(1)], {'Strike','End Trial'},'Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]
% else
    % legend(h1(1), 'End Trial','Position',[0.85 0.15 0.04 0.08],'Box','off');%[x, y, width, height]

% end

thr = prctile(BG_2plot(:),85);
clim([min(BG_2plot(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'D/B ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, nanmean(BG_2plot,1),'LineWidth', 1.5)
xlabel('Time from Bug Appearance (s)'); ylabel('avg. D2B')
% ylim([10 85])
xline(0,'g' ,'LineWidth', 1.5)
linkaxes([ax1,ax2],'x');


% power of frequency bands


% animal = "PV228";
% isAnimal = mat.Animal == animal;
% 
% T = mat(isAnimal, :);

mat_delta = cell2mat(T.d_band); mat_beta = cell2mat(T.b_band); mat_gamma = cell2mat(T.g_band);
mat_delta = mat_delta (idx_2plot,:); mat_beta =mat_beta(idx_2plot,:); mat_gamma =  mat_gamma(idx_2plot,:);
% rows = all(~isnan(BG(:,51:70)), 2); % ignore trials with nans after bug appearance
% mat_delta = mat_delta (rows,:); mat_beta =mat_beta(rows,:); mat_gamma =  mat_gamma(rows,:);
mat_delta(rows2remove,:)=[]; mat_beta(rows2remove,:)=[]; mat_gamma(rows2remove,:)=[];


mean_band_d = nanmean(mat_delta);
mean_band_b = nanmean(mat_beta);
mean_band_g = nanmean(mat_gamma);

% mean_band_d = nanmean(cell2mat(T.d_band));
% mean_band_b = nanmean(cell2mat(T.b_band));
% mean_band_g = nanmean(cell2mat(T.g_band));

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);
timeAxis = (((0:length(mean_band_d)-1) * tdiff)-win)/1000; %in s

% fs1 = 10;               % samples per second
% n = 240;  %120 sound old %460 screen diode   240 sound new   600 hunter      
% time1 = (0:n-1)/fs1;              
% time1_centered = time1 - 12;     % -5 hunter




figure(3);
subplot(3,1,1);
plot(timeAxis,mean_band_d,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Sound Onset (s)"); 
ylabel('Power');
% ylim([0 500]);
set(gca,'FontSize',18)
title('PV143: Delta = 0-5 Hz','FontSize',18);

subplot(3,1,2);
plot(timeAxis,mean_band_b,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Meadow Movement (s)");
ylabel('Power');
title('Beta = 10-30 Hz','FontSize',18);
set(gca,'FontSize',18)

subplot(3,1,3);
plot(timeAxis,mean_band_g,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
xlabel("Time from Bug Appearance (s)");  ylabel('Power'); %Time from Meadow Movement   Bug Appearance
title('Gamma = 60-80 Hz','FontSize',18);
set(gca,'FontSize',18)



%% 3. Hunter - Bug OFF: all
BG_2plot_off=[]; BG_2save_off=[]; mat_delta_off=[]; mat_beta_off=[]; mat_gamma_off=[];
for i=1:length(endTrialT_2plot)
    end_ind = find(timeAxis <= endTrialT_2plot(i)/1000, 1, 'last');
    % BG_2plot_off = [BG_2plot_off; BG_2plot(end_ind-50:end_ind+49)];
    BG_mc_save = BG_2save(i,end_ind-50:end_ind+49);
    BG_2save_off(i,1:length(BG_mc_save)) = BG_mc_save;
    BG_mc_plot = BG_2plot(i,end_ind-50:end_ind+49);
    BG_2plot_off(i,1:length(BG_mc_plot)) = BG_mc_plot;
    mat_delta_mc_off = mat_delta(i,end_ind-50:end_ind+49);
    mat_delta_off(i,1:length(mat_delta_mc_off)) = mat_delta_mc_off;
    mat_beta_mc_off = mat_beta(i,end_ind-50:end_ind+49);
    mat_beta_off(i,1:length(mat_beta_mc_off)) = mat_beta_mc_off;
    mat_gamma_mc_off = mat_gamma(i,end_ind-50:end_ind+49);
    mat_gamma_off(i,1:length(mat_gamma_mc_off)) = mat_gamma_mc_off;
end


figure(4);
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis1 = (((0:size(BG_2save_off,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis1,1:size(BG_2save_off),BG_2save_off)
axis xy % important: row 1 at top
xlabel("Time from Bug Disapperance (s)")
hold on
title("PV106 No strikes, circles");
% Plot start time horizontal line:
xline(0,'r' ,'LineWidth', 1.5, 'Label', 'Bug Disapperance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
xticks(-5:1:5);

thr = prctile(BG_2save_off(:),95);
clim([min(BG_2save_off(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'B/G ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis1, nanmean(BG_2save_off,1),'LineWidth', 1.5)
xlabel('Time from Bug Disapperance (s)'); ylabel('avg. B2G')
% ylim([10 85])
xline(0,'r' ,'LineWidth', 1.5)
xticks(-5:1:5);
linkaxes([ax1,ax2],'x');


figure(5);
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis1 = (((0:size(BG_2plot_off,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis1,1:size(BG_2plot_off),BG_2plot_off)
axis xy % important: row 1 at top
xlabel("Time from Bug Disapperance (s)")
hold on
title("PV106 No strikes, circles");
% Plot start time horizontal line:
xline(0,'r' ,'LineWidth', 1.5, 'Label', 'Bug Disapperance', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
xticks(-5:1:5);

thr = prctile(BG_2plot_off(:),85);
clim([min(BG_2plot_off(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'D/B ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis1, nanmean(BG_2plot_off,1),'LineWidth', 1.5)
xlabel('Time from Bug Disapperance (s)'); ylabel('avg. D2B')
% ylim([10 85])
xline(0,'r' ,'LineWidth', 1.5)
xticks(-5:1:5);
linkaxes([ax1,ax2],'x');



mean_band_d_off = nanmean(mat_delta_off);
mean_band_b_off = nanmean(mat_beta_off);
mean_band_g_off = nanmean(mat_gamma_off);

% mean_band_d = nanmean(cell2mat(T.d_band));
% mean_band_b = nanmean(cell2mat(T.b_band));
% mean_band_g = nanmean(cell2mat(T.g_band));

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);
timeAxis1 = (((0:length(mean_band_d_off)-1) * tdiff)-win)/1000; %in s

% fs1 = 10;               % samples per second
% n = 240;  %120 sound old %460 screen diode   240 sound new   600 hunter      
% time1 = (0:n-1)/fs1;              
% time1_centered = time1 - 12;     % -5 hunter



figure(6);
subplot(3,1,1);
plot(timeAxis1,mean_band_d_off,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'r','LineWidth',2 , 'Label', 'Bug Off', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Sound Onset (s)"); 
ylabel('Power');
% ylim([0 500]);
set(gca,'FontSize',18)
title('PV143: Delta = 0-5 Hz','FontSize',18);

subplot(3,1,2);
plot(timeAxis1,mean_band_b_off,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'r','LineWidth',2 , 'Label', 'Bug Off', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Meadow Movement (s)");
ylabel('Power');
title('Beta = 10-30 Hz','FontSize',18);
set(gca,'FontSize',18)

subplot(3,1,3);
plot(timeAxis1,mean_band_g_off,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'r','LineWidth',2 , 'Label', 'Bug Off', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
xlabel("Time from Bug Disappearance (s)");  ylabel('Power'); %Time from Meadow Movement   Bug Appearance
title('Gamma = 60-80 Hz','FontSize',18);
set(gca,'FontSize',18)




%% 4. Hunter - strike time: all
BG_2plot_strike=[]; BG_2save_strike=[]; mat_delta_strike=[]; mat_beta_strike=[]; mat_gamma_strike=[];
for i=1:length(strikeT_2plot)
    strike_ind = find(timeAxis >= strikeT_2plot(i)/1000, 1, 'first');
    % BG_2plot_off = [BG_2plot_off; BG_2plot(end_ind-50:end_ind+49)];
    BG_mc_save = BG_2save(i,strike_ind-50:strike_ind+49);
    BG_2save_strike(i,1:length(BG_mc_save)) = BG_mc_save;
    BG_mc_plot = BG_2plot(i,strike_ind-50:strike_ind+49);
    BG_2plot_strike(i,1:length(BG_mc_plot)) = BG_mc_plot;
    mat_delta_mc_strike = mat_delta(i,strike_ind-50:strike_ind+49);
    mat_delta_strike(i,1:length(mat_delta_mc_strike)) = mat_delta_mc_strike;
    mat_beta_mc_strike = mat_beta(i,strike_ind-50:strike_ind+49);
    mat_beta_strike(i,1:length(mat_beta_mc_strike)) = mat_beta_mc_strike;
    mat_gamma_mc_strike = mat_gamma(i,strike_ind-50:strike_ind+49);
    mat_gamma_strike(i,1:length(mat_gamma_mc_strike)) = mat_gamma_mc_strike;
end


figure(7);
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis2 = (((0:size(BG_2save_strike,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis2,1:size(BG_2save_strike),BG_2save_strike)
axis xy % important: row 1 at top
xlabel("Time from Bug Disapperance (s)")
hold on
title("PV106 No strikes, circles");
% Plot start time horizontal line:
xline(0,'k' ,'LineWidth', 1.5, 'Label', 'Strike', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
xticks(-5:1:5);

thr = prctile(BG_2save_strike(:),95);
clim([min(BG_2save_strike(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'B/G ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis2, nanmean(BG_2save_strike,1),'LineWidth', 1.5)
xlabel('Time from Strike (s)'); ylabel('avg. B2G')
% ylim([10 85])
xline(0,'k' ,'LineWidth', 1.5)
xticks(-5:1:5);
linkaxes([ax1,ax2],'x');


figure(8);
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis1 = (((0:size(BG_2plot_strike,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis2,1:size(BG_2plot_strike),BG_2plot_strike)
axis xy % important: row 1 at top
xlabel("Time from Bug Disapperance (s)")
hold on
title("PV106 No strikes, circles");
% Plot start time horizontal line:
xline(0,'k' ,'LineWidth', 1.5, 'Label', 'Strike', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
xticks(-5:1:5);

thr = prctile(BG_2plot_strike(:),85);
clim([min(BG_2plot_strike(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'D/B ratio'; % Set label text
c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis2, nanmean(BG_2plot_strike,1),'LineWidth', 1.5)
xlabel('Time from Strike (s)'); ylabel('avg. D2B')
% ylim([10 85])
xline(0,'k' ,'LineWidth', 1.5)
xticks(-5:1:5);
linkaxes([ax1,ax2],'x');



mean_band_d_strike = nanmean(mat_delta_strike);
mean_band_b_strike = nanmean(mat_beta_strike);
mean_band_g_strike = nanmean(mat_gamma_strike);

% mean_band_d = nanmean(cell2mat(T.d_band));
% mean_band_b = nanmean(cell2mat(T.b_band));
% mean_band_g = nanmean(cell2mat(T.g_band));

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);
timeAxis2 = (((0:length(mean_band_d_strike)-1) * tdiff)-win)/1000; %in s

% fs1 = 10;               % samples per second
% n = 240;  %120 sound old %460 screen diode   240 sound new   600 hunter      
% time1 = (0:n-1)/fs1;              
% time1_centered = time1 - 12;     % -5 hunter



figure(9);
subplot(3,1,1);
plot(timeAxis2,mean_band_d_strike,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'k','LineWidth',2 , 'Label', 'Strike', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Sound Onset (s)"); 
ylabel('Power');
% ylim([0 500]);
set(gca,'FontSize',18)
title('PV143: Delta = 0-5 Hz','FontSize',18);

subplot(3,1,2);
plot(timeAxis2,mean_band_b_strike,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'k','LineWidth',2 , 'Label', 'Strike', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Meadow Movement (s)");
ylabel('Power');
title('Beta = 10-30 Hz','FontSize',18);
set(gca,'FontSize',18)

subplot(3,1,3);
plot(timeAxis2,mean_band_g_strike,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'k','LineWidth',2 , 'Label', 'Strike', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
xlabel("Time from Strike (s)");  ylabel('Power'); %Time from Meadow Movement   Bug Appearance
title('Gamma = 60-80 Hz','FontSize',18);
set(gca,'FontSize',18)


%% 5. Hunter - all: Getting the matrices

dataName = sprintf('hunter_%s', animal);

data = struct();

data.(sprintf('BG_hunter_%s_others_rests', animal)) = BG_2save;
data.(sprintf('DB_hunter_%s_others_rests', animal)) = BG_2plot;
data.(sprintf('delta_hunter_%s_others_rests', animal)) = mat_delta;
data.(sprintf('beta_hunter_%s_others_rests', animal)) = mat_beta;
data.(sprintf('gamma_hunter_%s_others_rests', animal)) = mat_gamma;

data.(sprintf('BG_hunter_%s_others_rests_off', animal)) = BG_2save_off;
data.(sprintf('DB_hunter_%s_others_rests_off', animal)) = BG_2plot_off;
data.(sprintf('delta_hunter_%s_others_rests_off', animal)) = mat_delta_off;
data.(sprintf('beta_hunter_%s_others_rests_off', animal)) = mat_beta_off;
data.(sprintf('gamma_hunter_%s_others_rests_off', animal)) = mat_gamma_off;



% data.(sprintf('BG_hunter_%s_others_strikes_strike', animal)) = BG_2save_strike;
% data.(sprintf('DB_hunter_%s_others_strikes_strike', animal)) = BG_2plot_strike;
% data.(sprintf('delta_hunter_%s_others_strikes_strike', animal)) = mat_delta_strike;
% data.(sprintf('beta_hunter_%s_others_strikes_strike', animal)) = mat_beta_strike;
% data.(sprintf('gamma_hunter_%s_others_strikes_strike', animal)) = mat_gamma_strike;



eval(sprintf('%s = data;', dataName));

save('mats4stat.mat', dataName, '-append');

%% 6. Hunter - all: Getting the matrices - version 2

dataName = sprintf('hunter_%s', animal);

% Load the existing structure
S = load('mats4stat.mat', dataName);
data = S.(dataName);

% Add new variables
data.(sprintf('BG_hunter_%s_circles_strikes', animal)) = BG_2save;
data.(sprintf('DB_hunter_%s_circles_strikes', animal)) = BG_2plot;
data.(sprintf('delta_hunter_%s_circles_strikes', animal)) = mat_delta;
data.(sprintf('beta_hunter_%s_circles_strikes', animal)) = mat_beta;
data.(sprintf('gamma_hunter_%s_circles_strikes', animal)) = mat_gamma;

data.(sprintf('BG_hunter_%s_circles_strikes_off', animal)) = BG_2save_off;
data.(sprintf('DB_hunter_%s_circles_strikes_off', animal)) = BG_2plot_off;
data.(sprintf('delta_hunter_%s_circles_strikes_off', animal)) = mat_delta_off;
data.(sprintf('beta_hunter_%s_circles_strikes_off', animal)) = mat_beta_off;
data.(sprintf('gamma_hunter_%s_circles_strikes_off', animal)) = mat_gamma_off;



data.(sprintf('BG_hunter_%s_circles_strikes_strike', animal)) = BG_2save_strike;
data.(sprintf('DB_hunter_%s_circles_strikes_strike', animal)) = BG_2plot_strike;
data.(sprintf('delta_hunter_%s_circles_strikes_strike', animal)) = mat_delta_strike;
data.(sprintf('beta_hunter_%s_circles_strikes_strike', animal)) = mat_beta_strike;
data.(sprintf('gamma_hunter_%s_circles_strikes_strike', animal)) = mat_gamma_strike;


% Save the updated structure
eval(sprintf('%s = data;', dataName));
save('mats4stat.mat', dataName, '-append');


%% 1. Sound analysis: Plot ratios of frequencies per animal

mat = neutral_sound;
animal = "PV228";
isAnimal = mat.Animal == animal;

T = mat(isAnimal, :);

BG_2save=cell2mat(T.b2g);

figure;
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2save,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis,1:size(BG_2save),BG_2save)
axis xy % important: row 1 at top
xlabel("Time from Sound Onset (s)")
hold on
title(animal);
% Plot start time horizontal line:
xline(0,'r' ,'LineWidth', 1.5, 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
xticks(-5:1:5);

thr = prctile(BG_2save(:),95); %95
clim([min(BG_2save(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'B/G ratio'; % Set label text
% c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, nanmean(BG_2save,1),'LineWidth', 1.5)
xlabel('Time from Neutral Sound Onset (s)'); ylabel('avg. B2G')
% ylim([10 85])
xline(0,'r' ,'LineWidth', 1.5)
xticks(-5:1:5);
linkaxes([ax1,ax2],'x');



BG_2plot=cell2mat(T.d2b);

figure;
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2plot,2)-1) * tdiff)-win)/1000; %in s
imagesc(timeAxis,1:size(BG_2plot),BG_2plot)
axis xy % important: row 1 at top
xlabel("Time from Sound Onset (s)")
hold on
title(animal);
% Plot start time horizontal line:
xline(0,'r' ,'LineWidth', 1.5, 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
xticks(-5:1:5);

thr = prctile(BG_2plot(:),85); %95
clim([min(BG_2plot(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'D/B ratio'; % Set label text
% c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, nanmean(BG_2plot,1),'LineWidth', 1.5)
xlabel('Time from Neutral Sound Onset (s)'); ylabel('avg. D2B')
% ylim([10 85])
xline(0,'r' ,'LineWidth', 1.5)
xticks(-5:1:5);
linkaxes([ax1,ax2],'x');

%% 2. Plot power of frequency bands


animal = "PV228";
isAnimal = mat.Animal == animal;

T = mat(isAnimal, :);


mean_band_d = nanmean(cell2mat(T.d_band));
mean_band_b = nanmean(cell2mat(T.b_band));
mean_band_g = nanmean(cell2mat(T.g_band));

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);
timeAxis = (((0:length(mean_band_d)-1) * tdiff)-win)/1000; %in s

% fs1 = 10;               % samples per second
% n = 240;  %120 sound old %460 screen diode   240 sound new   600 hunter      
% time1 = (0:n-1)/fs1;              
% time1_centered = time1 - 12;     % -5 hunter




figure;
subplot(3,1,1);
plot(timeAxis,mean_band_d,'k','LineWidth',2); 
xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Sound Onset (s)"); 
ylabel('Power');
% ylim([0 500]);
set(gca,'FontSize',18)
title('PV228: Delta = 0-5 Hz','FontSize',18);

subplot(3,1,2);
plot(timeAxis,mean_band_b,'k','LineWidth',2); 
xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Meadow Movement (s)");
ylabel('Power');
title('Beta = 10-30 Hz','FontSize',18);
set(gca,'FontSize',18)

subplot(3,1,3);
plot(timeAxis,mean_band_g,'k','LineWidth',2); 
xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
xlabel("Time from Sound Onset (s)");  ylabel('Power'); %Time from Meadow Movement   Bug Appearance
title('Gamma = 60-80 Hz','FontSize',18);
set(gca,'FontSize',18)

%% 3. Getting the matrices

BG_neutral_sound_228 = BG_2save;
DB_neutral_sound_228 = BG_2plot;
delta_neutral_sound_228 = cell2mat(T.d_band);
beta_neutral_sound_228 = cell2mat(T.b_band);
gamma_neutral_sound_228 = cell2mat(T.g_band);

save('mats4stat.mat', ...
    'BG_neutral_sound_228', ...
    'DB_neutral_sound_228', ...
    'delta_neutral_sound_228', ...
    'beta_neutral_sound_228', ...
    'gamma_neutral_sound_228', ...
            '-append');

%% 1. Meadow analysis - stat: Plot ratios of frequencies per animal

T = reg_meadow228;
animal = "PV228";

BG_2save=cell2mat(T.b2g);


% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2save,2)-1) * tdiff)-win)/1000; %in s

figure;
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
imagesc(timeAxis,1:size(BG_2save),BG_2save)
axis xy % important: row 1 at top
xlabel("Time from Meadow Static (s)")
hold on
title(animal);
% Plot start time horizontal line:
xline(0,'r' ,'LineWidth', 1.5, 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
% xticks(-10:1:20);

% Plot meadow movement onset times
for i = 1:size(BG_2save, 1)
    h1 = plot((T.Move_ms(i)-T.Triger_ms(i)-win)/1000, i, '.', 'Color', 'k','MarkerSize',20); % 'ko' for black circles
end
%add legend and title:
legend([h1(1)], {'Meadow movement'},'Position',[0.89 0.15 0.04 0.08],'Box','off');%[x, y, width, height]


thr = prctile(BG_2save(:),95); %95
clim([min(BG_2save(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'B/G ratio'; % Set label text
% c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, nanmean(BG_2save,1),'LineWidth', 1.5)
xlabel('Time from Meadow Static (s)'); ylabel('avg. B2G')
% ylim([10 85])
xline(0,'r' ,'LineWidth', 1.5)
% xticks(-10:1:20);
linkaxes([ax1,ax2],'x');



BG_2plot=cell2mat(T.d2b);


% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2plot,2)-1) * tdiff)-win)/1000; %in s


figure;
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
imagesc(timeAxis,1:size(BG_2plot),BG_2plot)
axis xy % important: row 1 at top
xlabel("Time from Meadow Static (s)")
hold on
title(animal);
% Plot start time horizontal line:
xline(0,'r' ,'LineWidth', 1.5, 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
% xticks(-10:1:20);

% Plot meadow movement onset times
for i = 1:size(BG_2plot, 1)
    h1 = plot((T.Move_ms(i)-T.Triger_ms(i)-win)/1000, i, '.', 'Color', 'k','MarkerSize',20); % 'ko' for black circles
end
%add legend and title:
legend([h1(1)], {'Meadow movement'},'Position',[0.89 0.15 0.04 0.08],'Box','off');%[x, y, width, height]


thr = prctile(BG_2plot(:),85); %95
clim([min(BG_2plot(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'D/B ratio'; % Set label text
% c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis, nanmean(BG_2plot,1),'LineWidth', 1.5)
xlabel('Time from Meadow Static (s)'); ylabel('avg. D2B')
% ylim([10 85])
xline(0,'r' ,'LineWidth', 1.5)
% xticks(-10:1:20);
linkaxes([ax1,ax2],'x');

% Plot power of frequency bands
mat_delta = cell2mat(T.d_band);
mat_beta = cell2mat(T.b_band);
mat_gamma = cell2mat(T.g_band);

mean_band_d = nanmean(mat_delta);
mean_band_b = nanmean(mat_beta);
mean_band_g = nanmean(mat_gamma);

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);
timeAxis = (((0:length(mean_band_d)-1) * tdiff)-win)/1000; %in s

% fs1 = 10;               % samples per second
% n = 240;  %120 sound old %460 screen diode   240 sound new   600 hunter      
% time1 = (0:n-1)/fs1;              
% time1_centered = time1 - 12;     % -5 hunter




figure;
subplot(3,1,1);
plot(timeAxis,mean_band_d,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Sound Onset (s)"); 
ylabel('Power');
% ylim([0 500]);
set(gca,'FontSize',18)
title('PV228 white: Delta = 0-5 Hz','FontSize',18);

subplot(3,1,2);
plot(timeAxis,mean_band_b,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Meadow Movement (s)");
ylabel('Power');
title('Beta = 10-30 Hz','FontSize',18);
set(gca,'FontSize',18)

subplot(3,1,3);
plot(timeAxis,mean_band_g,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(20,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
xlabel("Time from Meadow Static (s)");  ylabel('Power'); %Time from Meadow Movement   Bug Appearance
title('Gamma = 60-80 Hz','FontSize',18);
set(gca,'FontSize',18)

%% 2. Meadow analysis - move: Plot ratios of frequencies per animal

% default from dataOrganizer
movWin = 2000; 
movOLWin = 1900; 
win = 5000; % for hunter

tdiff = (movWin-movOLWin);

timeAxis = (((0:size(BG_2save,2)-1) * tdiff)-win)/1000; %in s

BG_2plot_move=[]; BG_2save_move=[]; mat_delta_move=[]; mat_beta_move=[]; mat_gamma_move=[];
for i=1:size(BG_2plot,1)
    move_ind = find(timeAxis >= (T.Move_ms(i)-T.Triger_ms(i)-win)/1000, 1, 'first'); % To remove the win (--win-win=0)
    BG_mc_save = BG_2save(i,move_ind-100:move_ind+10*20-1);
    BG_2save_move(i,1:length(BG_mc_save)) = BG_mc_save;
    BG_mc_plot = BG_2plot(i,move_ind-100:move_ind+10*20-1);   
    BG_2plot_move(i,1:length(BG_mc_plot)) = BG_mc_plot;
    mat_delta_mc_move = mat_delta(i,move_ind-100:move_ind+10*20-1);
    mat_delta_move(i,1:length(mat_delta_mc_move)) = mat_delta_mc_move;
    mat_beta_mc_move = mat_beta(i,move_ind-100:move_ind+10*20-1);
    mat_beta_move(i,1:length(mat_beta_mc_move)) = mat_beta_mc_move;
    mat_gamma_mc_move = mat_gamma(i,move_ind-100:move_ind+10*20-1);
    mat_gamma_move(i,1:length(mat_gamma_mc_move)) = mat_gamma_mc_move;
end
timeAxis_move = (((0:size(BG_2save_move,2)-1) * tdiff)-win*2)/1000; %in s

figure;
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
imagesc(timeAxis_move,1:size(BG_2save_move),BG_2save_move)
axis xy % important: row 1 at top
xlabel("Time from Meadow Movement (s)")
hold on
title(animal);
% Plot start time horizontal line:
xline(0,'g' ,'LineWidth', 1.5, 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
% xticks(-10:1:20);

thr = prctile(BG_2save_move(:),95); %95
clim([min(BG_2save_move(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'B/G ratio'; % Set label text
% c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis_move, nanmean(BG_2save_move,1),'LineWidth', 1.5)
xlabel('Time from Meadow Movement (s)'); ylabel('avg. B2G')
% ylim([10 85])
xline(0,'g' ,'LineWidth', 1.5)
% xticks(-10:1:20);
linkaxes([ax1,ax2],'x');




figure;
%add subplot for the B/G ratio
ax1 = subplot('Position', [0.1, 0.3, 0.7, 0.6]);
imagesc(timeAxis_move,1:size(BG_2plot_move),BG_2plot_move)
axis xy % important: row 1 at top
xlabel("Time from Meadow Movement (s)")
hold on
title(animal);
% Plot start time horizontal line:
xline(0,'g' ,'LineWidth', 1.5, 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
ylabel('Trial #')
% xticks(-10:1:20);

thr = prctile(BG_2plot_move(:),85); %95
clim([min(BG_2plot_move(:)),thr]);

%add colorbar
c = colorbar;
c.Position = [0.82 0.3 0.02 0.6];
% Set vertical label for colorbar
c.Label.String = 'D/B ratio'; % Set label text
% c.Label.Position = [4, 80, 0]; % Adjust label position (relative to colorbar)
c.Label.Rotation = 90; % Set label rotation to 0 degrees for vertical alignment

% Plot Average
ax2 = subplot('Position', [0.1, 0.1, 0.7, 0.15]);
plot(timeAxis_move, nanmean(BG_2plot_move,1),'LineWidth', 1.5)
xlabel('Time from Meadow Movement (s)'); ylabel('avg. D2B')
% ylim([10 85])
xline(0,'g' ,'LineWidth', 1.5)
% xticks(-10:1:20);
linkaxes([ax1,ax2],'x');

mean_band_d_move = nanmean(mat_delta_move);
mean_band_b_move = nanmean(mat_beta_move);
mean_band_g_move = nanmean(mat_gamma_move);



% fs1 = 10;               % samples per second
% n = 240;  %120 sound old %460 screen diode   240 sound new   600 hunter      
% time1 = (0:n-1)/fs1;              
% time1_centered = time1 - 12;     % -5 hunter




figure;
subplot(3,1,1);
plot(timeAxis_move,mean_band_d_move,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Sound Onset (s)"); 
ylabel('Power');
% ylim([0 500]);
set(gca,'FontSize',18)
title('PV228 white: Delta = 0-5 Hz','FontSize',18);

subplot(3,1,2);
plot(timeAxis_move,mean_band_b_move,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
% xlabel("Time from Meadow Movement (s)");
ylabel('Power');
title('Beta = 10-30 Hz','FontSize',18);
set(gca,'FontSize',18)

subplot(3,1,3);
plot(timeAxis_move,mean_band_g_move,'k','LineWidth',2); 
% xline(0,'r','LineWidth',2 , 'Label', 'Sound onset', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'g','LineWidth',2 , 'Label', 'Bug On', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
xline(0,'g','LineWidth',2 , 'Label', 'Movement', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left');
% xline(0,'r','LineWidth',2 , 'Label', 'Static', 'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'Left')
xlabel("Time from Meadow Movement (s)");  ylabel('Power'); %Time from Meadow Movement   Bug Appearance
title('Gamma = 60-80 Hz','FontSize',18);
set(gca,'FontSize',18)

%% 3. Getting the matrices


dataName = sprintf('reg_meadow_%s', animal);

data = struct();

data.(sprintf('BG_reg_meadow_%s', animal)) = BG_2save;
data.(sprintf('DB_reg_meadow_%s', animal)) = BG_2plot;
data.(sprintf('delta_reg_meadow_%s', animal)) = mat_delta;
data.(sprintf('beta_reg_meadow_%s', animal)) = mat_beta;
data.(sprintf('gamma_reg_meadow_%s', animal)) = mat_gamma;

data.(sprintf('BG_reg_meadow_%s_move', animal)) = BG_2save_move;
data.(sprintf('DB_reg_meadow_%s_move', animal)) = BG_2plot_move;
data.(sprintf('delta_reg_meadow_%s_move', animal)) = mat_delta_move;
data.(sprintf('beta_reg_meadow_%s_move', animal)) = mat_beta_move;
data.(sprintf('gamma_reg_meadow_%s_move', animal)) = mat_gamma_move;




eval(sprintf('%s = data;', dataName));

save('mats4stat.mat', dataName, '-append');


%%
%% Spectrogram of a vertical step, using the same FFT settings as getBeta2GammaRatioOfir
clear; clc;

%% --- Parameters (copied from your script) ---
Fs            = 1000;      % Hz  (stands in for obj.filt.FFs - change to your real rate)
movWin        = 2000;      % ms  moving window
movOLWin      = 1900;      % ms  window overlap
segmentWelch  = 1000;      % ms  Welch segment
dftPointsWelch= 2^10;
OLWelch       = 0.5;
gammaLow = 60; gammaHigh = 80;
betaLow  = 10; betaHigh  = 30;

%% --- Step signal (vertical edge: 0 -> 1 at t = tStep) ---
dur   = 10;                % s
tStep = 5;                 % s
t = (0:1/Fs:dur-1/Fs);
x = double(t >= tStep);    % sampled step: edge is a single-sample jump

%% --- Derived values (same formulas as your script) ---
movWinSamples        = movWin/1000*Fs;
movOLWinSamples      = movOLWin/1000*Fs;
timeBin              = movWin - movOLWin;                 % ms
segmentWelchSamples  = round(segmentWelch/1000*Fs);
samplesOLWelch       = round(segmentWelchSamples*OLWelch);

%% --- Frequency axis and band indices ---
[~,f] = pwelch(randn(1,movWinSamples),segmentWelchSamples,samplesOLWelch,dftPointsWelch,Fs);
pGamma = find(f>=gammaLow & f<gammaHigh);
pBeta  = find(f>=betaLow  & f<betaHigh);

%% --- Moving-window Welch ---
xB = buffer(x, movWinSamples, movOLWinSamples, 'nodelay');   % one window per column
pValid = all(~isnan(xB));
[pxx,f] = pwelch(xB(:,pValid), segmentWelchSamples, samplesOLWelch, dftPointsWelch, Fs);

nWin  = size(xB,2);
t_ms  = (movWin/2) + (0:nWin-1)*timeBin;                     % window centers
t_ms  = t_ms(pValid);

beta2gamma = (mean(pxx(pBeta,:)) ./ mean(pxx(pGamma,:)))';
betaBand   = mean(pxx(pBeta,:))';
gammaBand  = mean(pxx(pGamma,:))';

%% --- Plots ---
figure('Color','w','Position',[100 100 900 750]);

subplot(3,1,1);
plot(t, x, 'LineWidth', 2); grid on;
ylim([-0.2 1.2]); xlabel('Time (s)'); ylabel('Amplitude');
title('Step signal');

subplot(3,1,2);
imagesc(t_ms/1000, f, 10*log10(pxx + eps)); axis xy;
ylim([0 120]); colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title(sprintf('Spectrogram (win %d ms, overlap %d ms, Welch seg %d ms, nfft %d)', ...
    movWin, movOLWin, segmentWelch, dftPointsWelch));
yline(betaLow,'w--'); yline(betaHigh,'w--');
yline(gammaLow,'k--'); yline(gammaHigh,'k--');

subplot(3,1,3);
plot(t_ms/1000, beta2gamma, 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Beta / Gamma');
title('Beta (10-30 Hz) to Gamma (60-80 Hz) ratio');