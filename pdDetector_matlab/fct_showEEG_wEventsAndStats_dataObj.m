function [ f1, thisStatsObj] = fct_showEEG_wEventsAndStats_dataObj(data_obj)
%fct_showEEG_wEvents Re-write of showEEG() to do overlay
Fs = 200;
f1 = figure('Position',[965 119 1023 1125],'Visible','off');
numPlots = 10;
% Load from data_obj
seg_bi = data_obj.seg_bi;
bipolarEvents = data_obj.pdEvents_bi_select;
thisLOC = data_obj.thisLOC ;
stats_obj = data_obj.stats_obj;

%%
Ax_events = subplot(numPlots,1,1);
%Ax_EEG = subplot(3,1,2,'position',[.05 .05 .93 .9]);
% imagesc(thisLOC);    % stopped here!! need to merge LOCs (or get actual LOC from the classification)
[pks,locs,w,p]= findpeaks(smooth(thisLOC,10),'MinPeakProminence',2,'MinPeakDistance',20);

plot(smooth(thisLOC,10));
text(locs,pks,'*','FontSize',14);
axis tight;
ylabel('Discharge groups');

Ax_events.XAxis.Visible = 'off';

%%  evolution in frequency
% locs_aug = [0 locs' length(thisLOC)];
% IDI_pks = diff(locs_aug)/Fs;

% figure; plot(indx_pks, IDI_pks.^-1);
if length(locs)>2

IDI_pks = diff(locs)/Fs;
indx_pks = 1:length(IDI_pks);
    
mdl=fitlm(indx_pks,IDI_pks.^-1,'y~x1');     % does not have good intuition for double peaks, weak peaks, and onset/offset
    evolution_intercept = mdl.Coefficients.Estimate(1);
    evolution_slope = mdl.Coefficients.Estimate(2);   % units is delta(Hz)/delta(step)
    evolution_SE = mdl.Coefficients.SE(2);  % standard error
    evolution_RMSE = mdl.RMSE;
    evolution_pVal = mdl.Coefficients.pValue(2);
    evolution_range = [evolution_intercept, evolution_intercept+evolution_slope*indx_pks(end)];
else

    evolution_slope = nan;
    evolution_SE = nan;
    evolution_RMSE = nan;
    evolution_pVal = nan;
    evolution_range = [nan, nan];
end

%% spatial extent evolution 
% figure; plot(indx_pks_amp, pks);
if length(pks)>1
indx_pks_amp = 1:length(pks);

    mdl_amp=fitlm(indx_pks_amp ,pks,'y~x1');     % does not have good intuition for double peaks, weak peaks, and onset/offset
    evolution_amp_intercept = mdl_amp.Coefficients.Estimate(1);
    evolution_amp_slope = mdl_amp.Coefficients.Estimate(2);   % units is delta(channels)/delta(step)
    evolution_amp_SE = mdl_amp.Coefficients.SE(2);  % standard error
    evolution_amp_RMSE = mdl_amp.RMSE;
    evolution_amp_pVal = mdl_amp.Coefficients.pValue(2);
    evolution_amp_range = [evolution_amp_intercept, evolution_amp_intercept+evolution_amp_slope*indx_pks_amp(end)];

else
    evolution_amp_intercept = nan;
    evolution_amp_slope = nan;
    evolution_amp_SE = nan;
    evolution_amp_RMSE = nan;
    evolution_amp_pVal = nan;
    evolution_amp_range = [nan, nan];

end 

%% 

%Ax_EEG = gca; 
Ax_EEG = subplot(numPlots,1,[2,3,4]); 

seg_bi= eegfilt(seg_bi,200,1,40,0,[],0,'fir1',[]);

gap = nan(1,size(seg_bi,2));

% Parameter settings %
zScale = 1/100;

% [B, A] = butter(3, [5 40]/(.5*Fs));
channel_withspace = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; '';'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; '';'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; '';'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; '';'Fz-Cz';'Cz-Pz'};
ii = 1;

set(Ax_EEG, 'Visible', 'on');

% dispay %
data = [seg_bi(  1: 4,:); gap; seg_bi(  5: 8,:); gap; seg_bi(  9:12,:); gap; seg_bi( 13:16,:); gap; seg_bi( 17:18,:)];

nCh = size(data, 1);
DCoff = repmat(flipud((1:nCh)'), 1, size(seg_bi , 2));

To = datetime('00:00:00');
to = 1;
t1 = size(seg_bi, 2);

t = to:t1;
timeStamps = datestr(To + seconds(round(to/Fs:2:t1/Fs)), 'hh:MM:ss');

set(f1,'CurrentAxes',Ax_EEG);cla(Ax_EEG)

hold(Ax_EEG, 'on')
%                 title(strrep(strrep(file, '_', ' '), '.mat', ''))

%% spikes 
%plot(t(loc_==1),  nCh+.2, 'rs', 'markerfacecolor', 'r')
%plot(t(loc_1==1), nCh+.4, 'gs', 'markerfacecolor', 'g')

plot(Ax_EEG, t, zScale*data+DCoff,'k');  box on;

set(Ax_EEG, 'ytick',1:nCh,'yticklabel',flipud(channel_withspace),'box','on', 'ylim', [0 nCh+1], 'xlim', [to t1], 'xtick',round(t(1):2*Fs:t(end)),'xticklabel',timeStamps)
xlabel(Ax_EEG, 'Time')

hold on;
%%
channelLabels_bipolar_withspace_ud= flipud({'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; '';'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; '';'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; '';'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; '';'Fz-Cz';'Cz-Pz'});

red_c = [1 0 0 0.5];
blue_c = [0 0 1 0.5];
green_c = [0 1 0 0.5];
cP = {red_c, blue_c, green_c};

channelColors = zeros(22,1);
channelColors([1:4,11:14]) = 1;
channelColors([6:9,16:19]) = 2;
channelColors([21,22]) = 3;
channelColors_ud = flipud(channelColors);

theseRect = cell(size(bipolarEvents,1),1);
for j=1:size(bipolarEvents,1)
    thisEvent = bipolarEvents(j,:);
    thisChanName = thisEvent.channel;
    chan= thisEvent.channelNum;
    thisBB = thisEvent.BoundingBox;
    indx_chan = find(ismember(channelLabels_bipolar_withspace_ud,thisChanName));
    thisBB(2) = indx_chan-0.5;
    
    thisRect = rectangle('Position',thisBB, 'FaceColor', cP{channelColors_ud(indx_chan)}, 'EdgeColor','none');
    
    theseRect{j} = thisRect;
end
  %% 1-sec marks.
  c_gray = [0.7451, 0.7451, 0.7451];
  theseTicks = [t(1):1*Fs:t(end)];
  vline(theseTicks);
    
%%
% Scale rulers % for sanity check %
dt = t1-to+1;
a = round(dt*4/5);

xa1 = to+[a a+Fs-1];
ya1 = [9 9];


xa2 = to+[a a];
ya2 = ya1+[0 100*zScale];

text(xa1(1)-.45*a/10,    mean(ya2), '100\muV','Color', 'b','FontSize',14);
text(mean(xa1)-0.2*a/10, 8.7, '1 sec','Color', 'b','FontSize',14);
line(xa1,ya1, 'LineWidth', 2, 'Color','b');
line(xa2,ya2, 'LineWidth', 2, 'Color','b');

hold(Ax_EEG, 'off')


%% Stats

% Row 1
Ax_stats = subplot(numPlots,1,5);
set(Ax_stats,'YTickLabel',{},'XTickLabel',{},'XTick',[],'YTick',[]);
Ax_stats.XAxis.Visible = 'off';
Ax_stats.YAxis.Visible = 'off';
% 
% text(0,0,'(0,0)','Units','normalized','FontSize',10);
% text(1,0,'(1,0)','Units','normalized','FontSize',10);
% text(0,1,'(0,1)','Units','normalized','FontSize',10);
% text(1,1,'(1,1)','Units','normalized','FontSize',10);
%  
pow_L = stats_obj.pow_L;
pow_R = stats_obj.pow_R;
pow_G = stats_obj.pow_G;
relPLvPG = stats_obj.relPLvPG;
relPRvPG = stats_obj.relPRvPG;
relPLvPR = log(pow_L/pow_R);
stats_obj.relPLvPR = log(pow_L/pow_R);   % overwrites previous calculation

stats_obj_t = struct2table(stats_obj);
% Get the table in string form.
TString = evalc('disp(stats_obj_t)');
% Use TeX Markup for bold formatting and underscores.
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
% Get a fixed-width font.
FixedWidth = get(0,'FixedWidthFontName');
% Output the table using the annotation command.
% annotation(Ax_stats,'Textbox','String',TString,'Interpreter','Tex',...
%     'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

Ax_stats.FontName = FixedWidth;
text(0.5,0.5,TString,'Units','normalized','HorizontalAlignment','center','FontName',FixedWidth,'FontSize',8)

% Row 2
Ax_stats2 = subplot(numPlots,1,6);
set(Ax_stats2,'YTickLabel',{},'XTickLabel',{},'XTick',[],'YTick',[]);
Ax_stats2.XAxis.Visible = 'off';
Ax_stats2.YAxis.Visible = 'off';

spatialExtent = data_obj.spatialExtent; 
freq = data_obj.freq;
channels = data_obj.channels;
interpretation = data_obj.interpretation;

stats_obj2 = struct;
stats_obj2.freq = freq;
% stats_obj2.evolution_slope = data_obj.thisEvolution_slope;
% stats_obj2.evolution_RMSE = data_obj.thisEvolution_RMSE;
%stats_obj2.evolution_intercept = evolution_intercept;
stats_obj2.evolution_slope = evolution_slope;  % take from this function, not data_obj
%stats_obj2.evolution_RMSE = evolution_RMSE; % take from this function, not data_obj
stats_obj2.evolution_pVal = evolution_pVal;
stats_obj2.evolution_range = evolution_range;
    
stats_obj2_t = struct2table(stats_obj2);

TString2 = evalc('disp(stats_obj2_t)');
% Use TeX Markup for bold formatting and underscores.
TString2 = strrep(TString2,'<strong>','\bf');
TString2 = strrep(TString2,'</strong>','\rm');
TString2 = strrep(TString2,'_','\_');

text(0.5,0.5,TString2,'Units','normalized','HorizontalAlignment','center','FontName',FixedWidth,'FontSize',8)

% Row 2B
Ax_stats2 = subplot(numPlots,1,7);
set(Ax_stats2,'YTickLabel',{},'XTickLabel',{},'XTick',[],'YTick',[]);
Ax_stats2.XAxis.Visible = 'off';
Ax_stats2.YAxis.Visible = 'off';

spatialExtent = data_obj.spatialExtent; 
interpretation = data_obj.interpretation;

stats_obj3 = struct;
stats_obj3.spatialExtent = spatialExtent; 
stats_obj3.interpretation = interpretation;

stats_obj3.evolution_amp_extent = evolution_amp_intercept/18;  % take from this function, not data_obj
stats_obj3.evolution_amp_slope = evolution_amp_slope;  % take from this function, not data_obj
stats_obj3.evolution_amp_pVal = evolution_amp_pVal;
stats_obj3.evolution_amp_range= evolution_amp_range;

stats_obj3_t = struct2table(stats_obj3);

TString2 = evalc('disp(stats_obj3_t)');
% Use TeX Markup for bold formatting and underscores.
TString2 = strrep(TString2,'<strong>','\bf');
TString2 = strrep(TString2,'</strong>','\rm');
TString2 = strrep(TString2,'_','\_');

text(0.5,0.5,TString2,'Units','normalized','HorizontalAlignment','center','FontName',FixedWidth,'FontSize',8)

% Row 3
Ax_stats3 = subplot(numPlots,1,8);
set(Ax_stats3,'YTickLabel',{},'XTickLabel',{},'XTick',[],'YTick',[]);
Ax_stats3.XAxis.Visible = 'off';
Ax_stats3.YAxis.Visible = 'off';
% print Y_human
Y_human = data_obj.Y_human;

TString2 = evalc('disp(Y_human)');
% Use TeX Markup for bold formatting and underscores.
TString2 = strrep(TString2,'<strong>','\bf');
TString2 = strrep(TString2,'</strong>','\rm');
TString2 = strrep(TString2,'_','\_');

text(0.5,0.5,TString2,'Units','normalized','HorizontalAlignment','center','FontName',FixedWidth,'FontSize',8)

% Row 4
% Ax_stats4 = subplot(numPlots,1,9);
% set(Ax_stats4,'YTickLabel',{},'XTickLabel',{},'XTick',[],'YTick',[]);
% Ax_stats4.XAxis.Visible = 'off';
% Ax_stats4.YAxis.Visible = 'off';
% print Y_human
P_model = data_obj.P_model;

P_model_vars = P_model.Properties.VariableNames;
P_model_rowVars = {'0-2s','2-4s','4-6s','6-8s','8-10s','10-12s','12-14s'};

P_model_t=cell2table(table2cell(P_model)','RowNames',P_model_vars,'VariableNames',P_model_rowVars);

TString3 = evalc('disp(P_model_t)');
% Use TeX Markup for bold formatting and underscores.
TString3 = strrep(TString3,'<strong>','\bf');
TString3 = strrep(TString3,'</strong>','\rm');
TString3 = strrep(TString3,'_','\_');

% text(0.5,0.5,TString3,'Units','normalized','HorizontalAlignment','center','FontName',FixedWidth,'FontSize',8)

Ax_stats4 = subplot(numPlots,1,9:10);
imagesc(cell2mat(table2cell(P_model_t)));
set(gca,'YMinorTick','on');
set(gca,'XTickLabel',P_model_rowVars);
set(gca,'YTickLabel',P_model_vars);
% set(gca,'YLabel','P_model');
colorbar;


%% 
suptitle([num2str(data_obj.num) '-' strrep(data_obj.filename,'_','\_')]);
set(f1, 'PaperPositionMode','auto', 'Visible','off');

thisStatsObj = cat(2, stats_obj_t, stats_obj2_t, stats_obj3_t);

end

