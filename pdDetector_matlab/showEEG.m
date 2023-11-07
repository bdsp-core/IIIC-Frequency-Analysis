  function f = showEEG(seg_bi) 
  
  f =gcf;
%   Ax_EEG = subplot('position',[.05 .05 .93 .9]);
  Ax_EEG = gca; 
  
  gap = nan(1,size(seg_bi,2));
  
  % Parameter settings %
  zScale = 1/100;
  Fs = 200;
  [B, A] = butter(3, [5 40]/(.5*Fs));
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
  
  
  
  set(f,'CurrentAxes',Ax_EEG);cla(Ax_EEG)
  
  hold(Ax_EEG, 'on')
  %                 title(strrep(strrep(file, '_', ' '), '.mat', ''))
  
  % spikes %
  %plot(t(loc_==1),  nCh+.2, 'rs', 'markerfacecolor', 'r')
  %plot(t(loc_1==1), nCh+.4, 'gs', 'markerfacecolor', 'g')
  %%
  plot(Ax_EEG, t, zScale*data+DCoff,'k');  box on;
  
  set(Ax_EEG, 'ytick',1:nCh,'yticklabel',flipud(channel_withspace),'box','on', 'ylim', [0 nCh+1], 'xlim', [to t1], 'xtick',round(t(1):2*Fs:t(end)),'xticklabel',timeStamps)
  xlabel(Ax_EEG, 'Time')
  
  %% 1-sec marks.
  c_gray = [0.7451, 0.7451, 0.7451];
  theseTicks = [t(1):1*Fs:t(end)];
  vline(theseTicks);
    
  %% Scale rulers % for sanity check %
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
    end