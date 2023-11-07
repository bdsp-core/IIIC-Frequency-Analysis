function freqSpike_GUI_v2
    
    % more controls on parameters %

    % Initialize - figure settings % 
    close all
    clc
    warning('off','all')
    
      
    % Register %
    %user_ID = inputdlg('Please enter your initials','user_ID'); % ask for user id %
      
    
    f = figure('units','normalized','outerposition',[0 0 1 1]);
    set(f,'HitTest','off')
    set(f,'WindowButtonDownFcn',@clicks_Callback);
    set(f,'ButtonDownFcn',@clicks_Callback);
    set(f,'KeyPressFcn',@keys_Callback);
    
    % Dirs %
    dataDir = [pwd, '\Data\GPD\'];
 
    files = struct2cell(dir([dataDir, '*.mat']))';
    files  = files(:, 1);

    % Buttons %
    startx = uicontrol(f,'style','pushbutton','units','normalized',  'position',[0.4690    0.58    0.0520    0.0250], 'string','Start','callback',@fcn_start);      
    Ax_EEG = subplot('position',[.05 .05 .93 .9]);
    
    winSize    = 5;              
    winSizeH   =  uicontrol('Style', 'edit', 'String',  '5','units','normalized','Position', [.05 .97  0.02 0.0250],'Callback', @fcn_WINDOWSIZE, 'HorizontalAlignment','right');
    winSizeStr =  uicontrol('Style', 'text', 'String',  's','units','normalized','Position', [.07 .97  0.01 0.0250],'HorizontalAlignment','left');
    
    stepSize = 20;              
    stepSizeH   =  uicontrol('Style', 'edit', 'String',  '20','units','normalized','Position', [.08 .97  0.02 0.0250],'Callback', @fcn_STEPSIZE, 'HorizontalAlignment','right');
    stepSizeStr =  uicontrol('Style', 'text', 'String',   '%','units','normalized','Position', [.10 .97  0.01 0.0250],'HorizontalAlignment','left');
    
    smoothWinSize = 120;       
    smoothWinSizeH   =  uicontrol('Style', 'edit', 'String',  '120','units','normalized','Position', [.11 .97  0.02 0.0250],'Callback', @fcn_SMWINSIZE, 'HorizontalAlignment','right');
    smoothWinSizeStr =  uicontrol('Style', 'text', 'String',   'ms','units','normalized','Position', [.13 .97  0.01 0.0250],'HorizontalAlignment','left');
    
    thrCh = 9;              
    thrChH   =  uicontrol('Style', 'edit', 'String',   '9','units','normalized','Position', [.14 .97  0.02 0.0250],'Callback', @fcn_THRCH, 'HorizontalAlignment','right');
    thrChStr =  uicontrol('Style', 'text', 'String',  'ch','units','normalized','Position', [.16 .97  0.01 0.0250],'HorizontalAlignment','left');
    
    thrDur = 40;              
    thrDurH   =  uicontrol('Style', 'edit', 'String',   '40','units','normalized','Position', [.17 .97  0.02 0.0250],'Callback', @fcn_THRDUR, 'HorizontalAlignment','right');
    thrDurStr =  uicontrol('Style', 'text', 'String',  'ms','units','normalized','Position', [.19 .97  0.01 0.0250],'HorizontalAlignment','left');
    
    eventFreq = NaN;              
    eventFreqH   =  uicontrol('Style', 'text', 'String',   'Hz','units','normalized','Position', [.22 .97  0.05 0.0250], 'HorizontalAlignment','right');

    greyBar = 0;
    greyBarH   =  uicontrol('Style', 'checkbox', 'String',   'Highlight', 'Value',0, 'units','normalized','Position', [.28 .97  0.05 0.0250], 'HorizontalAlignment','right','Callback', @fcn_GreyBar);

    % Parameter settings %
    zScale = 1/100;
    Fs = 200;
    [B, A] = butter(3, [5 40]/(.5*Fs));
    channel_withspace = {'Fp1-F7';'F7-T3';'T3-T5';'T5-O1'; '';'Fp2-F8';'F8-T4';'T4-T6';'T6-O2'; '';'Fp1-F3';'F3-C3';'C3-P3';'P3-O1'; '';'Fp2-F4';'F4-C4';'C4-P4';'P4-O2'; '';'Fz-Cz';'Cz-Pz'};
    ii = 1;

    set(Ax_EEG, 'Visible', 'off');
    set(winSizeH, 'Visible', 'off');set(winSizeStr, 'Visible', 'off');
    set(stepSizeH, 'Visible', 'off');set(stepSizeStr, 'Visible', 'off');
    set(smoothWinSizeH, 'Visible', 'off');set(smoothWinSizeStr, 'Visible', 'off');
    set(thrChH, 'Visible', 'off'); set(thrChStr, 'Visible', 'off');
    set(thrDurH, 'Visible', 'off');set(thrDurStr, 'Visible', 'off');
    set(eventFreqH, 'Visible', 'off');
    set(greyBarH, 'Visible', 'off');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uiwait(f);
    
    %try 
        while true  
            set(winSizeH, 'Enable', 'on');set(winSizeStr, 'Enable', 'on');
            set(stepSizeH, 'Enable', 'on');set(stepSizeStr, 'Enable', 'on');
            set(smoothWinSizeH, 'Enable', 'on');set(smoothWinSizeStr, 'Enable', 'on');
            set(thrChH, 'Enable', 'on'); set(thrChStr, 'Enable', 'on');
            set(thrDurH, 'Enable', 'on');set(thrDurStr, 'Enable', 'on');
            set(greyBarH, 'Enable', 'on');
            drawnow
            
            file = LUT{ii};
            tmp = load([dataDir, file]);
            
%             % band pass filter %
             seg = tmp.seg;
%             for ch = 1:size(seg, 1)
%                 seg(ch,:) = filter(B, A, seg(ch,:));
%             end
            
            
            gap = NaN(1, size(seg , 2));
            %seg_car = seg - repmat(mean(seg, 1), 19, 1); 
            seg_bi = fcn_LBipolar(seg);
            

            % gpd detection - train of 1s and 0s %
            LOC = fct_eventDetector(seg_bi, Fs, winSize*Fs, (stepSize/100*winSize)*Fs, round((smoothWinSize/1000)*Fs));
            loc_ = sum(LOC, 1);
            tmp_ = loc_;
            loc_(tmp_>=thrCh) = 1; loc_(tmp_ <thrCh) = 0;
            [loc_1, loc_2, counts, iOnsets, iOffsets] = fct_labelSmoother(loc_, 2, (thrDur/1000*Fs));
            
            if counts == 0
                eventFreq = NaN;
            else
                eventFreq = counts/(size(seg, 2)/Fs);
            end
            set(eventFreqH, 'string', ['Freq.: ', num2str(round(100*eventFreq)/100), ' Hz'])
            
            % dispay %
            data = [seg_bi(  1: 4,:); gap; seg_bi(  5: 8,:); gap; seg_bi(  9:12,:); gap; seg_bi( 13:16,:); gap; seg_bi( 17:18,:)];

            nCh = size(data, 1);
            DCoff = repmat(flipud((1:nCh)'), 1, size(seg , 2));

            To = datetime('00:00:00');
            to = 1;
            t1 = size(seg, 2);

            t = to:t1;
            timeStamps = datestr(To + seconds(round(to/Fs:2:t1/Fs)), 'hh:MM:ss');



            set(f,'CurrentAxes',Ax_EEG);cla(Ax_EEG)
    
            hold(Ax_EEG, 'on')   
                title(strrep(strrep(file, '_', ' '), '.mat', ''))

                % spikes %
                %plot(t(loc_==1),  nCh+.2, 'rs', 'markerfacecolor', 'r')
                %plot(t(loc_1==1), nCh+.4, 'gs', 'markerfacecolor', 'g')
                
                if greyBar
                    for i_ =1:length(iOnsets)
                        
                        fill([iOnsets(i_) iOffsets(i_) iOffsets(i_) iOnsets(i_)], [0 0 nCh+1 nCh+1], [.95 .95 .95], 'edgecolor', [.95 .95 .95])
                    end
                else
                    plot(t(loc_2==1), nCh+.6, 'bs', 'markerfacecolor', 'b')
                    
                end


                % line on every sec %
                for iSec = 1:round((t1-to+1)/Fs)
                    ta = to + Fs*(iSec-1);
                    line([ta ta],  [0 nCh+1], 'linestyle', '--', 'color', [.5 .5 .5])
                end

                plot(Ax_EEG, t, zScale*data+DCoff,'k');  box on;

                set(Ax_EEG, 'ytick',1:nCh,'yticklabel',flipud(channel_withspace),'box','on', 'ylim', [0 nCh+1], 'xlim', [to t1], 'xtick',round(t(1):2*Fs:t(end)),'xticklabel',timeStamps)
                xlabel(Ax_EEG, 'Time')


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
             uiwait(f);


        end
        
%     catch err
%         choice = questdlg('Exit?', ...
%                             'Warning', ...
%                             'No, I am good.','Yes, please!','No, I am good.');
% 
%         switch choice
%             case 'Yes, please!'
%                 try
%                     delete(f);
%                 catch err
%                 end
% 
%             otherwise
%                uiwait(f);
%         end
%     end
%     
    
    %%%% Callbacks %%%
    % Start %
    function fcn_start(varargin) 
        
        set(startx,'Visible','off','Enable', 'off');
        loading = uicontrol(f,'style','text','units','normalized','position',[0.4690    0.5000+.08    0.0520    0.0250],'string','Loading...','FontSize',11);
        drawnow
        
        LUT = files;
 
      
        delete(loading)

        % Reveal axis buttons etc. %
        set(Ax_EEG, 'Visible', 'on');    
           
        set(winSizeH, 'Visible', 'on');set(winSizeStr, 'Visible', 'on');
        set(stepSizeH, 'Visible', 'on');set(stepSizeStr, 'Visible', 'on');
        set(smoothWinSizeH, 'Visible', 'on');set(smoothWinSizeStr, 'Visible', 'on');
        set(thrChH, 'Visible', 'on'); set(thrChStr, 'Visible', 'on');
        set(thrDurH, 'Visible', 'on');set(thrDurStr, 'Visible', 'on');
        set(eventFreqH, 'Visible', 'on');
        set(greyBarH, 'Visible', 'on');
        
        uiresume(f);
    end

    % Keys %
     %%%% Keyboards %%%%
    function keys_Callback(f,varargin)

        k = get(f,'CurrentKey');
        
        switch k
            case 'rightarrow'
                ii = min(ii+1, size(LUT, 1));
                
            case 'leftarrow'
                ii = max(ii-1, 1);
                
            case 'uparrow'
                zScale =  zScale*1.5;
                
            case 'downarrow'
                zScale =  zScale/1.5;
                
            case {'1', 'numpad1', 'end'}
                score{ii, end} = 1;
                title(['No.', num2str(ii), ' [\color{red}TPW', '\color{black}] ', num2str(m), ' left'])
                fcn_next();
                
            case {'0', 'numpad0', 'insert'}
                score{ii, end} = 0;
                title(['No.', num2str(ii), ' [\color{blue}BG', '\color{black}] ', num2str(m), ' left'])
                fcn_next();
                
            otherwise
                disp(k)
       
        end
        uiresume(f);
    end

    function dataBipolar = fcn_LBipolar(data)
    %labels = {'Fp1';'F3';'C3';'P3';'F7';'T3';'T5';'O1';'Fz';'Cz';'Pz';'Fp2';'F4';'C4';'P4';'F8';'T4';'T6';'O2'};

        dataBipolar(9,:) = data(1,:) - data(2,:); % Fp1-F3
        dataBipolar(10,:) = data(2,:) - data(3,:); % F3-C3
        dataBipolar(11,:) = data(3,:) - data(4,:); % C3-P3
        dataBipolar(12,:) = data(4,:) - data(8,:); % P3-O1

        dataBipolar(13,:) = data(12,:) - data(13,:); % Fp2-F4
        dataBipolar(14,:) = data(13,:) - data(14,:); % F4-C4
        dataBipolar(15,:) = data(14,:) - data(15,:); % C4-P4
        dataBipolar(16,:) = data(15,:) - data(19,:); % P4-O2

        dataBipolar(1,:) = data(1,:) - data(5,:);  % Fp1-F7
        dataBipolar(2,:) = data(5,:) - data(6,:); % F7-T3
        dataBipolar(3,:) = data(6,:) - data(7,:); % T3-T5
        dataBipolar(4,:) = data(7,:) - data(8,:); % T5-O1

        dataBipolar(5,:) = data(12,:) - data(16,:); % Fp2-F8
        dataBipolar(6,:) = data(16,:) - data(17,:); % F8-T4
        dataBipolar(7,:) = data(17,:) - data(18,:); % T4-T6
        dataBipolar(8,:) = data(18,:) - data(19,:); % T6-O2

        dataBipolar(17,:) = data(9,:) - data(10,:);   % Fz-Cz
        dataBipolar(18,:) = data(10,:) - data(11,:); % Cz-Pz
    end

    function LOC = fct_eventDetector(eeg_bi, Fs, w, s, sw)

        %N = floor(size(eeg_bi, 2)/Fs)-4; % 
        %w = 5*Fs;    % window size %
        %s = 1*Fs;    % stepsize %
        %sw = .12*Fs; % smoothing window %
        LOC = NaN(size(eeg_bi));
        N = floor((size(eeg_bi,2) - (w-s))/s);

        for i = 1:size(eeg_bi, 1)
            x = eeg_bi(i,:);

            % nleo %
            z = fct_nleo(x);

            % smooth %
            z = smooth(z,  sw);

            for n = 1:N
                aa = (n-1)*s+1;
                bb = aa+w-1;

                seg_ = z(aa:bb);
                thr_ = .6*(std(seg_)+prctile(seg_, 75));

                LOC(i, aa:bb) = (seg_>=thr_);

            end
        end

        % callbacks %
        function y = fct_nleo(x)

            x = [0 0 0 x];
            y = NaN(1, length(x));
            for j = 4:length(x)
                y(j) = x(j-1)*x(j-2) - x(j)*x(j-3);
            end

            y = y(4:end);

        end
    end

    function [z, h, counts, iOnsets, iOffsets] = fct_labelSmoother(y, thr_sm, thr_d)
        
        % smooth %
        z = y;
        a_ = NaN; b_ = NaN;
        c = 1;
        while c <=length(y)
            z_ = z(c);

            if z_== 0 
                if isnan(a_) % the start
                    a_ = c;
                    b_ = c;
                else
                    b_ = c;
                end

            else % 1
                if ~isnan(a_)
                    d = b_-a_+1;

                    if d <= thr_sm
                        z(a_:b_) = 1;
                    end

                    a_ = NaN;
                    b_ = NaN;
                end
            end
            c = c+1;
        end

        h = z;
        a_ = NaN; b_ = NaN;
        c = 1;
        counts = 0;
        iOnsets = [];
        iOffsets = [];
        while c <=length(y)
            h_ = h(c);

            if h_== 1 
                if isnan(a_) % the start
                    a_ = c;
                    b_ = c;
                else
                    b_ = c;
                end

            else % 0
                if ~isnan(a_)
                    d = b_-a_+1;

                    if d<=thr_d
                        h(a_:b_) = 0;
                    else
                        counts = counts+1;
                        iOnsets = [iOnsets; a_];
                        iOffsets = [iOffsets; b_];
                    end

                    a_ = NaN;
                    b_ = NaN;
                end
            end
            c = c+1;
        end       
    end


    % interval edit box %
    function fcn_WINDOWSIZE(source, event)
        set(winSizeH, 'Enable', 'off');
        drawnow;

        val = source.String;
        winSize = str2double(val);

        uiresume(f);
    end

    function fcn_STEPSIZE(source, event)
        set(stepSizeH, 'Enable', 'off');
        drawnow;

        val = source.String;
        stepSize = str2double(val);

        uiresume(f);
    end

    function fcn_SMWINSIZE(source, event)
        set(smoothWinSizeH, 'Enable', 'off');
        drawnow;

        val = source.String;
        smoothWinSize = str2double(val);

        uiresume(f);
    end

    function fcn_THRCH(source, event)
        set(thrChH, 'Enable', 'off');
        drawnow;

        val = source.String;
        thrCh = str2double(val);

        uiresume(f);
    end

    function fcn_THRDUR(source, event)
        set(thrDurH, 'Enable', 'off');
        drawnow;

        val = source.String;
        thrDur = str2double(val);

        uiresume(f);
    end
        
    function fcn_GreyBar(source, event)
        set(greyBarH, 'Enable', 'off');
        drawnow;
        
        greyBar = get(greyBarH, 'value');

        uiresume(f);
    end


end