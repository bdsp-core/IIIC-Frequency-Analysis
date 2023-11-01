
theseDir = {'/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Hard/', ...
    '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Medium/',...
    '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Easy/'};

for j=1:length(theseDir)
    
     parentDir = theseDir{j};

    addpath(parentDir);
    
    formatOut = 'yy-mm-dd';
    thisDate = datestr(now,formatOut);
        % directory structure
    % >thisDir            a directory containing AVI files
    % >> avi files        to be processed
    % >> analysis_<date>    subDir to place CSV/XLS output from analysis
    % >>> mat               subDir to place MAT files
    % >>> images            subDir to place figures/AVI output
    
    thisAnalysisDir = [parentDir 'analysis_' thisDate '/'];
    thisImageDir = [thisAnalysisDir 'images/'];
    thisMATDir = [thisAnalysisDir 'mat/'];
    make_directories({thisAnalysisDir, thisImageDir, thisMATDir});
    addpath(thisAnalysisDir);
    addpath(thisImageDir);
    
    files = dir([parentDir '*.mat']);
    
    resultTable = table;
    
    tic
    for n=1:length(files)
        
        thisResult=table;
        S =load([files(n).folder '/' files(n).name]);
        % seg= S.data(1:19,:);  % remove EKG
        
        P_model = S.P_model;
        Y_human = S.Y_human;
        
        Y_human_numVotes = sum(Y_human{1,:});
        
        Y_human_p = array2table(Y_human{1,:}./sum(Y_human{1,:}),'VariableNames',Y_human.Properties.VariableNames);  % convert to percentage
        
        % determine which detection function to run based on human rating.
        Y_human_s = Y_human_p(:,Y_human_p{1,:}>=1/3); % all types >1/3
        theseTypes = Y_human_s.Properties.VariableNames;
        theseTypes_p = Y_human_s{1,:};
        
        thisResult.n = n;
        thisResult.file = {files(n).name};
        thisResult.y_human_types = {theseTypes};
        thisResult.y_human_types_prct = {theseTypes_p};
        thisResult.statsObj = {NaN};
        
        Y_human_s_m = Y_human_s(:,ismax(Y_human_s{1,:})); % top type
        thisType = Y_human_s_m.Properties.VariableNames;
        
        if ~isempty(thisType)
            
            switch thisType{:}
                case {'GRDA','LRDA'}
                    disp('grda,lrda');
                    resultTable = cat(1,resultTable,thisResult);
                    continue
                case {'LPD','GPD'}
                    disp('lpd, gpd');
                    
                    %      S =load([files(n).folder '/' files(n).name]);
                    seg= S.data(1:19,:);  % remove EKG
                    
                    [freq, channels, interpretation, data_obj] = fct_PD_detect(seg);
                    
                case 'Seizure'
                    disp('seizure');
                    seg= S.data(1:19,:);  % remove EKG
                    
                    [freq, channels, interpretation, data_obj] = fct_PD_detect(seg);
                    
                case 'Other'
                    disp('other');
                    resultTable = cat(1,resultTable,thisResult);
                    continue
            end
            
            %% Generate figure
            if isstruct(data_obj)
                data_obj.num = n;
                data_obj.filename = files(n).name;
                data_obj.theseTypes = theseTypes;
                data_obj.theseTypes_prctVotes = theseTypes_p;
                data_obj.Y_human = Y_human;
                data_obj.P_model = P_model;

                [f1, thisStatsObj] = fct_showEEG_wEventsAndStats_dataObj(data_obj);
%                 set(f1, 'Visible','on');
                % thisStatsObj contains the stats re: evolution 

                print(f1, '-dpng','-zbuffer','-r600', [thisImageDir files(n).name '.png']);
                close(f1);
                
                save([thisMATDir files(n).name '.mat'],'data_obj','thisStatsObj');
                thisResult.statsObj = {thisStatsObj};
                
            end
            
        end
        %%
        
        resultTable = cat(1,resultTable,thisResult);
        
        
    end
    toc
    save([thisMATDir  'resultTable.mat'],'resultTable');

end