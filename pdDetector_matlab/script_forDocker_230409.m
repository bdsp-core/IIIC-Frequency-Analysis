function [errorCode] = script_forDocker_230409(targetDir, outputDir, saveMAT, savePNG)
%% Script for Docker
% INPUT
%   targetDir, with a trailing slash
%   outputDir
%
%% Choose Target directory

% Display a dialog box for the user to select a directory
% selected_directory = uigetdir();
%
% Check if the user selected a valid directory
% if selected_directory ~= 0
%     fprintf('Selected directory: %s\n', selected_directory);
% else
%     fprintf('No directory was selected.\n');
% end
%
% theseDir = {'/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Hard/', ...
%     '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Medium/',...
%     '/data/01 Code/MATLAB/05 PD/04 Data from Jing 9-2022/IIICData4Chris/Easy/'};
%
% saveDir = '';

if strcmp(targetDir,'') || strcmp(outputDir,'') || strcmp(saveMAT,'') || strcmp(savePNG,'') 
disp('please call function with: targetDir, outputDir, saveMAT, savePNG');
disp('targetDir, the directory location containing .MAT files with EEG segments');
disp('outputDir, the directory to save analyses');
disp('saveMAT, indicate if MAT files should be saved (1= true or 0=false)');
disp('savePNG, indicate if PNG files should be saved (1= true or 0= false)');
exit  
else
disp(['Target directory: ', targetDir, ', Output directory: ', outputDir]);
end

%% Create directory structure
if ~strcmp(targetDir(end),'/') 
    targetDir = [targetDir '/'];
end

if ~strcmp(outputDir(end),'/') 
    outputDir = [outputDir '/'];
end

% addpath(targetDir);
% addpath(outputDir);

formatOut = 'yy-mm-dd';
thisDate = datestr(now,formatOut);
% directory structure
% >thisDir            a directory containing AVI files
% >> avi files        to be processed
% >> analysis_<date>    subDir to place CSV/XLS output from analysis
% >>> mat               subDir to place MAT files
% >>> images            subDir to place figures/AVI output

thisAnalysisDir = [outputDir 'analysis_' thisDate '/'];
thisImageDir = [thisAnalysisDir 'images/'];
thisMATDir = [thisAnalysisDir 'mat/'];
make_directories({thisAnalysisDir, thisImageDir, thisMATDir});
% addpath(thisAnalysisDir);
% addpath(thisImageDir);
% addpath(thisMATDir);

%% Iterate over files

files = dir([targetDir '*.mat']);

resultTable = table;

tic
for n=1:length(files)
    
    % for testing
    %if ~strcmp(files(n).name, 'sid1489_20140905_074518_10912.mat')
    %    continue
    %end

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
            
            if savePNG==0
%                 display(savePNG)
%                 print(f1, '-dpng','-zbuffer','-r600', [thisImageDir files(n).name '.png']);
                print(f1, '-dpng','-r600', [thisImageDir files(n).name '.png']);
                close(f1);
            else
                close(f1);
            end

            if saveMAT==0
%                 display(saveMAT)
                save([thisMATDir files(n).name '.mat'],'data_obj','thisStatsObj');
            end

            thisResult2 = cat(2, thisResult,thisStatsObj);
            writetable(thisResult2,[thisMATDir 'result_' files(n).name '.csv']);
        end
        
    end
    %%
    
%     resultTable = cat(1,resultTable,thisResult);
    
    
end
toc
% save([thisMATDir  'resultTable.mat'],'resultTable');

