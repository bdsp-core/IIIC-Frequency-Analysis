%%
% Take file, seg. 
% dataDirec = 'C:\Users\Chris McGraw\Dropbox (Partners HealthCare)\01 Data\2018\FrequencyOfDischarges_Chris_JJ\code\data\';
dataDirec = 'C:\Users\Chris McGraw\Dropbox (Partners HealthCare)\01 Data\2019\EEG data from Jing\'
outputDirec = [dataDirec 'output\'];
make_directories({outputDirec});

files = dir([dataDirec '*.mat']);

for i=1:length(files)
thisFile = files(i);
filename = thisFile.name;

[seg, seg_bi] = loadFile(files,i);
seg_avg = seg-mean(seg,1);

% save bipolar
fig = showEEG(seg_bi);
set(fig, 'Position', [-1508 905 893 641], 'Visible', 'off')
set(fig, 'PaperPositionMode','auto');
ax = gca;
set(ax, 'Position', [.05 .05 .93 .9]);
print(fig, '-dpng','-zbuffer','-r200', [outputDirec filename '_bi.png']);
close(fig);

% save average
fig = fct_showEEG_avg(seg_avg);
set(fig, 'Position', [-1508 905 893 641], 'Visible', 'off')
set(fig, 'PaperPositionMode','auto');
ax = gca;
set(ax, 'Position', [.05 .05 .93 .9]);
print(fig, '-dpng','-zbuffer','-r200', [outputDirec filename '_avg.png']);
close(fig);

end


% 
