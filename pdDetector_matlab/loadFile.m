function [ seg, seg_bi ] = loadFile( files, n )
%loadFile 
%             nameOfField = 'eeg';   %v2  need to specify the object in stored MAT that contains EEG data.
            nameOfField = 'seg'; %v1

            file = files(n);
            thisFilename=[file.folder '/' file.name];
            tmp = load(thisFilename);
            disp(['loadFile.Loading file: ' thisFilename]);
            
            
%             disp('loadFile.Demean and normalize by rms');
            seg = tmp.(nameOfField);
%             seg = (seg-mean(seg))./rms(seg);
            
            seg_bi = fcn_LBipolar(seg);
         
            

end

