function [ output_args ] = make_directories( theseDirecs )
%make_directories Check if directories have been created and if not, make
%them.

if ~isempty(theseDirecs)
    
    for i=1:length(theseDirecs)
        thisDirec = theseDirecs{i};
        
        if isempty(dir(thisDirec))
            disp([thisDirec ' created.']);
            mkdir(thisDirec);
        else
            disp([thisDirec ' already exists.']);
        end
    end
    
end


end

