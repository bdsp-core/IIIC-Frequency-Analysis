function [ thisTable ] = appendTblNames( thisTable, thisLabel)
%appendTblNames Add a fixed string to the end of table column names in table to
%differentiate from pre-existing column names
% INPUT
%
% OUTPUT
% 

        theseVarNames = thisTable.Properties.VariableNames;
           
           % append 'raw' to variable names from raw eeg (vs. values from
           % the transformed data)
           for k=1:length(theseVarNames)
              thisVarName = theseVarNames{k};
              theseVarNames(k) = {[thisVarName thisLabel]};           
           end
           
           thisTable.Properties.VariableNames = theseVarNames;

end

