function a = tableNames( thisTable)
%tableNames Return the variable names from a given table wo having to type
%out explicit reference to structure field Properties.VariableNames

a = thisTable.Properties.VariableNames;

end

