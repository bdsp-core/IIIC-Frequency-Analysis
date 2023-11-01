function [ newResult ] = makeNanTable(thisResult)
%makeNanTable When Matlab creates an empty table, strangely, you can't
%modify it (because it's empty and any operation requires that you maintain
% the same 'size' of the table (which is zero wrt rows). Workaround is to
% generate a new table with all of the same columns, but with nan inserted
% to give it size 1 in the 1st dimension. 

theseVars = tableNames(thisResult);

newResult = table;

for j=1:length(theseVars)
    newResult.(theseVars{j}) = nan;
end


end

