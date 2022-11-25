function subT = subtable(motherTable,namesSubtable)
% Allows you to form a subtable of a larger table, consisting
% of the names in the string array namesSubtable.

subT = table;
for i=1:length(namesSubtable)
    subT = [subT, table(motherTable.(namesSubtable{i}),'v',namesSubtable(i))]; %#ok<*AGROW>
end

subT.Properties.RowNames = motherTable.Properties.RowNames;

end