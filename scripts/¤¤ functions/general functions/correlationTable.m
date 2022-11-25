function [Tcorr,pVal] = correlationTable(T,variables,round2,calcpVal,includePvalInTable)
% Takes a table T and produces a table with correlations between each pair
% of variables in 'variables'.
% X = zeros(height(T), numel(variables));
if nargin==2
    round2 = [];
    calcpVal = false;
    includePvalInTable = false;
elseif nargin==3
    calcpVal = false;
    includePvalInTable = false;
elseif nargin==4
    includePvalInTable = false;
end

% extract subtable:
T = subtable(T,variables);
% convert to double:
T = convertvars(T,variables,'double');

if calcpVal
    [X,pVal] = mycorr(table2array(T),[],true);
else
    X = mycorr(table2array(T));
    pVal = [];
end

if not(isempty(round2))
    X = round(X,round2);
end

Tcorr = array2table(X,'V',variables,'R',variables);

if includePvalInTable
    Nrows = height(Tcorr);
    Ncols = width(Tcorr);
    T = cell(Nrows,Ncols);
    for i=1:Nrows
        for j=1:Ncols
            stars = getPvalStars(pVal(i,j));
            if contains(stars,' ')
                T{i,j} = sprintf('%g',Tcorr{i,j});
            else
                T{i,j} = sprintf('%g (%s)',Tcorr{i,j},stars);
            end
        end
    end
    Tcorr = cell2table(T,'Row',Tcorr.Properties.RowNames,'Var',Tcorr.Properties.VariableNames);
end

end