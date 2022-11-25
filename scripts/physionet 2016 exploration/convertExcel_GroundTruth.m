% in this script I import the ground truth data for the physionet 2016
% challenge training set into tables.

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["recordName", "class"];
opts.VariableTypes = ["string", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "recordName", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "recordName", "EmptyFieldRule", "auto");

% Import the data
Ya = readtable("C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\physionet 2016 training data\training-a\REFERENCE.csv", opts);
Yb = readtable("C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\physionet 2016 training data\training-b\REFERENCE.csv", opts);
Yc = readtable("C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\physionet 2016 training data\training-c\REFERENCE.csv", opts);
Yd = readtable("C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\physionet 2016 training data\training-d\REFERENCE.csv", opts);
Ye = readtable("C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\physionet 2016 training data\training-e\REFERENCE.csv", opts);
Yf = readtable("C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\physionet 2016 training data\training-f\REFERENCE.csv", opts);

Ygt.a = Ya;
Ygt.b = Yb;
Ygt.c = Yc;
Ygt.d = Yd;
Ygt.e = Ye;
Ygt.f = Yf;
%% Clear temporary variables
clear opts