function [date,time,b532,b700,b695] = raw_import(file1, inf)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [VARNAME1,VARNAME2,VARNAME4,VARNAME6,VARNAME8] = IMPORTFILE(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   [VARNAME1,VARNAME2,VARNAME4,VARNAME6,VARNAME8] = IMPORTFILE(FILENAME,
%   STARTROW, ENDROW) Reads data from rows STARTROW through ENDROW of text
%   file FILENAME.
%
% Example:
%   [VarName1,VarName2,VarName4,VarName6,VarName8] =
%   importfile('BB2FLWB-1314_May2017.raw',1, 62825);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/07/25 11:27:44

%% Initialize variables.
delimiter = '\t';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column4: double (%f)
%	column6: double (%f)
%   column8: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%*s%f%*s%f%*s%f%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(file1,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
date = dataArray{:, 1};
time = dataArray{:, 2};
b532 = dataArray{:, 3};
b700 = dataArray{:, 4};
b695 = dataArray{:, 5};


