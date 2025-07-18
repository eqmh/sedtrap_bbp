% % This script loads and appends raw data tables from ECO-Triplet sensors (Sea-Bird)
% mounted on sediment traps 

% % Writen by: Enrique Montes
% % Date: February 18 2025

clear all

cd('~/Documents/sedtrap_bbp');
dirPath = strcat(pwd, '/data_gom/raw');

% Get a list of files in the directory
fileList = dir(fullfile(dirPath, '*.raw'));

% Display the file names
for k = 1:length(fileList)
    disp(fileList(k).name);
end

% Initialize an empty table to store the combined data
conc_date = [];
conc_time = [];
conc_b532 = [];
conc_b700 = [];
conc_b695 = [];

% Loop over each file
for k = 1:length(fileList)
    % Get the full path of the file
    filePath = fullfile(fileList(k).folder, fileList(k).name);
    
    % Open the file
    fid = fopen(filePath, 'r');
    
    % Read the first row (string)
    firstRow = fgetl(fid);
    disp(['Date ranges in ', fileList(k).name, ': ', firstRow]);
    
    % Close the file
    fclose(fid);
    
    % Convert the cell array to a table
    [tmp_date,tmp_time,tmp_b532,tmp_b700,tmp_b695] = raw_import(filePath, inf);

    % Delete first and last rows in each array
    tmp_date = tmp_date(2:end-1);
    tmp_time = tmp_time(2:end-1);
    tmp_b532 = tmp_b532(2:end-1);
    tmp_b700 = tmp_b700(2:end-1);
    tmp_b695 = tmp_b695(2:end-1);
    
    % Concatenate each array
    conc_date = [conc_date; tmp_date];
    conc_time = [conc_time; tmp_time];
    conc_b532 = [conc_b532; tmp_b532];
    conc_b700 = [conc_b700; tmp_b700];
    conc_b695 = [conc_b695; tmp_b695];
    
end

%Creates datetime arrays for each dataset: merges dates and times and convert to datetime format
%Sensor 1
dateval = [];
for i=1:length(conc_date)
    conc_str = strcat(conc_date(i),{' '},conc_time(i));
    conc_num = datenum(conc_str);
    dateval = [dateval; conc_num];
end

date_time = datetime(dateval,'ConvertFrom','datenum');

combined_tbl = table(date_time, conc_b532, conc_b700, conc_b695);
combined_tbl = sortrows(combined_tbl, 'date_time');

figure()
subplot(2,1,1)
    plot(combined_tbl.date_time,combined_tbl.conc_b532)
    title('ECO 532nm - raw')
    ylabel('Counts','FontSize', 10)
    ylim([0 4500])
subplot(2,1,2)
    plot(combined_tbl.date_time,combined_tbl.conc_b700)
    title('ECO 700nm - raw')
    ylabel('Counts','FontSize', 10)
    ylim([0 4500])

% Change the column headers
% combined_tbl.Properties.VariableNames = {'date', 'time', 'bbp_532', 'bbp_700', 'chl_695'};
% writetable(combined_tbl, 'eco_counts_time_series.csv');

% % %%%%%%%%%%%%% FOR Bbp data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a list of all CSV files in the current directory
cd('~/Documents/sedtrap_bbp/data_sbb/');
bbp_files = dir('bbp_*.tsv');

% Initialize an empty table
bbp_tbl = table();

% Loop through each file and append its contents
for i = 1:length(bbp_files)
    
    % --- Read the current TSV file ---
    % First, detect the import options for a delimited text file.
    opts = detectImportOptions(bbp_files(i).name, 'FileType', 'text');
    
    % Now, specify that the delimiter is a tab ('\t').
    opts.Delimiter = '\t';
    
    % Read the table using these specific options.
    tempTable = readtable(bbp_files(i).name, opts);
    
    % Append to the combined table
    bbp_tbl = [bbp_tbl; tempTable]; % Concatenation
end

bbp_tbl_sorted = sortrows(bbp_tbl, 'datetime');

figure(2)
subplot(2,1,1)
    plot(bbp_tbl_sorted{:,1},bbp_tbl_sorted{:,2})
    title('532nm')
    ylabel('b_{bp} (m^{-1})','FontSize', 22)
    ylim([0 0.2])
subplot(2,1,2)
    plot(bbp_tbl_sorted{:,1},bbp_tbl_sorted{:,3})
    title('700nm')
    ylabel('b_{bp} (m^{-1})','FontSize', 22)
    ylim([0 0.2])

% --- Write Data to CSV (Optional) ---
% disp('Writing data to bbp_*.tsv...');
% writetable(bbp_tbl_sorted, 'eco_bbp_time_series.tsv', 'FileType', 'text', 'Delimiter', '\t');
% disp('Done.');
