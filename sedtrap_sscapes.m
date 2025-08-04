% % This program extracts seascape class from the location of the sediment
% trap deployment in the Gulf of Mexico

clear all; close all; 

% Load data and mapping dependencies:
cd('~/Documents/sedtrap_bbp/data_gom/');
addpath(genpath('/Users/enrique.montes/Library/CloudStorage/GoogleDrive-enriquemontes01@gmail.com/My Drive/GDrive/software/matlab/m_map'));

%Loads data tables
fileID = 'trap_flux_gom.csv';
data = readtable(fileID); 

% Set time window
buffer = 1; % number of months before and after start and end dates
start_year = num2str(year(min(data.Date)));
start_month = num2str(month(min(data.Date)) - buffer);
end_year = num2str(year(max(data.Date)));
end_month = num2str(month(max(data.Date)) + buffer);

start_date = strcat(start_year,'-',start_month);
end_date = strcat(end_year,'-',end_month);

% Sediment trap coordinates
sedtrap_lon = -89;
sedtrap_lat = 28;

% Study area limits
n_lat = '32'; s_lat = '24'; e_lon = '-85'; w_lon = '-93';

% % Getting the data from ERDDAP
% 8-DAY
% urlfile=strcat('https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_8DAY/SEASCAPES.nc?var=CLASS&var=P&north=',n_lat,'&west=',w_lon,'&east=',e_lon,'&south=',s_lat,'&disableProjSubset=on&horizStride=1&time_start=',start_date,'-15T12%3A00%3A00Z&time_end=',end_date,'-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf');

% % Getting the data from ERDDAP
% % MONTHLY
urlfile=strcat('https://cwcgom.aoml.noaa.gov/thredds/ncss/SEASCAPE_MONTH/SEASCAPES.nc?var=CLASS&var=P&north=',n_lat,'&west=',w_lon,'&east=',e_lon,'&south=',s_lat,'&disableProjSubset=on&horizStride=1&time_start=',start_date,'-15T12%3A00%3A00Z&time_end=',end_date,'-15T12%3A00%3A00Z&timeStride=1&addLatLon=true&accept=netcdf');

options = weboptions('CertificateFilename', '', 'Timeout',180); % this command is necessary in how matlab looks for bigger data files
websave('dummy.nc',urlfile,options); % note in matlab this saves the url exactly as the file that it is, not a text or html
lat=double(ncread('dummy.nc','lat'));
lon=double(ncread('dummy.nc','lon'));
CLASS=double(ncread('dummy.nc','CLASS'));
P=double(ncread('dummy.nc','P'));
time=double(ncread('dummy.nc', 'time'));
[lat1,lon1]=meshgrid(lat,lon);

%time is seconds since Jan 01, 1970 in Coastwatch world
dnum=datenum([1970,1,1])+time./(60*60*24) ; % need to convert seconds to days
TIMEVEC=datevec(dnum);
doy=nan(length(time),1);
for q=1:length(doy);
    if TIMEVEC(q,2)==1;
        doy(q,1)=TIMEVEC(q,3);
    else;
        doy(q,1)=sum(eomday(TIMEVEC(q,1),1:TIMEVEC(q,2)-1),2)+TIMEVEC(q,3); 
        %eomday provides the total number days in a month and accounts for
        %leap years. This is an easy way to calculate the day of the year
        %so that we can calculate a series date for easy plotting
    end;
end;

%%create the series date
SDATE=TIMEVEC(:,1)+doy./366;

% % Extract seascapes at sediment trap location
seascape_list = [];
n_pix = 3; %CHANGE HERE
pix_res = 24; %CHANGE HERE: USE 110 WITH 1-KM AND 24 WITH 4-KM SEASCAPES
inx = find(distance(sedtrap_lat, sedtrap_lon, lat1(:), lon1(:))<(n_pix/pix_res));  
seaslat = lat1(inx);
seaslon = lon1(inx);
DDIST=distance(sedtrap_lat,sedtrap_lon, lat1(inx),lon1(inx)); %Distance from station
seascape_count_matrix = zeros(size(TIMEVEC, 1), 33);
for j=1:size(TIMEVEC, 1)

    % extract seascape classes surrounding stations
    selected_seascape = CLASS(:,:,j);
    selected_probability = P(:,:,j);
    seastype = selected_seascape(inx);
    seas_prob = selected_probability(inx);
    
    % Calculate seascape mode 
    seascape_mode = mode(seastype(~isnan(seastype)));
    
    % Calculate the mean seas_prob using rows containing seascape_mode values only
    valid_indices = seastype == seascape_mode; % Find indices where seastype matches seascape_mode
    mean_seas_prob = mean(seas_prob(valid_indices), 'omitnan'); % Calculate mean, ignoring NaN values

    seascape_list = [seascape_list; TIMEVEC(j,1), TIMEVEC(j,2), TIMEVEC(j,3), SDATE(j), seascape_mode, mean_seas_prob];

    % Count the number of pixels for each seascape class (1 to 33)
    for class_id = 1:33
        seascape_count_matrix(j, class_id) = sum(seastype(:) == class_id, 'omitnan');
    end

end

% Use this chunk to check pixels from which seascape class is extracted
% per station (need to generate seascape map first). NOTE: Cross is placed in
% the lower left corner of the pixel.
cmap1 = csvread('cmap1.csv'); % color palette for seascapes
figure();
    m_proj('robinson','lon',[str2num(w_lon) str2num(e_lon)],'lat',[str2num(s_lat) str2num(n_lat)]); 
    m_pcolor(lon1,lat1, CLASS(:,:,9));
    colormap(cmap1);
    m_coast('patch',[.7 .7 .7],'edgecolor','none');
    m_grid('tickdir','out','linewi',2);
    h=colorbar;
    scale = [1 33];
    caxis(scale);
hold on
    seas_box = find(distance(sedtrap_lat, sedtrap_lon, lat1(:), lon1(:)) < (n_pix/pix_res)); 
    for q=1:size(seas_box,1)
        xx2=lon1(seas_box(q));
        yy2=lat1(seas_box(q));
        [X,Y]=m_ll2xy(xx2,yy2);
        line(X,Y,'marker','+','markersize',7,'linewidth',1,'color','k');
    end
    [X,Y]=m_ll2xy(sedtrap_lon, sedtrap_lat);
    line(X,Y,'marker','o','markersize', 7,'linewidth', 0.5,'color','k','MarkerFaceColor','white');
hold off


% % Plot number of pixels per class within box surrounding trap location
% Normalize the seascape count matrix to a maximum of 20
normalized_count_matrix = (seascape_count_matrix ./ length(DDIST));

% --- Step 1: Aggregate your data into daily sums ---

% Create a table to make grouping easy. Let's create generic variable names.
varNames = "Class-" + (1:33);
T = array2table(normalized_count_matrix, 'VariableNames', varNames);

% Create the T.Date column from the first three columns of TIMEVEC
T.Date = datetime(TIMEVEC(:,1), TIMEVEC(:,2), TIMEVEC(:,3));

% --- Step 2: Create the stacked bar chart ---

figure();
% --- Left Y-Axis (Primary) ---
yyaxis left; % Activate the left axis FIRST
% Create the stacked bar chart
h = bar(T.Date, normalized_count_matrix, 'stacked', 'BarWidth', 1, 'EdgeColor', 'k', 'LineWidth', 0.5);
ylabel('Class fraction within box');

% Loop through each of the 33 bar series and apply the correct color.
for i = 1:length(h)
    h(i).FaceColor = cmap1(i, :);
end

% --- Right Y-Axis (Secondary) ---
yyaxis right; % NOW, activate the right axis
% Plot the total mass flux
plot(data.Date, data.org_c_flux, 'LineWidth', 2.5, 'Color', 'w'); % Adjust color as needed
% Set the label for the right axis
ylabel('Org. C flux');

% Set the axis colors to avoid confusion
ax = gca;
ax.YColor = 'w'; % Set right y-axis color to match the line plot
yyaxis left; % Switch back to left
ax = gca;
ax.YColor = 'w'; % Set left y-axis color to black (default)

% Add a legend. Using the original variable names for clarity.
legend(varNames, 'Location', 'eastoutside');

% Optional: If you have many dates, improve the x-axis tick labels
ax = gca;
ax.XTick = dateshift(min(T.Date), 'start', 'month'):calmonths(1):max(T.Date);
ax.XTickLabelRotation = 45; % Angle the labels to prevent overlap

% Convert seascape_list to a table with specified column headers
% seascape_table = array2table(seascape_list, 'VariableNames', {'year', 'month', 'day', 'date_vector', 'class', 'probability'});
% % Save table to a TSV file 
% writetable(seascape_table, 'seascape_data_8day.tsv', 'FileType', 'text', 'Delimiter', '\t');