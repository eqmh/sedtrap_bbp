% % This script process raw data from ECO-Triplet sensors (Sea-Bird)
% for raw count data to Bbp and creates plots for data visualization and
% exploration.

% % Writen by: Enrique Montes
% % Date: July 17 2025

clear all

% % 1. Configuration and Setup
% % Define constants and parameters 

% --- Input File ---
filename = 'BB2FLWB-1315_May2017.raw';

% --- Data Filtering Parameters ---
% Note: The original script uses 'MM/dd/yy' format.
startDate = '05/10/17'; 
endDate   = '11/03/17'; 

% --- Calibration Coefficients (for BB2FLWB-1315, 03/20/2023) ---
% BB2FLWB calibration coefficients:
% 06/25/2015:
% BB2FLWB-1314: SF(532) = 6.999e-6 / Dark counts = 55
%               SF(700) = 3.037e-6 / Dark counts = 50
% BB2FLWB-1315: SF(532) = 6.946e-6 / Dark counts = 50
%               SF(700) = 3.002e-6 / Dark counts = 43
% 03/20/2023:
% BB2FLWB-1314: SF(532) = 7.419e-6 / Dark counts = 55
%               SF(700) = 3.211e-6 / Dark counts = 53 

cal.sf_532 = 6.946e-6;   % Scale factor for 532nm
cal.off_532 = 50;        % Dark counts for 532nm
cal.sf_700 = 3.002e-6;   % Scale factor for 700nm
cal.off_700 = 43;        % Dark counts for 700nm

% --- Physical Parameters for Water Scattering Calculation ---
phys.salinity = 34;      % Salinity (S)
phys.tempC = 9;         % Seawater temperature (Celsius): 9 C in SBB @200m (Glodap: https://explore.webodv.awi.de/ocean/hydrography/glodap/glodapv2.2023/)
phys.angle = 124;        % Instrument scattering angle (degrees)
phys.delta = 0.039;      % Depolarization ratio (constant)
phys.chi_factor = 1.1;   % Proportionality factor for bbp calculation (from "WET Labs ECO BB User guide.pdf" page 16)

%% 2. Data Loading and Pre-processing
% Load the raw data and perform all filtering in a single, efficient step.

% --- Load Data ---
% The 'inf' argument suggests reading the entire file.
addpath('~/Documents/sedtrap_bbp')
[dateStr, timeStr, b532, b700, b695] = raw_import(filename, inf);

% --- Create Datetime Vector ---
% The format 'MM/dd/yy' is inferred from your original 'ismember' calls.
fullDateTime = datetime(strcat(dateStr, " ", timeStr), 'InputFormat', 'MM/dd/yy HH:mm:ss');

% --- Create a Single Logical Mask for Filtering ---
% Build one index of valid data points
isValid = true(size(fullDateTime));

% Mark first and last rows for removal (replicates original script's logic)
isValid(1) = false;
isValid(end) = false;

% Mark rows with NaNs in any of the measurement columns for removal.
% 'any' operates along the second dimension (rows) to find any row with a NaN.
nanRows = any(isnan([b532, b700, b695]), 2);
isValid(nanRows) = false;

% Mark rows outside the desired date range for removal.
% This is much cleaner and more reliable than using find/ismember.
startDt = datetime(startDate, 'InputFormat', 'MM/dd/yy');
endDt   = datetime(endDate, 'InputFormat', 'MM/dd/yy');
% By default, 'isbetween' creates an inclusive interval [startDt, endDt].
isValid(~isbetween(fullDateTime, startDt, endDt)) = false;

% --- Apply the Filter Mask ---
% Apply the single logical mask to all data arrays at once. This is far
% more efficient than deleting elements from arrays sequentially.
date_final = fullDateTime(isValid);
b532 = b532(isValid);
b700 = b700(isValid);
% b695 is also filtered in case it's needed later.
b695 = b695(isValid);

%% 3. Raw Data Visualization (Optional Check)
% Plot the cleaned, raw data to ensure filtering worked as expected.

figure();
clf; % Clear the figure before plotting

subplot(2, 1, 1);
plot(date_final, b532);
title('532nm - Raw Counts');
ylabel('Counts', 'FontSize', 10);
ylim([0 4500]);
grid on;

subplot(2, 1, 2);
plot(date_final, b700);
title('700nm - Raw Counts');
ylabel('Counts', 'FontSize', 10);
ylim([0 4500]);
grid on;

%% 4. Physical Calculations
% All calculations are vectorized and benefit from the clean data.

% --- Apply Calibration ---
% Calculate volume scattering coefficient (beta)
B532 = cal.sf_532 * (b532 - cal.off_532);
B700 = cal.sf_700 * (b700 - cal.off_700);

% --- Calculate Volume Scattering of Pure Seawater (Bw) ---
% Using the betasw_ZHH2009 function (assumed to be on MATLAB path).
[Bw_532, ~, ~] = betasw_ZHH2009(532, phys.tempC, phys.angle, phys.salinity, phys.delta);
[Bw_700, ~, ~] = betasw_ZHH2009(700, phys.tempC, phys.angle, phys.salinity, phys.delta);

% --- Calculate Particle Scattering (Bp) ---
Bp532 = B532 - Bw_532;
Bp700 = B700 - Bw_700;

% --- Calculate Particle Backscattering (bbp) ---
bbp_532 = 2 * pi * phys.chi_factor * Bp532;
bbp_700 = 2 * pi * phys.chi_factor * Bp700;

%% 5. Final Visualization and Output

% --- Plot Final Processed Data ---
figure();
clf; % Clear the figure

subplot(2, 1, 1);
plot(date_final, bbp_532);
title('Particle Backscattering at 532nm');
ylabel('b_{bp} (m^{-1})', 'FontSize', 10);
ylim([0 0.2]);
grid on;

subplot(2, 1, 2);
plot(date_final, bbp_700);
title('Particle Backscattering at 700nm');
ylabel('b_{bp} (m^{-1})', 'FontSize', 10);
ylim([0 0.1]);
grid on;

% --- Write Data to TSV (Optional) ---
disp('Writing data to bbp_*.tsv...');
bbp_table = table(date_final, bbp_532, bbp_700, ...
    'VariableNames', {'datetime', 'bbp_532nm', 'bbp_700nm'});
writetable(bbp_table, 'bbp_1315_20170510.tsv', 'FileType', 'text', 'Delimiter', '\t');
disp('Done.');
