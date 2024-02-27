
% STEP 3
%% PLEASE NOTE:
% THIS CODE REQUIRES a Licensed BARON SOLVER!
% BARON itself uses IBM CPLEX which also requires a license. 
% While BARON Academic license can be acquired for a fee,
% CPLEX solver does provide academic free license easily obtainable.
%%
% Clear the workspace, console, and close all figures
clear, clc, close all
% Load the network data from a MAT file
load('../model/sphingolipid_network.mat')
% Set parameters for DMFA analysis
lambdaOpt = 1e-4;
N = 2;
sdFracScale = 0;
tolerance = 1e-3;
robustFit = false;
for datasetId = [3, 4] % Choose either 3(C) or 4(D) 
% Determine the dataset identifier for 'v'
if datasetId == 3
    v = 'C';
elseif datasetId == 4
    v = 'D';
end
% Set the random number generator to its default state
rng default
% Specify the input Excel file
inputFilename = '../data/N15_dynamic_data.xlsx';
% Determine the MATLAB version and read data accordingly
[version_, date] = version;
if year(date) >= 2020
    fileSheetNames = sheetnames(inputFilename);
elseif year(date) >= 2011
    fileSheetNames = xl_xlsfinfo(fullfile(pwd, inputFilename));
else
    fprintf('Code only works for versions of MATLAB released 2011b or after.\n')
end
% Select the sheet based on the dataset identifier
selectedSheet = fileSheetNames{datasetId};
% Read data from the selected sheet
if year(date) >= 2016
    tmpTbl = readtable(inputFilename, 'Sheet', selectedSheet);
else
    fprintf('Code only works for versions of MATLAB released 2016b or after.\n')
end
% Extract measured metabolites and concentration data
measuredMetabolites = table2cell(tmpTbl(:, 2));
nonEmptyIdx = find(not(cellfun('isempty', measuredMetabolites)));
measuredMetabolites = measuredMetabolites(nonEmptyIdx(1:end));
numMetabolites = numel(measuredMetabolites);
concentrationData = table2cell(tmpTbl(:, 3:end));
concentrationData_ = concentrationData;
nonEmptyIdx = find(not(cellfun('isempty', concentrationData_)));
concentrationValues = cell2mat(concentrationData_(nonEmptyIdx));
minConcentration = min(min(concentrationValues));
maxConcentration = max(max(concentrationValues));
concentrationValues = reshape(normalize(concentrationValues, 'range'), size(concentrationData_)) + 1e-8;
concRawData = concentrationValues;
stdMeasurementRaw = concRawData(:, 8:end);
concentrationData = concRawData(:, 1:7);
stdMeasurement = stdMeasurementRaw;
if datasetId == 4
    stdMeasurement = stdMeasurement(:, 1:6);
    concentrationData = concentrationData(:, 1:6);
end
% Initialize timing
StartRunTime_all = tic;
numSamples = 1;
% Preallocate arrays to store results
flux_data_n_all = cell(1, numSamples);
flux_data_sd_n_all = cell(1, numSamples);
flux_data_all = cell(1, numSamples);
flux_data_sd_all = cell(1, numSamples);
EndRunTime_all = zeros(1, numSamples);
nDMFA_all = zeros(1, numSamples);
data_all = cell(1, numSamples);
tn_all = cell(1, numSamples);
fit_ssr_all = zeros(1, numSamples);
accum_rate_all = cell(1, numSamples);
SampleName_all = cell(1, numSamples);
% Initialize valid sample indices
validSampleIndices = [];
sample = 1;
% Perform DMFA analysis and store results
output_sample = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minConcentration, maxConcentration, N, lambdaOpt, sdFracScale, robustFit, tolerance, []);
% Store reference value for nDMFA
reference_nDMFA = output_sample.nDMFA;
% Store the results for the current sample
flux_data_n_all{sample} = output_sample.flux_data_n;
flux_data_sd_n_all{sample} = output_sample.flux_data_sd_n;
flux_data_all{sample} = output_sample.flux_data;
flux_data_sd_all{sample} = output_sample.flux_data_sd;
EndRunTime_all(sample) = output_sample.EndRunTime;
nDMFA_all(sample) = output_sample.nDMFA;
data_all{sample} = output_sample.data;
tn_all{sample} = output_sample.tn;
fit_ssr_all(sample) = output_sample.fit_ssr;
accum_rate_all{sample} = model.S * output_sample.flux_data;
SampleName_all{sample} = ['WT_' v];
% Add the sample index to the validSampleIndices array
validSampleIndices = [validSampleIndices, sample];
% Create a cell array to store data
dataCell = {flux_data_n_all, flux_data_sd_n_all, flux_data_all, flux_data_sd_all, ...
    EndRunTime_all, nDMFA_all, data_all, tn_all, fit_ssr_all, accum_rate_all, SampleName_all, validSampleIndices'};
% Define field names for the structure
fieldNames = {'flux_data_n_all', 'flux_data_sd_n_all', 'flux_data_all', 'flux_data_sd_all', ...
    'EndRunTime_all', 'nDMFA_all', 'data_all', 'tn_all', 'fit_ssr_all', 'accum_rate_all', 'SampleName_all', 'ValidSampleIndices'};
% Convert the cell array to a structure
resultData = cell2struct(dataCell, fieldNames, 2);
% Extract relevant data for bounds calculation
conc_data_ = output_sample.meas_raw_x_vals;
meas_est_x_vals = output_sample.meas_est_x_vals;
meas_est_x_vals_sd = output_sample.meas_est_x_vals_sd;
idxsDays = output_sample.indicesDays;
% Calculate concentration bounds
meas_est_x_vals_ub = meas_est_x_vals + meas_est_x_vals_sd;
meas_est_x_vals_lb = meas_est_x_vals - meas_est_x_vals_sd;
meas_est_x_vals_lb(meas_est_x_vals_lb < 0) = 1e-6;
conc_pred_bounds = [min(meas_est_x_vals_lb(:, idxsDays(3:end)), [], 2) max(meas_est_x_vals_ub(:, idxsDays(3:end)), [], 2)];
% Save the sampling settings and results to a MAT file
EndRunTime = toc(StartRunTime_all) / 60;
outputfilename = sprintf('../results/%s_DynFit_WT.mat', v);
save(outputfilename, 'output_sample', 'model', 'resultData', 'numSamples', 'lambdaOpt', 'N', 'datasetId', 'sdFracScale', 'tolerance', 'robustFit', 'EndRunTime', 'v', ...
    'measuredMetabolites', 'conc_pred_bounds', 'meas_est_x_vals', 'meas_est_x_vals_sd', 'meas_est_x_vals_lb', 'meas_est_x_vals_ub', 'idxsDays', '-v7.3');
% Display a message indicating successful data saving
fprintf('\nData saved to file: %s\n', outputfilename);

end
%%
clear
load('../results/C_DynFit_WT.mat')
WT_flux_C = resultData.flux_data_n_all{1, 1}  ;
load('../results/D_DynFit_WT.mat')
WT_flux_D = resultData.flux_data_n_all{1, 1}  ;
save('../results/WT_fluxes_C_D.mat','WT_flux_C', 'WT_flux_D');
close all
%
% load('../results/WT_fluxes_C_D.mat')
WT_flux_C = WT_flux_C(:, 1:end-1);
% Initialize matrices to hold the means and standard deviations
mean_profiles = zeros(size(WT_flux_C)); % 40 rows for items, 7 columns for time points or conditions
std_profiles = zeros(size(WT_flux_C));
% Loop over each set of 40 rows to calculate the mean and standard deviation
for i = 1:size(WT_flux_C, 1) % 
    % Calculate mean and standard deviation for the current block for WT_flux_C
    mean_profiles(i, :) = mean([WT_flux_C(i, :); WT_flux_D(i, :)], 1);
    std_profiles(i, :) = std([WT_flux_C(i, :); WT_flux_D(i, :)], 0, 1);
    if ~any(std_profiles(i, :))
        std_profiles(i, :) = .025*ones(size(std_profiles(i, :)));
    end
end
% Now mean_profiles and std_profiles contain the means and standard deviations
% of WT_flux_C and WT_flux_D combined for each block of 40 rows.
% Load the network data from a MAT file
load('../model/sphingolipid_network.mat')
load('color_bar_id.mat')
rxnIDs = model.rxnIDs ;
figNo = 0;
V_predFitted = mean_profiles;
V_predFitted_sd = std_profiles;
time = [0:5];
rxn_names = [];
for j = 1:numel(rxnIDs)
    tmp_labels = rxnIDs(j);
    h_tmp = figure(figNo + 1);
    tmp_rxn_data = V_predFitted(j, :);
    tmp_rxn_data_sd = V_predFitted_sd(j, :);
    %  Bar graph
    errlow = tmp_rxn_data_sd;
    errhigh = tmp_rxn_data_sd;
    figure(figNo+1)
    bar(time,tmp_rxn_data, 'FaceColor', color_id_rxn2{j})
    hold on
    plot(time, tmp_rxn_data, 'ko', 'MarkerFaceColor', 'k')
    hold on
    er = errorbar(time,tmp_rxn_data,errlow,errhigh);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    er.LineWidth = 3;
    hold on    
    xline(2-0.5,'--g', 'linewidth', 2.0) ;
    xline(5+0.5,'--g', 'linewidth', 2.0) ;
    hold off
    set(gca,'box','off') 
    a_x = get(gca,'XTickLabel');
    fontsize = 36;
    set(gca,'XTickLabel',a_x,'FontName','Times','fontsize',fontsize)
    a_y = get(gca,'XTickLabel');
    fontSize = 30; 
    set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
    fig_name_plot = fullfile('../results','CodeFigures', 'Fig_4_NormFluxBarPlots_C_D_Avg', [char(rxnIDs(j)) '.png']);
%     fig_name_plot = fullfile(pwd, 'NormFluxStats_Barplots', 'C_D_Avg', [char(rxnIDs(j)) '.svg']);
    folder_path = fileparts(fig_name_plot);  % Extract the folder path from the full file name
    if exist(folder_path, 'dir') ~= 7
        mkdir(folder_path);  % Create the folder if it doesn't exist
    end
    saveas(gcf, fig_name_plot);
    rxn_names = [rxn_names; rxnIDs(j)];
end



