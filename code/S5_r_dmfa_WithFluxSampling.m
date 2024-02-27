% STEP 5

%% PLEASE NOTE:
% THIS CODE REQUIRES a Licensed BARON SOLVER!
% BARON itself uses IBM CPLEX which also requires a license. 
% While BARON Academic license can be acquired for a fee,
% CPLEX solver does provide academic free license easily obtainable.
%%
% Clear workspace, command window, and close all figures
clear, clc, close all
% Load data from 'sphingo_net_new3_intra.mat'
load('../model/sphingolipid_network.mat')
% Define parameters
numSamples =  1e6 ;
lambdaOpt = 1e-4;
N = 2;
% i_d = 3;  % = 4 ;
sdFracScale = 1; % 0.5;
tolerance = 1e-3;
robustFit = false;
useParallel = true; % Set to true to use parallel processing, false otherwise
% Input data filename
inputFilename = "../data/N15_dynamic_data.xlsx";
for i_d = [3, 4]  % Choose either 3(C) or 4(D) 
% Determine the dataset identifier for 'v'
if i_d  == 3
    v = 'C';
elseif i_d  == 4
    v = 'D';
end
% Get sheet names from the input Excel file
[version_, date] = version;
if year(date) >= 2020
    fileSheetsNames = sheetnames(inputFilename);
elseif year(date) >= 2011
    fileSheetsNames = xl_xlsfinfo(fullfile(pwd, inputFilename));
else
    fprintf('Code only works for versions of MATLAB released 2011b or after.\n')
end
tmpSheet = fileSheetsNames{i_d};
if year(date) >= 2016
    tmpTbl = readtable(inputFilename, 'Sheet', tmpSheet);
else
    fprintf('Code only works for versions of MATLAB released 2016b or after.\n')
end
StartRunTime_all = tic;
% Set random number generator to default
rng default
% Extract measured metabolites and concentration data
measuredMetabolites = table2cell(tmpTbl(:, 2));
tmpId = find(not(cellfun('isempty', measuredMetabolites)));
measuredMetabolites = measuredMetabolites(tmpId(1:end));
noMetsMeasured = numel(measuredMetabolites);
concentrationData = table2cell(tmpTbl(:, 3:end));
concentrationData_ = concentrationData;
tmpId = find(not(cellfun('isempty', concentrationData_)));
concentrationValues = cell2mat(concentrationData_(tmpId));
minX = min(min(concentrationValues));
maxX = max(max(concentrationValues));
concentrationValues = reshape(normalize(concentrationValues, 'range'), size(concentrationData_)) + 1e-8;
concRawData = concentrationValues;
stdMeasurementRaw = concRawData(:, 8:end);
concentrationData = concRawData(:, 1:7);
stdMeasurement = stdMeasurementRaw;
% Check if Parallel Computing Toolbox is available
if useParallel
    if ~matlab.internal.parallel.isPCTInstalled
        error('Parallel Computing Toolbox is not installed or licensed.');
    end
    numWorkers = 20;
    % Set up parallel pool with all available workers
    parpool('local', numWorkers);
end
% Initialize arrays to store sample results
flux_data_n_all = cell(1, numSamples);
flux_data_sd_n_all = cell(1, numSamples);
flux_data_all = cell(1, numSamples);
flux_data_sd_all = cell(1, numSamples);
EndRunTime_all = zeros(1, numSamples);
nDMFA_all = zeros(1, numSamples);
data_all = cell(1, numSamples);
tn_all = cell(1, numSamples);
fit_ssr_all = zeros(1, numSamples);
% Store reference nDMFA value
sample = 1;
output_sample = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance, []);
reference_nDMFA = output_sample.nDMFA;
% Store results for the reference sample
flux_data_n_all{sample} = output_sample.flux_data_n;
flux_data_sd_n_all{sample} = output_sample.flux_data_sd_n;
flux_data_all{sample} = output_sample.flux_data;
flux_data_sd_all{sample} = output_sample.flux_data_sd;
EndRunTime_all(sample) = output_sample.EndRunTime;
nDMFA_all(sample) = output_sample.nDMFA;
data_all{sample} = output_sample.data;
tn_all{sample} = output_sample.tn;
fit_ssr_all(sample) = output_sample.fit_ssr;
% Valid sample indices
validSampleIndices = [sample];
% Sampling loop
while numel(validSampleIndices) < numSamples
    if useParallel
        % Use parallel processing
        parfor sample = 2:numSamples
            output_sample = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance, []);
            % Check if nDMFA matches the reference value
            if output_sample.nDMFA == reference_nDMFA
                % Store sample results
                flux_data_n_all{sample} = output_sample.flux_data_n;
                flux_data_sd_n_all{sample} = output_sample.flux_data_sd_n;
                flux_data_all{sample} = output_sample.flux_data;
                flux_data_sd_all{sample} = output_sample.flux_data_sd;
                EndRunTime_all(sample) = output_sample.EndRunTime;
                nDMFA_all(sample) = output_sample.nDMFA;
                data_all{sample} = output_sample.data;
                tn_all{sample} = output_sample.tn;
                fit_ssr_all(sample) = output_sample.fit_ssr;

                % Add sample index to validSampleIndices
                validSampleIndices = [validSampleIndices, sample];
            end
            % Print progress
            fprintf('Progress: %.2f%%\n', (sample / numSamples) * 100);
        end
    else
        % Use serial processing
        for sample = 2:numSamples
            output_sample = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance, []);
	        % Check if nDMFA matches the reference value
            if output_sample.nDMFA == reference_nDMFA
                % Store sample results
                flux_data_n_all{sample} = output_sample.flux_data_n;
                flux_data_sd_n_all{sample} = output_sample.flux_data_sd_n;
                flux_data_all{sample} = output_sample.flux_data;
                flux_data_sd_all{sample} = output_sample.flux_data_sd;
                EndRunTime_all(sample) = output_sample.EndRunTime;
                nDMFA_all(sample) = output_sample.nDMFA;
                data_all{sample} = output_sample.data;
                tn_all{sample} = output_sample.tn;
                fit_ssr_all(sample) = output_sample.fit_ssr;
                % Add sample index to validSampleIndices
                validSampleIndices = [validSampleIndices, sample];
            end
            % Print progress
            fprintf('Progress: %.2f%%\n', (sample / numSamples) * 100);
        end
    end
end
% Trim arrays to match the desired number of valid samples
flux_data_n_all = flux_data_n_all(validSampleIndices)';
flux_data_sd_n_all = flux_data_sd_n_all(validSampleIndices)';
flux_data_all = flux_data_all(validSampleIndices)';
flux_data_sd_all = flux_data_sd_all(validSampleIndices)';
EndRunTime_all = EndRunTime_all(validSampleIndices)';
nDMFA_all = nDMFA_all(validSampleIndices)';
data_all = data_all(validSampleIndices)';
tn_all = tn_all(validSampleIndices)';
fit_ssr_all = fit_ssr_all(validSampleIndices)';
% Create a cell array with the data
dataCell = {flux_data_n_all, flux_data_sd_n_all, flux_data_all, flux_data_sd_all, ...
    EndRunTime_all, nDMFA_all, data_all, tn_all, fit_ssr_all};
% Field names for the structure
fieldNames = {'flux_data_n_all', 'flux_data_sd_n_all', 'flux_data_all', 'flux_data_sd_all', ...
    'EndRunTime_all', 'nDMFA_all', 'data_all', 'tn_all', 'fit_ssr_all'};
% Convert cell array to structure array
resultData = cell2struct(dataCell, fieldNames, 2);
% Save the sampling settings to a .mat file
EndRunTime = toc(StartRunTime_all) / 60;
outputfilename = sprintf('../results/%s__SampledDynFluxes_1_%s.mat', v, datestr(now, 'HH-MM-dd-mmm-yyyy'));
save(outputfilename, 'model', 'resultData', 'numSamples', 'lambdaOpt', 'N', 'i_d', 'sdFracScale', 'tolerance', 'robustFit', 'useParallel', 'EndRunTime', '-v7.3');
fprintf('\nData saved to file: %s\n', outputfilename);
% Close the parallel pool (optional)
if useParallel
    delete(gcp);
end


end