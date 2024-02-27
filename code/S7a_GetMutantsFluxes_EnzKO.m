% STEP 7a
%%
%%
% Clear the workspace, console, and close all figures
clear, clc, close all
results_dir_name = '../results' ;
% Load the network data from a MAT file
load('../model/sphingolipid_network.mat')

%
lambdaOpt = 1e-4;
N = 2;
% i_d = 3;
sdFracScale = 0; %0.5;
tolerance = 1e-3;
robustFit = false;
% Set the flag for parallel processing
useParallel = false; % Set to true to use parallel processing, false otherwise
uniqueEnz = unique(model.EnzLabel);
indexCellArray = cell(length(uniqueEnz), 1);
for i = 1:length(uniqueEnz)
    indexCellArray{i} = find(strcmp(model.EnzLabel, uniqueEnz{i}))';
end
all_combinations = indexCellArray;
%
numSamples = length(all_combinations); %5; % Set the number of samples you want to perform

for i_d = [3, 4] % Choose either 3(C) or 4(D) 
% Determine the dataset identifier for 'v'
if i_d == 3
    v = 'C';
elseif i_d == 4
    v = 'D';
end
%
% Specify the input Excel file
inputFilename = '../data/N15_dynamic_data.xlsx';
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
% Set the random number generator to its default state
rng default
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
%
if i_d==4
    stdMeasurement = stdMeasurement(:, 1:6);
    concentrationData = concentrationData(:, 1:6) ;
end
%
StartRunTime_all = tic ;
% Preallocate arrays to store the values of the selected fields for all samples
flux_data_n_all = cell(1, numSamples);
flux_data_sd_n_all = cell(1, numSamples);
flux_data_all = cell(1, numSamples);
flux_data_sd_all = cell(1, numSamples);
EndRunTime_all = zeros(1, numSamples);
nDMFA_all = zeros(1, numSamples);
data_all = cell(1, numSamples);
tn_all = cell(1, numSamples);
fit_ssr_all = zeros(1, numSamples);
accum_rate_all = cell(1, numSamples) ;
SampleName_all = cell(1, numSamples) ;
% Keep track of valid sample indices where nDMFA matches the reference value
validSampleIndices = [];
sample = 1 ;
% Store reference value for nDMFA
% output_sample = DMFASampler(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance);
output_sample = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance, []);
% keep desired reference value
reference_nDMFA = output_sample.nDMFA; 
% Store the values of the selected fields for the current sample
flux_data_n_all{sample} = output_sample.flux_data_n;
flux_data_sd_n_all{sample} = output_sample.flux_data_sd_n;
flux_data_all{sample} = output_sample.flux_data;
flux_data_sd_all{sample} = output_sample.flux_data_sd;
EndRunTime_all(sample) = output_sample.EndRunTime;
nDMFA_all(sample) = output_sample.nDMFA;
data_all{sample} = output_sample.data;
tn_all{sample} = output_sample.tn;
fit_ssr_all(sample) = output_sample.fit_ssr;
accum_rate_all{sample} = model.S*output_sample.flux_data ; 
SampleName_all{sample} = 'WT' ;
% Add the sample index to the validSampleIndices array
validSampleIndices = [validSampleIndices, sample];
% Loop until the desired number of valid samples is reached
% while numel(validSampleIndices) < numSamples
    % Use parallel processing (parfor) or serial processing (for loop) based on the flag
    if useParallel
        % Use parallel processing
        parfor sample = 2:numSamples
            blockedRxnsIdx = all_combinations{sample} ;
            output_sample = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance, blockedRxnsIdx);
           % Check if nDMFA matches the reference value
                % Store the values of the selected fields for the current sample
                flux_data_n_all{sample} = output_sample.flux_data_n;
                flux_data_sd_n_all{sample} = output_sample.flux_data_sd_n;
                flux_data_all{sample} = output_sample.flux_data;
                flux_data_sd_all{sample} = output_sample.flux_data_sd;
                EndRunTime_all(sample) = output_sample.EndRunTime;
                nDMFA_all(sample) = output_sample.nDMFA;
                data_all{sample} = output_sample.data;
                tn_all{sample} = output_sample.tn;
                fit_ssr_all(sample) = output_sample.fit_ssr;
                accum_rate_all{sample} = model.S*output_sample.flux_data.*(output_sample.flux_data_n(1)) ; 
                SampleName_all{sample} = strjoin(model.rxnIDs(blockedRxnsIdx), ",") ;
                % Add the sample index to the validSampleIndices array
                validSampleIndices = [validSampleIndices, sample];
            % Print progress to screen
            fprintf('Progress: %.2f%%\n', (sample / numSamples) * 100);
        end
    else
        % Use serial processing
        for sample = 2:numSamples
            blockedRxnsIdx = all_combinations{sample} ;
            output_sample = DMFASampler_constr(model, measuredMetabolites, concentrationData, stdMeasurement, minX, maxX, N, lambdaOpt, sdFracScale, robustFit, tolerance, blockedRxnsIdx);
            % Check if nDMFA matches the reference value
                % Store the values of the selected fields for the current sample
                flux_data_n_all{sample} = output_sample.flux_data_n;
                flux_data_sd_n_all{sample} = output_sample.flux_data_sd_n;
                flux_data_all{sample} = output_sample.flux_data;
                flux_data_sd_all{sample} = output_sample.flux_data_sd;
                EndRunTime_all(sample) = output_sample.EndRunTime;
                nDMFA_all(sample) = output_sample.nDMFA;
                data_all{sample} = output_sample.data;
                tn_all{sample} = output_sample.tn;
                fit_ssr_all(sample) = output_sample.fit_ssr;
                accum_rate_all{sample} = model.S*output_sample.flux_data.*(output_sample.flux_data_n(1)) ; 
                SampleName_all{sample} = strjoin(model.rxnIDs(blockedRxnsIdx), ",") ;
                % Add the sample index to the validSampleIndices array
                validSampleIndices = [validSampleIndices, sample];
            % Print progress to screen
            fprintf('Progress: %.2f%%\n', (sample / numSamples) * 100);
        end
    end
% end
% Trim the arrays to match the desired number of valid samples
flux_data_n_all = flux_data_n_all(validSampleIndices)';
flux_data_sd_n_all = flux_data_sd_n_all(validSampleIndices)';
flux_data_all = flux_data_all(validSampleIndices)';
flux_data_sd_all = flux_data_sd_all(validSampleIndices)';
EndRunTime_all = EndRunTime_all(validSampleIndices)';
nDMFA_all = nDMFA_all(validSampleIndices)';
data_all = data_all(validSampleIndices)';
tn_all = tn_all(validSampleIndices)';
fit_ssr_all = fit_ssr_all(validSampleIndices)';
accum_rate_all = accum_rate_all(validSampleIndices)';
SampleName_all = SampleName_all(validSampleIndices)' ;
% Create a cell array with the data
dataCell = {flux_data_n_all, flux_data_sd_n_all, flux_data_all, flux_data_sd_all, ...
            EndRunTime_all, nDMFA_all, data_all, tn_all, fit_ssr_all, accum_rate_all,  SampleName_all, validSampleIndices'};
% Specify the field names for the structure
fieldNames = {'flux_data_n_all', 'flux_data_sd_n_all', 'flux_data_all', 'flux_data_sd_all', ...
              'EndRunTime_all', 'nDMFA_all', 'data_all', 'tn_all', 'fit_ssr_all', 'accum_rate_all', 'SampleName_all', 'ValidSampleIndices' };
% Convert the cell array to a structure array
resultData = cell2struct(dataCell, fieldNames, 2);
% Save the sampling settings
EndRunTime = toc(StartRunTime_all) / 60;
fprintf('\nDone with culture: %s \n', v);
%%
% Generate a filename and save the sampling settings to a .mat file
outputfilename = sprintf('%s_DynMOMAsingleEnzs.mat', v);
save(fullfile(results_dir_name, outputfilename), 'model', 'resultData', 'numSamples', 'lambdaOpt', 'N', 'i_d', 'sdFracScale', 'tolerance', 'robustFit', 'useParallel' , 'EndRunTime', 'v', '-v7.3');
fprintf('\nData saved to file: %s\n', outputfilename);

end


