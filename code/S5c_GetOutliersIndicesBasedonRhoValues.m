
% STEP 5c
%%
% Clear the workspace, console, and close all figures
clear, clc, close all

filename = 'ConcCorrelations_C_1E2.mat' ;
%%
results_dir_name = '../results' ;
fileName = fullfile(results_dir_name,  filename ) ;
sample = filename(end-8) 
tmp_dataMatrix =  load(fileName, 'concData_all_rhoValues') ;
data = tmp_dataMatrix.concData_all_rhoValues' ;
clear tmp_dataMatrix
% Initialize a logical matrix to mark outliers
outliers = false(size(data)); % Same size as data, initially all false
% Loop through each variable (column) in the dataset
for i = 1:size(data, 2) % Loop through columns
    % Identify outliers in the current variable
    columnOutliers = isoutlier(data(:, i));
    % Mark the outliers in the matrix
    outliers(:, i) = columnOutliers;
end
% To find samples that are outliers in any of the variables
samplesWithOutliers = any(outliers, 2); % True for rows with any outlier
% Identify the indices of samples with outliers
sampleWithOutliersIndices = find(samplesWithOutliers);
sampleWithNoOutliersIndices = find(~samplesWithOutliers);
% Display the indices of samples with outliers
% disp('Total Number of Samples with outliers in any variable:');
% disp(size(sampleWithOutliersIndices, 1));
save(fullfile(results_dir_name, ['ConcCorrOutlierIndices', filename(end-10:end)]), 'sampleWithOutliersIndices', 'sampleWithNoOutliersIndices')
% load('measMets.mat')
% MetsRhoMeanMedianOutliers = [measMets, num2cell(nanmean(data)'), num2cell(nanmedian(data)'), num2cell(sum(outliers)')];
%%
filename = 'ConcCorrelations_D_1E2.mat';
%%
results_dir_name = '../results' ;
fileName = fullfile(results_dir_name,  filename ) ;
sample = filename(end-8) 
tmp_dataMatrix =  load(fileName, 'concData_all_rhoValues') ;
data = tmp_dataMatrix.concData_all_rhoValues' ;
clear tmp_dataMatrix
% Initialize a logical matrix to mark outliers
outliers = false(size(data)); % Same size as data, initially all false
% Loop through each variable (column) in the dataset
for i = 1:size(data, 2) % Loop through columns
    % Identify outliers in the current variable
    columnOutliers = isoutlier(data(:, i));
    % Mark the outliers in the matrix
    outliers(:, i) = columnOutliers;
end
% To find samples that are outliers in any of the variables
samplesWithOutliers = any(outliers, 2); % True for rows with any outlier
% Identify the indices of samples with outliers
sampleWithOutliersIndices = find(samplesWithOutliers);
sampleWithNoOutliersIndices = find(~samplesWithOutliers);
% Display the indices of samples with outliers
% disp('Total Number of Samples with outliers in any variable:');
% disp(size(sampleWithOutliersIndices, 1));
save(fullfile(results_dir_name, ['ConcCorrOutlierIndices', filename(end-10:end)]), 'sampleWithOutliersIndices', 'sampleWithNoOutliersIndices')
% load('measMets.mat')
% MetsRhoMeanMedianOutliers = [measMets, num2cell(nanmean(data)'), num2cell(nanmedian(data)'), num2cell(sum(outliers)')];