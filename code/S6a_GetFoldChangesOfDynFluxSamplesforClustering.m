
% STEP 6a
%%
% Clear the workspace, console, and close all figures
clear, clc, close all

% filename = 'dataMatrix_fluxN_C_1E5.mat' ;
% %filename = 'dataMatrix_fluxN_D_1E5.mat' ;
% filename2 = 'ConcCorrOutlierIndicess_C_1E5.mat' ;
% %filename2 = 'ConcCorrOutlierIndicess_D_1E5.mat'


% filename = 'dataMatrix_fluxN_C_1E6.mat' ;
% %filename = 'dataMatrix_fluxN_D_1E6.mat' ;
% filename2 = 'ConcCorrOutlierIndicess_C_1E6.mat' ;
% %filename2 = 'ConcCorrOutlierIndicess_D_1E6.mat'

filename = 'dataMatrix_fluxN_C_1E2.mat' ;
% filename = 'dataMatrix_fluxN_D_1E2.mat' ;

filename2 = 'ConcCorrOutlierIndicess_C_1E2.mat' ;
% filename2 = 'ConcCorrOutlierIndicess_D_1E2.mat' ;
%%
results_dir_name = '../results' ;
fileName = fullfile(results_dir_name,  filename ) ;
sample = filename(end-8) ;
cultureName = [sample '_1E2'];  % 1E5 means samples 1E5 flux profiles.
% Step 2: Concatenate rows of each 40x7 double array into a single vector
% of 1x280 per sample... 
tmp_dataMatrix = load(fileName, 'dataMatrix_fluxN');
dataMatrix = tmp_dataMatrix.dataMatrix_fluxN ;
clear tmp_dataMatrix
load(fullfile(results_dir_name,[sample '_SampledDynFluxes_ref_1.mat'])) ;% C_SampledDynFluxes_ref_1.mat;
% Load the network data from a MAT file
load('../model/sphingolipid_network.mat')
% Selecting the columns corresponding to timepoints T2-T5
selectedTimepoints = 81:240;
dataMatrix = dataMatrix(:, selectedTimepoints);
dataMatrix_ref_1 = dataMatrix_ref(:, selectedTimepoints);
% Parameters for reshaping the data
nRxns = 40;
nSel_timepoints = 4;
nReplicates = 50;
%
% Define output directory
dir_name = fullfile(results_dir_name, 'CodesFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
% Reshaping the data matrix
reshapedData = cellfun(@(row) reshape(row, [nRxns, nSel_timepoints]), num2cell(dataMatrix,2), 'UniformOutput', false);
reshapedData_ref = reshape(dataMatrix_ref_1(:, :), [nRxns, nSel_timepoints]);
refArray = reshapedData_ref;
% Preallocating result matrix for DTW distances
fcMatrixCell = cell(size(dataMatrix, 1),1);
% Computing the foldchanges for reshaped data
for i = 1:numel(reshapedData)
    tempArray = reshapedData{i};
    % Compute fold change between each row in tempArray and the corresponding row in refArray
    tmp_log2fcArray = sign(tempArray(:,:)./refArray(:,:)) .* log2(abs(tempArray(:,:)./refArray(:,:)));
    fcMatrixCell{i} =  tmp_log2fcArray ;
end
% Sample data for demonstration
n_reactions = nRxns;
n_timepoints  = nSel_timepoints;
N = length(fcMatrixCell); % Number of samples per flux per epoch
% Create x-axis labels for epochs
x_labels = arrayfun(@(x) strcat('T_', num2str(x)), 2:2+n_timepoints-1, 'UniformOutput', false);
% Initialize a matrix to store the combined data
heatmap_data = zeros(nRxns, n_timepoints * N);
% Iterate through the cell array and reshape the data
for timepoint = 1:n_timepoints
    for sn = 1:N
        heatmap_data(:, (timepoint-1)*N + sn) = fcMatrixCell{sn}(:, timepoint);
    end
end
% Create y-axis labels for reactions
y_labels = model.rxnIDs';%
% Calculate the median of each row (along the second dimension, dim=2)
median_of_each_row = median(heatmap_data, 2);
figure;
boxplot(heatmap_data', y_labels, 'Orientation','horizontal', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median_of_each_row, 1:n_reactions, 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold off;
grid on
xlabel('log_2(FC)');
ylabel('Reaction IDs');
% title('Distribution of fluxes log_2(FC)');
fontSize = 12; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0.1 0.1 0.9 0.8]; % Set the new figure position (adjust the values as needed)
fig_name_plot = fullfile(dir_name, ['Fig_6_FC_boxplot_' cultureName  '.png']);
saveas(gcf, fig_name_plot);
%
tmp_filename = '_File.csv';
writematrix(heatmap_data', fullfile(results_dir_name, [sample '_FCdata_' tmp_filename]))
writecell(model.EnzLabel, fullfile(results_dir_name, [sample '_Enzymelabel_' tmp_filename]))
writecell(model.rxnIDs, fullfile(results_dir_name, [sample '_Rxnlabel_' tmp_filename]))
% Iterate through the cell array and reshape the data
for timepoint = 1:n_timepoints
    for sn = 1:N
        dataByRxns(:, (timepoint-1)*N + sn) = reshapedData{sn}(:, timepoint);
    end
end
writematrix(dataByRxns', fullfile(results_dir_name, fullfile(results_dir_name,[sample '_data_' tmp_filename])) ) 
fprintf('written to file! \n')
%
% Generate a filename and save the sampling settings to a .mat file
outputfilename = sprintf( [sample '_heatmapdata_%s.mat', datestr(now, 'HH-MM-dd-mmm-yyyy')]);
% save(fullfile(results_dir_name,outputfilename), '-v7.3');
% 
fprintf('\nData saved to file: %s\n', outputfilename);


%% WITHOUT OUTLIERS (NO OUTLIERS)
%%
fprintf('\nGenerating FC files with no outliers...');

fileName2 = fullfile(results_dir_name,  filename2 ) ;
load(fileName2)
%%
sample = filename2(end-8) ;
cultureName = [sample '_1E2'];  % 1E5 means samples 1E5 flux profiles.
% Step 2: Concatenate rows of each 40x7 double array into a single vector
% of 1x280 per sample... 
tmp_dataMatrix = load(fileName, 'dataMatrix_fluxN');
dataMatrix = tmp_dataMatrix.dataMatrix_fluxN ;
clear tmp_dataMatrix
% Selecting the columns corresponding to timepoints T2-T5
selectedTimepoints = 81:240;
dataMatrix = dataMatrix(sampleWithNoOutliersIndices, selectedTimepoints);
dataMatrix_ref_1 = dataMatrix_ref(:, selectedTimepoints);
% Parameters for reshaping the data
nRxns = 40;
nSel_timepoints = 4;
nReplicates = 50;
%
% Define output directory
dir_name = fullfile(results_dir_name, 'CodesFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
% Reshaping the data matrix
reshapedData = cellfun(@(row) reshape(row, [nRxns, nSel_timepoints]), num2cell(dataMatrix,2), 'UniformOutput', false);
reshapedData_ref = reshape(dataMatrix_ref_1(:, :), [nRxns, nSel_timepoints]);
refArray = reshapedData_ref;
% Preallocating result matrix for DTW distances
fcMatrixCell = cell(size(dataMatrix, 1),1);
% Computing the foldchanges for reshaped data
for i = 1:numel(reshapedData)
    tempArray = reshapedData{i};
    % Compute fold change between each row in tempArray and the corresponding row in refArray
    tmp_log2fcArray = sign(tempArray(:,:)./refArray(:,:)) .* log2(abs(tempArray(:,:)./refArray(:,:)));
    fcMatrixCell{i} =  tmp_log2fcArray ;
end
% Sample data for demonstration
n_reactions = nRxns;
n_timepoints  = nSel_timepoints;
N = length(fcMatrixCell); % Number of samples per flux per epoch
% Create x-axis labels for epochs
x_labels = arrayfun(@(x) strcat('T_', num2str(x)), 2:2+n_timepoints-1, 'UniformOutput', false);
% Initialize a matrix to store the combined data
heatmap_data = zeros(nRxns, n_timepoints * N);
% Iterate through the cell array and reshape the data
for timepoint = 1:n_timepoints
    for sn = 1:N
        heatmap_data(:, (timepoint-1)*N + sn) = fcMatrixCell{sn}(:, timepoint);
    end
end
% Create y-axis labels for reactions
y_labels = model.rxnIDs';%
% Calculate the median of each row (along the second dimension, dim=2)
median_of_each_row = median(heatmap_data, 2);
figure;
boxplot(heatmap_data', y_labels, 'Orientation','horizontal', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median_of_each_row, 1:n_reactions, 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold off;
grid on
xlabel('log_2(FC)');
ylabel('Reaction IDs');
% title('Distribution of fluxes log_2(FC)');
fontSize = 12; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0.1 0.1 0.9 0.8]; % Set the new figure position (adjust the values as needed)
fig_name_plot = fullfile(dir_name, ['Fig_6_FC_boxplot_NoOutlier_' cultureName  '.png']);
saveas(gcf, fig_name_plot);
%
tmp_filename = '_FileNoOutlier.csv';
writematrix(heatmap_data', fullfile(results_dir_name, [sample '_FCdata_' tmp_filename]))
writecell(model.EnzLabel, fullfile(results_dir_name, [sample '_Enzymelabel_' tmp_filename]))
writecell(model.rxnIDs, fullfile(results_dir_name, [sample '_Rxnlabel_' tmp_filename]))
% Iterate through the cell array and reshape the data
for timepoint = 1:n_timepoints
    for sn = 1:N
        dataByRxns(:, (timepoint-1)*N + sn) = reshapedData{sn}(:, timepoint);
    end
end
writematrix(dataByRxns', fullfile(results_dir_name, fullfile(results_dir_name,[sample '_data_' tmp_filename])) ) 
fprintf('written to file! \n')
%
% Generate a filename and save the sampling settings to a .mat file
outputfilename = sprintf( [sample '_heatmapdata_NoOutlier_%s.mat', datestr(now, 'HH-MM-dd-mmm-yyyy')]);
% save(fullfile(results_dir_name,outputfilename), '-v7.3');
% 
fprintf('\nData saved to file: %s\n', outputfilename);
%%
%% ONLY OUTLIERS (OUTLIERS ONLY)
%%
fprintf('\nGenerating FC files for ONLY outliers...');
%%
sample = filename2(end-8) ;
cultureName = [sample '_1E2'];  % 1E5 means samples 1E5 flux profiles.
% Step 2: Concatenate rows of each 40x7 double array into a single vector
% of 1x280 per sample... 
tmp_dataMatrix = load(fileName, 'dataMatrix_fluxN');
dataMatrix = tmp_dataMatrix.dataMatrix_fluxN ;
clear tmp_dataMatrix
% Selecting the columns corresponding to timepoints T2-T5
selectedTimepoints = 81:240;
dataMatrix = dataMatrix(sampleWithOutliersIndices, selectedTimepoints);
dataMatrix_ref_1 = dataMatrix_ref(:, selectedTimepoints);
% Parameters for reshaping the data
nRxns = 40;
nSel_timepoints = 4;
nReplicates = 50;
%
% Define output directory
dir_name = fullfile(results_dir_name, 'CodesFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
% Reshaping the data matrix
reshapedData = cellfun(@(row) reshape(row, [nRxns, nSel_timepoints]), num2cell(dataMatrix,2), 'UniformOutput', false);
reshapedData_ref = reshape(dataMatrix_ref_1(:, :), [nRxns, nSel_timepoints]);
refArray = reshapedData_ref;
% Preallocating result matrix for DTW distances
fcMatrixCell = cell(size(dataMatrix, 1),1);
% Computing the foldchanges for reshaped data
for i = 1:numel(reshapedData)
    tempArray = reshapedData{i};
    % Compute fold change between each row in tempArray and the corresponding row in refArray
    tmp_log2fcArray = sign(tempArray(:,:)./refArray(:,:)) .* log2(abs(tempArray(:,:)./refArray(:,:)));
    fcMatrixCell{i} =  tmp_log2fcArray ;
end
% Sample data for demonstration
n_reactions = nRxns;
n_timepoints  = nSel_timepoints;
N = length(fcMatrixCell); % Number of samples per flux per epoch
% Create x-axis labels for epochs
x_labels = arrayfun(@(x) strcat('T_', num2str(x)), 2:2+n_timepoints-1, 'UniformOutput', false);
% Initialize a matrix to store the combined data
heatmap_data = zeros(nRxns, n_timepoints * N);
% Iterate through the cell array and reshape the data
for timepoint = 1:n_timepoints
    for sn = 1:N
        heatmap_data(:, (timepoint-1)*N + sn) = fcMatrixCell{sn}(:, timepoint);
    end
end
% Create y-axis labels for reactions
y_labels = model.rxnIDs';%
% Calculate the median of each row (along the second dimension, dim=2)
median_of_each_row = median(heatmap_data, 2);
figure;
boxplot(heatmap_data', y_labels, 'Orientation','horizontal', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median_of_each_row, 1:n_reactions, 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold off;
grid on
xlabel('log_2(FC)');
ylabel('Reaction IDs');
% title('Distribution of fluxes log_2(FC)');
fontSize = 12; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0.1 0.1 0.9 0.8]; % Set the new figure position (adjust the values as needed)
fig_name_plot = fullfile(dir_name, ['Fig_6_FC_boxplot_Outliers_' cultureName  '.png']);
saveas(gcf, fig_name_plot);
%
tmp_filename = '_FileOutliers.csv';
writematrix(heatmap_data', fullfile(results_dir_name, [sample '_FCdata_' tmp_filename]))
writecell(model.EnzLabel, fullfile(results_dir_name, [sample '_Enzymelabel_' tmp_filename]))
writecell(model.rxnIDs, fullfile(results_dir_name, [sample '_Rxnlabel_' tmp_filename]))
% Iterate through the cell array and reshape the data
for timepoint = 1:n_timepoints
    for sn = 1:N
        dataByRxns(:, (timepoint-1)*N + sn) = reshapedData{sn}(:, timepoint);
    end
end
writematrix(dataByRxns', fullfile(results_dir_name, fullfile(results_dir_name,[sample '_data_' tmp_filename])) ) 
fprintf('written to file! \n')
%
% Generate a filename and save the sampling settings to a .mat file
outputfilename = sprintf( [sample '_heatmapdata_Outliers_%s.mat', datestr(now, 'HH-MM-dd-mmm-yyyy')]);
% save(fullfile(results_dir_name,outputfilename), '-v7.3');
% 
fprintf('\nData saved to file: %s\n', outputfilename);