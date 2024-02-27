
% STEP 5b
%%
% Clear the workspace, console, and close all figures
clear, clc, close all
results_dir_name = '../results' ;
fileName = fullfile(results_dir_name,  'C_SampledDynFluxes_1_14-15-06-Feb-2024.mat' ) ;
% fileName = fullfile(results_dir_name,  'C_SampledDynFluxes_1_23-18-13-Feb-2024.mat' );
%%
chunksize = 1e2;  % Note this is a random 1E2 sample but actual paper used 1E5;
cultureName = [fileName(12) '_1E2'];  % 1E5 means samples 1E5 flux profiles.
coloridx = 2 ;
figFactor = 3;
colors = {[0.8902    0.4667    0.7608], [1.0000    0.4980    0.0549]} ;
load(fileName,  'numSamples')
% Set the seed for reproducibility
rng('default');%%
% Step 1: Randomly sample 1e4 entries from your cell array%load('C_SampledDynFluxes_1_14-15-06-Feb-2024.mat', 'resultData') ;
indices = randperm(numSamples, chunksize);% 
load(fileName, 'resultData') ;
% Load metadata and relevant data in chunks
measMets = {resultData.data_all{1, 1}.met}';
tmp_data = resultData.flux_data_n_all(indices); 
dataMatrix_fluxN = single(cell2mat(cellfun(@(x) x(:)', tmp_data , 'UniformOutput', false)));
tmp_fileName = fullfile(results_dir_name, ['dataMatrix_fluxN_',cultureName,'.mat'] );
save(tmp_fileName, 'dataMatrix_fluxN', '-v7.3')
clear dataMatrix_fluxN
tmp_data = resultData.res_EstConc_data_all(indices); 
dataMatrix_resid = single(cell2mat(cellfun(@(x) x(:)', tmp_data , 'UniformOutput', false)));
tmp_fileName = fullfile(results_dir_name, [ 'dataMatrix_resid_' cultureName '.mat']) ;
save(tmp_fileName, 'dataMatrix_resid', '-v7.3')
dataCell_mConc = resultData.measConc_data_all(indices); 
dataCell_eConc  = resultData.EstConc_data_all(indices); 
clear resultData
% Spearman correlation
fprintf('\n Computing the Spearman Correlations...')
[concData_all_rhoValues, concData_all_pValues, correlationRatios] = ...
 cellfun(@(x, y) GetSpearmanCorrelations(x, y), dataCell_mConc , dataCell_eConc   , 'UniformOutput', false);
concData_all_rhoValues = cell2mat(concData_all_rhoValues');
concData_all_pValues = cell2mat(concData_all_pValues');
correlationRatios = cell2mat(correlationRatios) ;
tmp_fileName = fullfile(results_dir_name, ['ConcCorrelations_' cultureName '.mat' ]);
save(tmp_fileName, 'concData_all_rhoValues', 'concData_all_pValues', 'correlationRatios', '-v7.3');
clear concData_all_rhoValues concData_all_pValues correlationRatios
%
dataMatrix_eConc= single(cell2mat(cellfun(@(x) x(:)',  dataCell_eConc , 'UniformOutput', false)));
dataMatrix_mConc = single(cell2mat(cellfun(@(x) x(:)', dataCell_mConc , 'UniformOutput', false)));
tmp_fileName = fullfile(results_dir_name, [ 'dataMatrix_emConc_' cultureName '.mat']) ;
save(tmp_fileName, 'dataMatrix_eConc', 'dataMatrix_mConc', '-v7.3')
clear dataMatrix_eConc dataMatrix_mConc dataCell_eConc dataCell_mConc
%
% Define output directory
dir_name = fullfile(results_dir_name, 'CodeFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
fprintf('\n Plotting the Correlations rho values...')
tmp_fileName = fullfile(results_dir_name, ['ConcCorrelations_' cultureName '.mat' ]);
load(tmp_fileName, 'concData_all_rhoValues')
figure(1);
boxplot(concData_all_rhoValues', measMets, 'Orientation','vertical', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median(concData_all_rhoValues, 2), 1:size(median(concData_all_rhoValues, 2), 1), 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold on;
% Get handles to the line graphics objects
h = findobj(gca, 'Tag', 'Box');
for j = 1:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
for j = 2:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
xtickangle(45)
ax= gca ;
ax.YLim = [0 1.2] ;
ax.Color = 'none';
hold off;
fontSize = 16; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0.1 0.1 0.9 0.8]/figFactor; % Set the new figure position (adjust the values as needed)
% Set the figure and axes background to none (transparent)
% fig.Color = 'none';
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
fig_name_plot_png = fullfile(dir_name,  ['SuppFig_4_CorrPlot_' cultureName '.png'] );
saveas(gcf, fig_name_plot_png);
%
clear concData_all_rhoValues
fprintf('\n Plotting the Correlations P-values...')
tmp_fileName = fullfile(results_dir_name, ['ConcCorrelations_' cultureName '.mat' ]);
load(tmp_fileName, 'concData_all_pValues')
figure(2);
boxplot(concData_all_pValues', measMets, 'Orientation','vertical', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median(concData_all_pValues, 2), 1:size(median(concData_all_pValues, 2), 1), 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold on;
% Get handles to the line graphics objects
h = findobj(gca, 'Tag', 'Box');
for j = 1:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
for j = 2:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
xtickangle(45)
fontSize = 16; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0 0 0.9 0.9]; 
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
ax= gca ;
ax.YLim = [0 0.1] ;
ax.Color = 'none';
hold off;
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
fig_name_plot_png = fullfile(dir_name,   ['SuppFig_4_CorrPValuesPlot_' cultureName '.png']);
saveas(gcf, fig_name_plot_png);
%
clear concData_all_pValues
fprintf('\n Reshaping the Residuals for Box Plots...')
tmp_fileName = fullfile(results_dir_name, ['dataMatrix_resid_' cultureName '.mat']) ;
load(tmp_fileName, 'dataMatrix_resid')
if ~isa(dataMatrix_resid, 'single')
    dataMatrix_resid = single(dataMatrix_resid);
end
nMets = 36;
% Step 1: Reshape to a 1nx36x7 matrix
% Step 2: Permute dimensions to move the 3rd dimension to the 2nd
% Step 3: Reshape to a 7nx36 matrix
dataMatrix_resid = reshape(permute(reshape(dataMatrix_resid, size(dataMatrix_resid, 1), nMets, []), [1, 3, 2]), [], nMets);
tmp_fileName = fullfile(results_dir_name, ['dataMatrix_resid4HistPlot_' cultureName '.mat'] );
save(tmp_fileName, 'dataMatrix_resid', '-v7.3')
fprintf('\n Plotting the Residuals Histograms...')
numColumns = size(dataMatrix_resid, 2); % Number of columns in your matrix
numRows = size(dataMatrix_resid, 1); % Number of rows, should be 10000
% Create a figure
figure(3);
% Loop through each column
for i = 1:numColumns
    subplot(6, 6, i); % Adjust the grid size if you want a different layout
    h = histogram(dataMatrix_resid(:, i));
	% Change the color of the bins to a specific color, e.g., a shade of orange
	h.FaceColor = colors{coloridx};
    title(measMets{i}); % Set the title for each subplot
	set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
	ax= gca ;
	ax.Color = 'none';
	hold on
end
% Optional: Improve layout
set(gcf, 'Position', get(0, 'Screensize')); % Make figure full screen
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
hold off;
fig_name_plot_png = fullfile(dir_name,  ['SuppFig_3_ResidualsHistogramPlot_' cultureName '.png']);
saveas(gcf, fig_name_plot_png);
%
fprintf('\n Done with culture: %s !', cultureName)

%%
fileName = fullfile(results_dir_name, 'D_SampledDynFluxes_1_14-16-06-Feb-2024.mat' );
% fileName = fullfile(results_dir_name,  'D_SampledDynFluxes_1_00-26-14-Feb-2024.mat') ;
cultureName = [fileName(12) '_1E2'];  % 1E5 means samples 1E5 flux profiles.
coloridx = 1 ;
load(fileName,  'numSamples')
% Set the seed for reproducibility
rng('default');%%
% Step 1: Randomly sample 1e4 entries from your cell array%load('C_SampledDynFluxes_1_14-15-06-Feb-2024.mat', 'resultData') ;
indices = randperm(numSamples, chunksize);% 
load(fileName, 'resultData') ;
% Load metadata and relevant data in chunks
measMets = {resultData.data_all{1, 1}.met}';
tmp_data = resultData.flux_data_n_all(indices); 
dataMatrix_fluxN = single(cell2mat(cellfun(@(x) x(:)', tmp_data , 'UniformOutput', false)));
tmp_fileName = fullfile(results_dir_name, ['dataMatrix_fluxN_',cultureName,'.mat'] );
save(tmp_fileName, 'dataMatrix_fluxN', '-v7.3')
clear dataMatrix_fluxN
tmp_data = resultData.res_EstConc_data_all(indices); 
dataMatrix_resid = single(cell2mat(cellfun(@(x) x(:)', tmp_data , 'UniformOutput', false)));
tmp_fileName = fullfile(results_dir_name, [ 'dataMatrix_resid_' cultureName '.mat']) ;
save(tmp_fileName, 'dataMatrix_resid', '-v7.3')
dataCell_mConc = resultData.measConc_data_all(indices); 
dataCell_eConc  = resultData.EstConc_data_all(indices); 
clear resultData
% Spearman correlation
fprintf('\n Computing the Spearman Correlations...')
[concData_all_rhoValues, concData_all_pValues, correlationRatios] = ...
 cellfun(@(x, y) GetSpearmanCorrelations(x, y), dataCell_mConc , dataCell_eConc   , 'UniformOutput', false);
concData_all_rhoValues = cell2mat(concData_all_rhoValues');
concData_all_pValues = cell2mat(concData_all_pValues');
correlationRatios = cell2mat(correlationRatios) ;
tmp_fileName = fullfile(results_dir_name, ['ConcCorrelations_' cultureName '.mat' ]);
save(tmp_fileName, 'concData_all_rhoValues', 'concData_all_pValues', 'correlationRatios', '-v7.3');
clear concData_all_rhoValues concData_all_pValues correlationRatios
%
dataMatrix_eConc= single(cell2mat(cellfun(@(x) x(:)',  dataCell_eConc , 'UniformOutput', false)));
dataMatrix_mConc = single(cell2mat(cellfun(@(x) x(:)', dataCell_mConc , 'UniformOutput', false)));
tmp_fileName = fullfile(results_dir_name, [ 'dataMatrix_emConc_' cultureName '.mat']) ;
save(tmp_fileName, 'dataMatrix_eConc', 'dataMatrix_mConc', '-v7.3')
clear dataMatrix_eConc dataMatrix_mConc dataCell_eConc dataCell_mConc
%
% Define output directory
dir_name = fullfile(results_dir_name, 'CodeFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
fprintf('\n Plotting the Correlations rho values...')
tmp_fileName = fullfile(results_dir_name, ['ConcCorrelations_' cultureName '.mat' ]);
load(tmp_fileName, 'concData_all_rhoValues')
figure(4);
boxplot(concData_all_rhoValues', measMets, 'Orientation','vertical', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median(concData_all_rhoValues, 2), 1:size(median(concData_all_rhoValues, 2), 1), 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold on;
% Get handles to the line graphics objects
h = findobj(gca, 'Tag', 'Box');
for j = 1:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
for j = 2:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
xtickangle(45)
ax= gca ;
ax.YLim = [0 1.2] ;
ax.Color = 'none';
hold off;
fontSize = 16; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0.1 0.1 0.9 0.8]/figFactor; % Set the new figure position (adjust the values as needed)
% Set the figure and axes background to none (transparent)
% fig.Color = 'none';
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
fig_name_plot_png = fullfile(dir_name,  ['SuppFig_4_CorrPlot_' cultureName '.png'] );
saveas(gcf, fig_name_plot_png);
%
clear concData_all_rhoValues
fprintf('\n Plotting the Correlations P-values...')
tmp_fileName = fullfile(results_dir_name, ['ConcCorrelations_' cultureName '.mat' ]);
load(tmp_fileName, 'concData_all_pValues')
figure(5);
boxplot(concData_all_pValues', measMets, 'Orientation','vertical', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median(concData_all_pValues, 2), 1:size(median(concData_all_pValues, 2), 1), 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold on;
% Get handles to the line graphics objects
h = findobj(gca, 'Tag', 'Box');
for j = 1:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
for j = 2:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{coloridx}, 'FaceAlpha', 1);
end
xtickangle(45)
fontSize = 16; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0 0 0.9 0.9]; 
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
ax= gca ;
ax.YLim = [0 0.1] ;
ax.Color = 'none';
hold off;
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
fig_name_plot_png = fullfile(dir_name,   ['SuppFig_4_CorrPValuesPlot_' cultureName '.png']);
saveas(gcf, fig_name_plot_png);
%
clear concData_all_pValues
fprintf('\n Reshaping the Residuals for Box Plots...')
tmp_fileName = fullfile(results_dir_name, ['dataMatrix_resid_' cultureName '.mat']) ;
load(tmp_fileName, 'dataMatrix_resid')
if ~isa(dataMatrix_resid, 'single')
    dataMatrix_resid = single(dataMatrix_resid);
end
nMets = 36;
% Step 1: Reshape to a 1nx36x7 matrix
% Step 2: Permute dimensions to move the 3rd dimension to the 2nd
% Step 3: Reshape to a 7nx36 matrix
dataMatrix_resid = reshape(permute(reshape(dataMatrix_resid, size(dataMatrix_resid, 1), nMets, []), [1, 3, 2]), [], nMets);
tmp_fileName = fullfile(results_dir_name, ['dataMatrix_resid4HistPlot_' cultureName '.mat'] );
save(tmp_fileName, 'dataMatrix_resid', '-v7.3')
fprintf('\n Plotting the Residuals Histograms...')
numColumns = size(dataMatrix_resid, 2); % Number of columns in your matrix
numRows = size(dataMatrix_resid, 1); % Number of rows, should be 10000
% Create a figure
figure(6);
% Loop through each column
for i = 1:numColumns
    subplot(6, 6, i); % Adjust the grid size if you want a different layout
    h = histogram(dataMatrix_resid(:, i));
	% Change the color of the bins to a specific color, e.g., a shade of orange
	h.FaceColor = colors{coloridx};
    title(measMets{i}); % Set the title for each subplot
	set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off', 'LineWidth', 1.7);
	ax= gca ;
	ax.Color = 'none';
	hold on
end
% Optional: Improve layout
set(gcf, 'Position', get(0, 'Screensize')); % Make figure full screen
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
hold off;
fig_name_plot_png = fullfile(dir_name,  ['SuppFig_3_ResidualsHistogramPlot_' cultureName '.png']);
saveas(gcf, fig_name_plot_png);

%
fprintf('\n Done with culture: %s !', cultureName)

%%
function [rho, pValue, CorrelationRatios] = GetSpearmanCorrelations(WTProfile, mutantProfile)
% GetSpearmanCorrelations: Computes the Spearman's rank correlation coefficients for each reaction 
% over all time points between the given mutant and WT profiles.
%
% Inputs:
%   WTProfile: Matrix representing the wild-type profile, where each row is a reaction and each column is a time point.
%   mutantProfile: Matrix representing the mutant profile, with the same structure as WTProfile.
%
% Outputs:
%   rho: Vector containing the Spearman's rank correlation coefficients for each reaction.
%   pValue: Vector containing the p-values corresponding to each correlation coefficient.
%   positiveCorrelationRatio: Proportion of reactions that have rho > 0.5 and p < 0.05.
%   negativeCorrelationRatio: Proportion of reactions that have rho < -0.5 and p < 0.05.
% Determine the number of reactions from the mutant profile.
numReactions = size(mutantProfile, 1);
rho = zeros(numReactions, 1);
pValue = zeros(numReactions, 1);
% Loop through each reaction to compute the correlation coefficient.
for r = 1:numReactions
    [rho(r), pValue(r)] = corr(mutantProfile(r, :)', WTProfile(r, :)', 'Type', 'Spearman');
end
% Calculate the proportion of reactions showing positive correlation (rho > 0.5) with statistical significance (p < 0.05).
positivelyCorrelated = (rho > 0.5) & (pValue < 0.05);
positiveCorrelationRatio = sum(positivelyCorrelated) / numReactions;
% Calculate the proportion of reactions showing negative correlation (rho < -0.5) with statistical significance (p < 0.05).
negativelyCorrelated = (rho < -0.5) & (pValue < 0.05);
negativeCorrelationRatio = sum(negativelyCorrelated) / numReactions;
CorrelationRatios = [positiveCorrelationRatio, negativeCorrelationRatio] ;
end