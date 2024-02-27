
% STEP 7c
%%
%%
% Clear the workspace, console, and close all figures
clear, clc, close all
results_dir_name = '../results' ;
% Load the network data from a MAT file
load('../model/sphingolipid_network.mat')
%%
sample = 'C' ;
load(fullfile(results_dir_name, [sample '_DynMOMAsingleEnzs.mat']), 'resultData')
[resultData_C, zd_log_norm_C, mutant_labels_C]  = processReactionData(resultData, model);
sample = 'D' ;
load(fullfile(results_dir_name, [sample '_DynMOMAsingleEnzs.mat']), 'resultData')
[resultData_D, zd_log_norm_D, mutant_labels_D]  = processReactionData(resultData, model);
% Define output directory
dir_name = fullfile(results_dir_name, 'CodeFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
trunc_id = 0 ; start_id = 2;
% Sample data (replace with actual values)
posRatio1 = cell2mat(resultData_C.flux_data_n_all_posRatio); 
negRatio1 = cell2mat(resultData_C.flux_data_n_all_negRatio); 
posRatio2 = cell2mat(resultData_D.flux_data_n_all_posRatio); 
negRatio2 = cell2mat(resultData_D.flux_data_n_all_negRatio); 
reactionNames = mutant_labels_C;
% Combine the data for stacked bar plot
allData1 = [negRatio1, posRatio1]; allData2 = [negRatio2, posRatio2];
reactionNames = reactionNames(start_id:end-trunc_id) ;
allData1 = allData1(start_id:end-trunc_id, :) ; allData2 = allData2(start_id:end-trunc_id, :) ;
% Create the stacked horizontal bar plots % Create a larger figure for clarity
figure('Position', [300, 110, 1200, 200]);
bgap = 3; % Increase this if you need more space between bars of different reactions
% Define y-positions for the bars
y_idx1 = 1:bgap:(bgap*length(reactionNames));
y_idx2 = y_idx1 + 1; % Offset for the second culture
% Plot for the first culture with custom bar width
bar(y_idx1, allData1, 0.4, 'stacked', 'EdgeColor', 'none'); 
hold on;
% Plot for the second culture with custom bar width and offset
bar(y_idx2, allData2, 0.4, 'stacked', 'EdgeColor', 'none');
% Adjust the y-axis labels to be the reaction names
xticks((y_idx1 + y_idx2) / 2); % Center the y-ticks between the bars of the two cultures
xticklabels(reactionNames);
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
% fig.Color = 'none';
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
ax= gca ;
% ax.Color = 'none';
fontSize = 16; 
set(gca,'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off');
fig_name_plot_png = fullfile(dir_name,  'Fig_7b_EnzKO_RxnCorrPosNeg.png');
saveas(gcf, fig_name_plot_png);
% Plot Z index
trunc_id = 0 ; reactionNames = mutant_labels_C;
% Combine the data for stacked bar plot
allData1 = zd_log_norm_C; allData2 = zd_log_norm_D;
reactionNames = reactionNames(start_id:end-trunc_id) ;
allData1 = allData1(start_id:end-trunc_id, :) ;
allData2 = allData2(start_id:end-trunc_id, :) ;
% Create the stacked horizontal bar plots
figure('Position', [300, 110, 1200, 200]);
bgap = 3; % Increase this if you need more space between bars of different reactions
% Define y-positions for the bars
y_idx1 = 1:bgap:(bgap*length(reactionNames));
y_idx2 = y_idx1 + 1; % Offset for the second culture
% Plot for the first culture with custom bar width
bar(y_idx1, allData1, 0.4, 'stacked', 'EdgeColor', 'none'); 
hold on;
% Plot for the second culture with custom bar width and offset
bar(y_idx2, allData2, 0.4, 'stacked', 'EdgeColor', 'none');
% Adjust the y-axis labels to be the reaction names
xticks((y_idx1 + y_idx2) / 2); % Center the y-ticks between the bars of the two cultures
xticklabels(reactionNames);
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
% fig.Color = 'none';
fig.InvertHardcopy = 'off'; % This ensures the background remains transparent
ax= gca ;
% ax.Color = 'none';
fontSize = 16; 
set(gca,'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Helvetica', 'Box', 'off');
fig_name_plot_png = fullfile(dir_name,  'Fig_7a_EnzKO_Zindex.png');
saveas(gcf, fig_name_plot_png);

%%
function cvArray = computeCV(mean_data, std_data)
    cvArray= abs(std_data./mean_data);
    cvArray(isnan(cvArray)) = 1;
    cvArray(cvArray<eps) = 1;
end
function normalized_data = normalizeData(data)
    % Handle -inf and +inf
    data(data == -inf) = 0; % Replace -inf with the smallest finite value
    data(data == inf) = max(data(isfinite(data)))+std(data(isfinite(data)), 'omitnan');  % Replace inf with the largest finite value
    % Normalize to [0, 1]
    normalized_data = (data - min(data)) / (max(data) - min(data));
end
function [rho, pValue, positiveCorrelationRatio, negativeCorrelationRatio] = GetSpearmanCorrelations(WTProfile, mutantProfile)
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
end
function [resultData, zd_log_norm, mutant_labels] = processReactionData(resultData, model)
    % Extract necessary data and labels from resultData structure
    mutant_labels = unique(model.EnzLabel);
    mutant_labels{1} = 'WT'; 
    mutant_labels = mutant_labels(resultData.ValidSampleIndices);
    resultData.flux_data_n_all_CV = cellfun(@(x, y) computeCV(x, y), resultData.flux_data_n_all, resultData.flux_data_sd_n_all, 'UniformOutput', false);
    dataCell_mut_nFlux = cellfun(@(x) x(:, 3:6), resultData.flux_data_n_all, 'UniformOutput', false);
    dataCell_mut_nFlux_CV = cellfun(@(x) x(:, 3:6), resultData.flux_data_n_all_CV, 'UniformOutput', false);
    dataMatrix_mut_nFlux = cell2mat(cellfun(@(x) x(:)', dataCell_mut_nFlux, 'UniformOutput', false));
    dataMatrix_mut_nFlux_CV = cell2mat(cellfun(@(x) x(:)', dataCell_mut_nFlux_CV, 'UniformOutput', false));
    % Z index calculation
    WT_flux = dataMatrix_mut_nFlux(1, :);
    V_VWT_VWT = abs((dataMatrix_mut_nFlux - WT_flux) ./ WT_flux);
    z_components = V_VWT_VWT ./ dataMatrix_mut_nFlux_CV;
    z_d = sum(z_components, 2, 'omitnan');
    zd_log = log(z_d);
    zd_log_norm = normalizeData(zd_log);
    % Spearman correlation
    fixedWTData = resultData.flux_data_n_all{1};
    [resultData.flux_data_n_all_rhoValues, resultData.flux_data_n_all_pValues, ...
     resultData.flux_data_n_all_posRatio, resultData.flux_data_n_all_negRatio] = ...
     cellfun(@(y) GetSpearmanCorrelations(fixedWTData, y), resultData.flux_data_n_all, 'UniformOutput', false);
    % Additional plotting or operations can be added here if needed
end
