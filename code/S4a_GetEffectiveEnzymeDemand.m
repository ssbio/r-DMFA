
%%
% Clear the workspace, console, and close all figures
clear, clc, close all
% Load the network data from a MAT file
load('../model/sphingolipid_network.mat')
% Load the kcat distribution from mat file
load('../results/kcat_dist.mat')
HighConfidenceData_idx = model.HighConfidenceData_idx ;
met_list_X_er = model.met_list_X_er;
tmp_data1 = C_k_data_Perhr;
tmp_data2 = D_k_data_Perhr ;
tmp_data1 = tmp_data1(HighConfidenceData_idx, :)  ;
tmp_data2 = tmp_data2(HighConfidenceData_idx, :) ;
met_list_X_ = met_list_X_er(HighConfidenceData_idx, 2)  ;
met_list_X_enzClass = met_list_X_er(HighConfidenceData_idx, 1)  ;
[met_list_X_enzClass_unique,~,idx_ec] = unique(met_list_X_enzClass) ;
tmp_rxnIdxs = [] ;
% Loop through each element in cellArray2
for i = 1:numel(met_list_X_)
    element = met_list_X_{i};
    tmp_rxnIdxs(i,1) = find(cellfun(@(x) isequal(x, element), model.rxnIDs)) ;
end
load('../results/C_DynFit_WT.mat','measuredMetabolites', 'conc_pred_bounds', 'meas_est_x_vals_lb', 'meas_est_x_vals_ub', 'idxsDays', 'resultData')
tmp_fluxNorm = resultData.flux_data_n_all{1, 1}(tmp_rxnIdxs, 1:end-1) ;
tmp_fluxNorm_sd = resultData.flux_data_sd_n_all{1, 1}(tmp_rxnIdxs, 1:end-1) ;
tmp_fluxNorm_lb = tmp_fluxNorm - tmp_fluxNorm_sd ;
tmp_fluxNorm_ub = tmp_fluxNorm + tmp_fluxNorm_sd ;
% dimensions
nvars = size(tmp_fluxNorm, 1);
nSamples = size(tmp_data1, 2);
nTimepoints = size(tmp_fluxNorm, 2);
% Initialize an empty matrix C
Econc = [] ;
% Loop through each time point in C
for t = 1:nTimepoints
    % Perform element-wise division B/A for the current time point
    tmp_Econc_a = tmp_fluxNorm_lb(:, t) ./ tmp_data1;
    tmp_Econc_b = tmp_fluxNorm_ub(:, t) ./ tmp_data1;
    Econc = [Econc tmp_Econc_a tmp_Econc_b] ;
end
Econc_1 = Econc ;
ECONC{1} = Econc ;
load('../results/D_DynFit_WT.mat','measuredMetabolites', 'conc_pred_bounds', 'meas_est_x_vals_lb', 'meas_est_x_vals_ub', 'idxsDays', 'resultData')
tmp_fluxNorm = resultData.flux_data_n_all{1, 1}(tmp_rxnIdxs, 1:end) ;
tmp_fluxNorm_sd = resultData.flux_data_sd_n_all{1, 1}(tmp_rxnIdxs, 1:end) ;
tmp_fluxNorm_lb = tmp_fluxNorm - tmp_fluxNorm_sd ;
tmp_fluxNorm_ub = tmp_fluxNorm + tmp_fluxNorm_sd ;
% dimensions
nvars = size(tmp_fluxNorm, 1);
nSamples = size(tmp_data1, 2);
nTimepoints = size(tmp_fluxNorm, 2);
% Initialize an empty matrix C
Econc = [] ;
% Loop through each time point in D
for t = 1:nTimepoints
    % Perform element-wise division B/A for the current time point
    tmp_Econc_a = tmp_fluxNorm_lb(:, t) ./ tmp_data1;
    tmp_Econc_b = tmp_fluxNorm_ub(:, t) ./ tmp_data1;
    Econc = [Econc tmp_Econc_a tmp_Econc_b] ;
end
Econc_2 = Econc ;
ECONC{2} = Econc ;
p_values_e = [];
p_values_e_ = [];
h_values_e = [];
met_list_X2_e = {};
k_data2_e = [];
rho_values_e = [] ;
chi2_values_e = [];
for iii = 1:numel(met_list_X_enzClass_unique)
        met_list_X2_e = [met_list_X2_e ; [met_list_X_enzClass_unique{iii, 1} '_C']; [met_list_X_enzClass_unique{iii, 1} '_D']] ;
        tmp_x =  trimOutliers(sum(Econc_1(idx_ec==iii,:), 1)) ;
        tmp_y =  trimOutliers(sum(Econc_2(idx_ec==iii,:), 1))  ;
        pd_x = fitdist(tmp_x, 'Normal');
        pd_y = fitdist(tmp_y, 'Normal');
        nSamp = 5e2;
        rng('default')  % For reproducibility
        tmp_x = random(pd_x, nSamp,1);
        rng('default')  % For reproducibility
        tmp_y = random(pd_y,nSamp,1);
        overlap(iii,1)= bhattacharyyCoef(pd_x , pd_y) ;
        tmp_x(tmp_x<=0) = NaN;
        tmp_y(tmp_y<=0) = NaN;
        k_data2_e = [k_data2_e; tmp_x'; tmp_y'] ;
        median_data_e(iii, :) = [median(tmp_x), median(tmp_y)] ;
      [h_values_e(iii, 1), p_values_e(iii, 1)] = kstest2(tmp_x, tmp_y , 'Alpha', 0.01) ;
end
overlap_mean = mean(overlap)
overlap_std = std(overlap)
%
% Define output directory
dir_name = fullfile('../results', 'CodeFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
% Create y-axis labels for reactions
y_labels2 = met_list_X2_e';
% Calculate the median of each row (along the second dimension, dim=2)
median_of_each_row = median(k_data2_e, 2);
figure (1);
boxplot(k_data2_e', y_labels2, 'Orientation','vertical', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median_of_each_row, 1:size(median_of_each_row, 1), 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold off;
ylabel('Apparent Enzyme Costs');
xlabel('Enzyme IDs');
title('Distribution of Relative Enzyme Costs');
% Get handles to the line graphics objects
h = findobj(gca, 'Tag', 'Box');
colors = { [173/255, 216/255, 230/255], [25/255, 25/255, 112/255]};
for j = 1:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{1}, 'FaceAlpha', 1);
end
for j = 2:2:length(h)
    patch(get(h(j), 'XData'), get(h(j), 'YData'), colors{2}, 'FaceAlpha', 1);
end
% Set custom x-tick locations
desired_ticks = 1.5:2:length(y_labels2) ;
xticks(desired_ticks) ;
xticklabels(met_list_X_enzClass_unique);
xtickangle(45)
% Define custom x-positions for vertical grid lines
grid_x_positions = 2.5:2:length(y_labels2);
% Add vertical lines at the specified positions
for x_val = grid_x_positions
    line([x_val x_val], get(gca, 'YLim'), 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
end
hold off;  % Release the plot
fontSize = 18; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'LineWidth', 2);
dummyC = patch(NaN, NaN, colors{2}, 'FaceAlpha', 1);
dummyD = patch(NaN, NaN, colors{1}, 'FaceAlpha', 1);
legend([dummyC, dummyD], {'Bio. Rep.- 1', 'Bio. Rep.- 2'}, 'Box', 'off');
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0.1 0.1 0.9 0.8]; % Set the new figure position (adjust the values as needed)
fig_name_plot_png = fullfile(dir_name, 'Fig_5a_EffEnzCostDist.png');
saveas(gcf, fig_name_plot_png);
%
function trimmed_data = trimOutliers(data)
    % This function trims outliers from the data using the IQR method.
    %
    % Input:
    %   - data: A column vector of data points
    %
    % Output:
    %   - trimmed_data: Data with outliers removed based on the IQR method
    % Ensure data is a column vector
    if size(data, 2) > 1
        data = data';
    end
    % Compute Q1, Q3, and IQR
    Q1 = quantile(data, 0.25);
    Q3 = quantile(data, 0.75);
    IQR = Q3 - Q1;
    % Compute lower and upper bounds
    LB = Q1 - 1.5 * IQR;
    UB = Q3  + 1.5 * IQR;
    % Filter data to keep values within bounds
    trimmed_data = data(data >= LB & data <= UB);
end
function bc = bhattacharyyCoef(pdf1, pdf2)
    pdf1 = makedist('Normal','mu',pdf1.mu,'sigma',pdf1.sigma);
    pdf2 = makedist('Normal','mu',pdf2.mu,'sigma',pdf2.sigma);
    fun = @(x) sqrt(pdf(pdf1,x) .* pdf(pdf2,x)) ;
    bc = integral(fun,0,inf);
end