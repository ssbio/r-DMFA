clear, clc
close all
%%
% Get the enrichment profiles, The ratios
Get_EnrichmentProfiles
%%
% sample = 'C';
samples = ['C'; 'D' ];
for sample_i = 1:numel(samples)
    sample=samples(sample_i);

%%
if strcmp(sample, 'C')
    kk_vec = [9:10] ;
elseif strcmp(sample, 'D')
    kk_vec = [11:12] ;
end
nBootstraps = 1e4;
%%
N15_mean_stacked = [];
N15_std_stacked = [];
for j = 1:no_mets_measured  % Replace num_j_values with the actual number of j values you have
    tmp_keep = [];
    for kk = kk_vec
        tmp_keep = [tmp_keep; R_S2{kk}(j, :)];
    end
    N15_tmp_keep_mean = mean(tmp_keep, 1);
    N15_mean_stacked = [N15_mean_stacked; N15_tmp_keep_mean];
    N15_tmp_keep_std = std(tmp_keep, [], 1) ;
    N15_std_stacked = [N15_std_stacked; N15_tmp_keep_std];
end
%%
time_vec_days = time_interp;
time_vec_hrs = time_interp*24;
time_vec_mins = time_interp*24*60;
time_vec_secs = time_interp*24*60*60;
%%
N15_mean_stacked = N15_mean_stacked(:, selectedDays(1:end-1)); 
N15_std_stacked = N15_std_stacked(:, selectedDays(1:end-1)); 
time_vec_days = time_interp(1, selectedDays(1:end-1));
time_vec_hrs = time_interp(1, selectedDays(1:end-1))*24;
time_vec_mins = time_interp(1, selectedDays(1:end-1))*24*60;
time_vec_secs = time_interp(1, selectedDays(1:end-1))*24*60*60;
%%
standard_errors = N15_std_stacked ;
[k_values_Persecs_b, k_stddev_Persecs_b, k_bootstraps_Persecs_b]  = estimate_k_with_bootstrap(time_vec_secs, N15_mean_stacked, standard_errors, nBootstraps);
[k_values_Permin_b, k_stddev_Permin_b, k_bootstraps_Permin_b]  =  estimate_k_with_bootstrap(time_vec_mins, N15_mean_stacked,  standard_errors, nBootstraps) ;
[k_values_Perhr_b, k_stddev_Perhr_b, k_bootstraps_Perhr_b]  =  estimate_k_with_bootstrap(time_vec_hrs, N15_mean_stacked,  standard_errors, nBootstraps) ;
[k_values_Perday_b, k_stddev_Perday_b, k_bootstraps_Perday_b]  =  estimate_k_with_bootstrap(time_vec_days, N15_mean_stacked,  standard_errors, nBootstraps); 
disp('Estimated k values for each metabolite using Bootstraping');
k_values_table_b = table(met_list_X, k_values_Persecs_b, k_stddev_Persecs_b, k_values_Permin_b, k_stddev_Permin_b,...
            k_values_Perhr_b, k_stddev_Perhr_b, k_values_Perday_b, k_stddev_Perday_b, ...
            'VariableNames', {'MetIDs', 'k (1/s)', 'std (1/s)', 'k (1/min)', 'std (1/min)', 'k (1/hr)', 'std (1/hr)', 'k (1/day)', 'std (1/day)'}) ;
% k_data = k_bootstraps_Persecs_b;
% k_data = k_bootstraps_Permin_b;
% k_data = k_bootstraps_Perhr_b;
k_data = k_bootstraps_Perday_b;
% Create y-axis labels for reactions
y_labels = met_list_X';
% Calculate the median of each row (along the second dimension, dim=2)
median_of_each_row = median(k_data, 2);

% figure;
% boxplot(k_data', y_labels, 'Orientation','horizontal', 'Symbol','.')
% hold on;
% % Add a scatter plot of the median values
% scatter(median_of_each_row, 1:no_mets_measured, 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
% hold off;
% grid on
% xlabel('k_{cat} values');
% ylabel('Metabolites');
% title('Distribution of k_{cat} values');
% fontSize = 12; 
% set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
% % Automatically resize the figure before saving
% fig = gcf; % Get the current figure handle
% fig.Units = 'normalized'; % Set the units to normalized
% % fig.Position = [0.1 0.1 0.8 0.6]; % Set the new figure position (adjust the values as needed)
% fig.Position = [0.1 0.1 0.9 0.8]; % Set the new figure position (adjust the values as needed)
%%
filename_save = '../results/kcat_dist.mat';
varName = [sample, '_k_data_Perhr'];
eval([varName, ' = k_data;']);
if exist(filename_save, 'file') == 2
    disp('File exists.');
    save(filename_save, varName, '-append')
else
    disp('File does not exist.');
    save(filename_save, varName)
end
end
%%
load('../model/sphingolipid_network.mat')
%
HighConfidenceData_idx = model.HighConfidenceData_idx ;
met_list_X_er = model.met_list_X_er;
tmp_data1 = C_k_data_Perhr;
tmp_data2 = D_k_data_Perhr ;
tmp_data1 = tmp_data1(HighConfidenceData_idx, :)  ;
tmp_data2 = tmp_data2(HighConfidenceData_idx, :) ;
% met_list_X_ = met_list_X(HighConfidenceData_idx)  ;
met_list_X_ = met_list_X_er(HighConfidenceData_idx, 2)  ;
%%
%
R = 8.314e-3;  % kJ/(mol K)
T = 295.15;  % Kelvin, assuming 22°C
% Example concentration matrix (rows are compounds, first column is min concentration, second column is max concentration)
conc_matrix = [100 1000; 1.7696 1.7696];  % Fill in with your data
load('../results/C_DynFit_WT.mat','measuredMetabolites', 'conc_pred_bounds', 'meas_est_x_vals_lb', 'meas_est_x_vals_ub', 'idxsDays')
conc_pred_bounds = [min(meas_est_x_vals_lb(:, idxsDays(3:end)), [], 2)  max(meas_est_x_vals_ub(:, idxsDays(3:end)), [], 2)] ;
conc_matrix = [conc_matrix; conc_pred_bounds] ;
conc_matrix(conc_matrix<0.1) = 0.1 ;
gamma(:, 1) = GetThermoGamma(measuredMetabolites, conc_matrix, model, R, T) ;
%
conc_matrix = [100 1000; 1.7696 1.7696];  % Fill in with your data
load('../results/D_DynFit_WT.mat','measuredMetabolites', 'conc_pred_bounds', 'meas_est_x_vals_lb', 'meas_est_x_vals_ub', 'idxsDays')
conc_pred_bounds = [min(meas_est_x_vals_lb(:, idxsDays(3:end)), [], 2)  max(meas_est_x_vals_ub(:, idxsDays(3:end)), [], 2)] ;
conc_matrix = [conc_matrix; conc_pred_bounds] ;
conc_matrix(conc_matrix<0.1) = 0.1 ; 
gamma(:, 2) = GetThermoGamma(measuredMetabolites, conc_matrix, model, R, T) ;
%% Rough plot of thermodynamic saturation
figure;
plot(gamma(:, 1), 'bo-', 'LineWidth', 2)
hold on
plot(gamma(:, 2), 'r+--', 'LineWidth', 2)
hold off
xticks(1:length(gamma(:,1)))
xticklabels(model.rxnIDs(1:length(gamma(:,1))))
grid on
ylabel('\gamma values');
xlabel('Reaction IDs');
title('Distribution of {\gamma} values');
fontSize = 12; 
set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'LineWidth', 2);
% Automatically resize the figure before saving
fig = gcf; % Get the current figure handle
fig.Units = 'normalized'; % Set the units to normalized
fig.Position = [0.1 0.1 0.9 0.8]; % Set the new figure position (adjust the values as needed)
legend('Biological Replicate - 1', 'Biological Replicate - 2', 'Box', 'off');

%%
p_values = [];
h_values = [];
met_list_X2 = {};
k_data2 = [];
for i = 1:size(met_list_X_, 1)
    if ~all(isnan(tmp_data1(i,:))) && ~all(isnan(tmp_data2(i,:)))
        met_list_X2 = [met_list_X2 ; [met_list_X_{i, 1} '_C']; [met_list_X_{i, 1} '_D']] ;
        tmp_idx = find(cellfun(@(x) strcmp(x, met_list_X_{i}), model.rxnIDs), 1) ;
        tmp_x =  trimOutliers(tmp_data1(i,:))./gamma(tmp_idx, 1) ;
        tmp_y =  trimOutliers(tmp_data2(i,:))./gamma(tmp_idx, 2) ;
        pd_x = fitdist(tmp_x, 'Normal') ;
        pd_y = fitdist(tmp_y, 'Normal') ;
        overlap(i,1)= bhattacharyyCoef(pd_x , pd_y) ;
        nSamp = 5e2;
        rng('default')  % For reproducibility
        tmp_x = random(pd_x, nSamp,1);
        rng('default')  % For reproducibility
        tmp_y = random(pd_y,nSamp,1);
        tmp_x(tmp_x<=0) = NaN;
        tmp_y(tmp_y<=0) = NaN;
       k_data2 = [k_data2; tmp_x'; tmp_y'] ;
       median_data(i, :) = [median(tmp_x), median(tmp_y)] ;
      [h_values(i, 1), p_values(i, 1)] = kstest2(tmp_x, tmp_y , 'Alpha', 0.01) ;
    else 
%         i
    end
end
% Display which variables have distributions that appear to be different
overlap_mean = mean(overlap)
overlap_std = std(overlap)
%%
% Define output directory
dir_name = fullfile('../results/', 'CodeFigures');
if ~isfolder(dir_name)
    mkdir(dir_name);
end
%%
% Create y-axis labels for reactions
y_labels2 = met_list_X2';
% Calculate the median of each row (along the second dimension, dim=2)
median_of_each_row = median(k_data2, 2);
figure;
% boxplot(k_data2', y_labels2, 'Orientation','horizontal', 'Symbol','.')
boxplot(k_data2', y_labels2, 'Orientation','vertical', 'Symbol','.')
hold on;
% Add a scatter plot of the median values
scatter(median_of_each_row, 1:size(median_of_each_row, 1), 5, 'r', 'filled', 'MarkerEdgeColor', 'k');
hold off;
% grid on
ylabel('Apparent k_{cat^{+}}  values');
xlabel('Reaction IDs');
title('Apparent k_{cat^{+}} distribution');
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
xticklabels(met_list_X_);
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
fig_name_plot_png = fullfile(dir_name,  'Fig_3_kcatDist.png');
saveas(gcf, fig_name_plot_png);
%%
%%
%%
%%
%%
%% Function to estimate the kcat/gamma values
function k_values = estimate_k_for_all_mets(timepoints, incorp_ratios)
    % Number of metabolites
    nMets = size(incorp_ratios, 1);
    % Initialize a vector to store k values for all metabolites
    k_values = zeros(nMets, 1);
    % Iterate over each metabolite and estimate k
    for i = 1:nMets
        linear_form = log(1 - incorp_ratios(i, :));
        fit_params = polyfit(timepoints, linear_form, 1);
        k_values(i) = -fit_params(1);
    end
    % Optional: plot the results
    for i = 1:nMets
        figure(2);
        linear_form = log(1 - incorp_ratios(i, :));
        plot(timepoints, linear_form, 'o', timepoints, polyval([-k_values(i), 0], timepoints), '-');
        xlabel('Time');
        ylabel('ln[1 - Incorporation Ratio]');
        title(['Metabolite ', num2str(i), ' with k = ', num2str(k_values(i))]);
        grid on;
        legend({'data', 'fit'})
    end
end
%%
function [k_values, k_stddev, k_bootstraps] = estimate_k_with_bootstrap(timepoints, incorp_ratios, standard_errors, nBootstraps)
    % Number of metabolites
    nMets = size(incorp_ratios, 1);
    % Initialize matrices to store bootstrapped k values
    k_bootstraps = zeros(nMets, nBootstraps);
    rng default % for reproducibility
    % Perform bootstrapping
    for b = 1:nBootstraps
        % Generate a bootstrapped dataset
        bootstrap_ratios = incorp_ratios + standard_errors .* randn(size(incorp_ratios));
        
        % Estimate k for each metabolite using bootstrapped dataset
        for i = 1:nMets
            linear_form = log(1 - bootstrap_ratios(i, :));
            fit_params = polyfit(timepoints, linear_form, 1);
            k_bootstraps(i, b) = -real(fit_params(1));
        end
    end
    % Calculate the mean and standard deviation of the bootstrapped k values
    k_values = mean(k_bootstraps, 2);
    k_stddev = std(k_bootstraps, 0, 2);
end
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
    UB = Q3 + 1.5 * IQR;
    % Filter data to keep values within bounds
    trimmed_data = data(data >= LB & data <= UB);
end
% Get thermodynamic factor
function gamma = GetThermoGamma(measuredMetabolites, conc_matrix, model, R, T)
% Example stoichiometric matrix (rows are compounds, columns are reactions)
stoich_matrix = model.S;  % Fill in with your data
deltaG_0 = model.SGFE;  % Fill in with your standard free energy changes for each reaction
num_reactions = size(stoich_matrix, 2)-3;
deltaG_prime = zeros(num_reactions, 1);
deltaG_prime_save = zeros(num_reactions, 1);
for col = 3:num_reactions
    reaction = stoich_matrix(:, col);
    reactants_idx = find(reaction < 0);  % Indices where reactants are
    products_idx = find(reaction > 0);   % Indices where products are
    reactants_idx_ = strcmp(measuredMetabolites, model.rawMetabName(reactants_idx));   % index of substrate metabolite
    products_idx_ = strcmp(measuredMetabolites, model.rawMetabName(products_idx));   % index of product metabolite
    % Get minimum concentrations of reactants and maximum concentrations of products
    reactant_min_concs = conc_matrix(reactants_idx_, 1);
    product_max_concs = conc_matrix(products_idx_, 2);
    % Calculate the ratio for the ln term
    ratio = prod(product_max_concs) / prod(reactant_min_concs);
    % Compute ΔG' for the reaction
    deltaG = deltaG_0(col) + R * T * log(ratio);
    deltaG_prime_save(col) = deltaG;
    if deltaG>0 
        deltaG = -1;
    end
    deltaG_prime(col) = deltaG;
end
gamma = 1-exp(deltaG_prime./( R * T)) ;
gamma(1:2) = 1;
end

function bc = bhattacharyyCoef(pdf1, pdf2)
    pdf1 = makedist('Normal','mu',pdf1.mu,'sigma',pdf1.sigma);
    pdf2 = makedist('Normal','mu',pdf2.mu,'sigma',pdf2.sigma);
    fun = @(x) sqrt(pdf(pdf1,x) .* pdf(pdf2,x)) ;
    bc = integral(fun,0,inf);
end
