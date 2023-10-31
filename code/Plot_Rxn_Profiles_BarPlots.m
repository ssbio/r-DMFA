

figNo = 1;
EPSImage = false;
load('color_bar_id.mat')

% Extract relevant data
rxnIDs = model.rxnIDs;
simdata = output_sample.simdata;
nDMFA = output_sample.nDMFA;
K = output_sample.K;
c0 = output_sample.c0;
U = output_sample.U;
tn = output_sample.tn;
tm = output_sample.tm;
time_interp = output_sample.tm;
time = output_sample.time;
data = output_sample.data;
fit = output_sample.fit;
fit_ssr = output_sample.fit_ssr;
SEL_ = output_sample.selection;
figNo = output_sample.figNo;

% Define variables for concentration data
conc_data_ = output_sample.meas_raw_x_vals;
meas_est_x_vals = output_sample.meas_est_x_vals;
meas_est_x_vals_sd = output_sample.meas_est_x_vals_sd;
idxsDays = output_sample.indicesDays;

% Set 'v' based on 'i_d'
switch i_d
    case 3
        v = 'C';
    case 4
        v = 'D';
end

% Create a directory for saving plots
dir_name = [pwd '\Metab_fit\' v];
if ~(exist(dir_name,'dir') == 7)
    mkdir(dir_name);
end

% Calculation of bounds
meas_est_x_vals_ub = meas_est_x_vals + meas_est_x_vals_sd;
meas_est_x_vals_lb = meas_est_x_vals - meas_est_x_vals_sd;
meas_est_x_vals_lb(meas_est_x_vals_lb < 0) = 1e-6;
err_sd = meas_est_x_vals_sd;
c_errlow = meas_est_x_vals_lb;
c_errhigh = meas_est_x_vals_ub;

% Initialize variables for metabolite names and mean data
met_names = [];
Mean_data = [];

% Plot concentration profiles
for j = 1:numel(met_list_X)
    if EPSImage
        % Higher resolution and editable in Adobe Illustrator images (eps)
        figure(figNo + 1)
        er = errorbar(time_interp(idxsDays), meas_est_x_vals(j, idxsDays), err_sd(j, idxsDays), 'o-b', 'linewidth', 2.0);
        hold on
        xline(2, '--g', 'linewidth', 3.0);
        xline(5, '--g', 'linewidth', 3.0);
        hold off
        xlabel('Time [days]', 'FontSize', 20)
        ylabel('Measurements [nmol/gdw]', 'FontSize', 20)
        set(gca, 'box', 'off')
        set(gca, 'FontName', 'Times', 'FontSize', 18)
        title(char(met_list_X(j)))
        set(gcf, 'color', 'none');
        set(gca, 'color', 'none');
        exportgraphics(gcf, 'transparent.eps', 'ContentType', 'vector', 'BackgroundColor', 'none')
        fig_name_plot = [dir_name '\' char(met_list_X(j)) '.eps'];
        exportgraphics(gcf, fig_name_plot, 'ContentType', 'vector', 'BackgroundColor', 'none');
    end
    
    % Non-editable format
    figure(figNo + 1)
    plot(time_interp, meas_est_x_vals(j, :), '-b', 'linewidth', 2.0)
    hold on
    scatter(time, conc_data_(j, :), 'r*')
    patch([time_interp fliplr(time_interp)], [meas_est_x_vals_ub(j, :) fliplr(meas_est_x_vals_lb(j, :))], 'g', 'FaceAlpha', 0.05)
    xline(2, '--g', 'linewidth', 3.0)
    xline(5, '--r', 'linewidth', 3.0)
    hold off
    legend({'predicted', 'rawdata', 'CI_{2.5 - 97.5%}'})
    legend('Box', 'off')
    grid on
    xlabel('Time [days]')
    ylabel('Concentration [nmol/gDW]')
    title([char(met_list_X(j))])
    fontSize = 12;
    set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
    fig_name_plot = [dir_name '\' char(met_list_X(j)) '.png'];
    saveas(gcf, fig_name_plot);
    
    met_names = [met_names; met_list_X(j)];
end

figNo = figNo + 1;

% TODO: Handle saving Mean_DATA if needed

% Create folders for saving flux plots
folderPaths = {
    fullfile(pwd, 'MainRxnFluxes', 'NormFluxStats');
    fullfile(pwd, 'MainRxnFluxes')
};

for i = 1:numel(folderPaths)
    if exist(folderPaths{i}, 'dir') ~= 7
        mkdir(folderPaths{i});
    end
end

% Calculate time-dependent parameters gamma and kappa
if length(tn) > 1  % DMFA analysis
    [g_t, k_t] = dynpars(tm, tn);
else  % MFA analysis
    g_t = tm - tm(1);
    k_t = ones(size(tm));
end

Nr = size(rxnIDs, 1);  % 78 reactions

% Normalize rates to the first reaction in the network
tmpCellRxnFluxes = fit(nDMFA).V * k_t;
tmpCellRxnFluxes_sd = fit(nDMFA).V_sd * k_t;

tmpCellRxnFluxes_ = tmpCellRxnFluxes;
tmpCellRxnFluxes_sd_ = tmpCellRxnFluxes_sd;

% Compute normalized fluxes
V_pred = tmpCellRxnFluxes_(1:Nr, :);
V_predNorm = V_pred ./ (V_pred(1, :));
V_pred_sd = tmpCellRxnFluxes_sd_(1:Nr, :);
V_pred_sd(abs(V_pred_sd) >= abs(V_pred)) = 0.25 * abs(V_pred(abs(V_pred_sd) >= abs(V_pred)));
V_predNorm_sd = V_pred_sd;

rxn_names = [];
Mean_data_rxn = [];

% Plot reaction fluxes
for j = 1:numel(rxnIDs)
    tmp_labels = rxnIDs(j);
    [~, tExp_idx_1] = min(abs(time_interp - 2));
    [~, tExp_idx_2] = min(abs(time_interp - 5));
    time_exphase = time_interp(tExp_idx_1:tExp_idx_2);
    tmp_mean_1 = mean(tmpCellRxnFluxes_(j, 3:6));
    tmp_mean_sd = mean(tmpCellRxnFluxes_sd_(j, 3:6));
    
    % Create bar graphs
    figure(figNo + 2);
    bar(time, tmpCellRxnFluxes_(j, :), 'FaceColor', color_id_rxn2{j})
    hold on
    er = errorbar(time, tmpCellRxnFluxes_(j, :), tmpCellRxnFluxes_sd_(j, :), 'o', 'color', 'black', 'linewidth', 1.5);
    er.LineStyle = 'none';
    er.CapSize = 2;
    hold off
    title(char(rxnIDs(j)));
    xlim([0, 7])
    xticks([1, 2, 5, 6])
    xticklabels({'1', '2', '5', '6'})
    xlabel('Time [days]')
    ylabel('Flux [mmol/(gDW*h)]')
    set(gca, 'box', 'off')
    fontSize = 10;
    set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
    
    if j == 1
        legend(tmp_labels, 'Location', 'northeast')
        legend('boxoff')
    end
    
    fig_name_plot = [dir_name '\MainRxnFluxes\' char(rxnIDs(j)) '.png'];
    saveas(gcf, fig_name_plot);
    rxn_names = [rxn_names; rxnIDs(j)];
    Mean_data_rxn = [Mean_data_rxn; tmp_mean_1 tmp_mean_sd];
end

figNo = figNo + 1;

% Prepare reaction flux data for saving to a file
Mean_DATA_rxn = [cellstr(rxn_names), num2cell(Mean_data_rxn)];
tmp_Mean_DATA_rxn = [cellstr('MetNames'), cellstr('PredictedMean'), cellstr('StandardDev'); Mean_DATA_rxn];

% Define the filename for saving
% tmp_filename_rxn = "RXNSMeanExpPhaseValues_1.xlsx";
% fprintf('\n Writing Analysis to %s file...\n', tmp_filename_rxn);
% if year(date) < 2019
%     xlswrite(tmp_filename_rxn, tmp_Mean_DATA_rxn, v);
% else
%     writecell(tmp_Mean_DATA_rxn, tmp_filename_rxn, 'Sheet', v, 'UseExcel', false);
% end

%% Saving to file

% Extract relevant flux data for saving
flux_data_n = V_predNorm(:, idxsDays);
flux_data_sd_n = V_predNorm_sd(:, idxsDays);

flux_data = V_pred(:, idxsDays);
flux_data_sd = V_pred_sd(:, idxsDays);

% flux_data = V_predNorm(:, idxsDays(3:5));
% flux_data_sd = V_predNorm_sd(:, idxsDays(3:5));
% save('dyn_ss_fluxes_C.mat', 'flux_data', 'flux_data_sd', 'flux_data_n', 'flux_data_sd_n')

fprintf('\n Done.  \n ');

% EndRunTime = toc(StartRunTime)/60
MatFileName = [v '_DMFA_fit_N_', char(string(N)), '_', datestr(now, 'HH-MM-dd-mmm-yyyy'), '.mat'];
% save(MatFileName)

%% Calculate time-dependent parameters gamma, kappa, and lambda
function [g, k, l, dg, dk, dg2, dk2] = dynpars(t, tn)
    % DYNPARS Calculate time-dependent parameters gamma, kappa, and lambda.
    % These parameters are used for dynamic metabolic flux analysis (DMFA).
    %
    %   [G, K, L, DG, DK, DG2, DK2] = DYNPARS(T, TN) computes time-dependent
    %   parameters gamma (G), kappa (K), lambda (L), and the first and second
    %   order derivatives of G and K with respect to TN, at the requested time(s) T. 
    %   TN are the DMFA time points. T should be within the time domain
    %   TN(1)...TN(end), which includes N-2 inflection points:
    %   TN = (t_begin, t_infection_1, t_infection_2, ..., t_end). 

    % Number of requested time points (m), and number of DMFA time points (n)
    m = length(t);   % number of requested time points
    n = length(tn);  % number of DMFA time points

    % Initialize matrices for parameters
    g = zeros(n, m);
    k = zeros(n, m);
    l = zeros(n, m);
    dg = zeros(n, m, n);  % first order derivatives
    dk = zeros(n, m, n);
    dg2 = zeros(n, m, n);  % second order derivatives
    dk2 = zeros(n, m, n);

    % Compute parameter values
    for j = 1:m  % loop over requested time points
        for i = 1:n  % loop over DMFA time points
            [g(i, j), k(i, j), l(i, j), dg(i, j, :), dk(i, j, :), dg2(i, j, :), dk2(i, j, :)] = ...
                calcpars(t(j), tn, i);
        end
    end
end

% Calculate time-dependent parameters gamma, kappa, and lambda for a specific DMFA time point
function [g, k, l, dg, dk, dg2, dk2] = calcpars(t, tn, i)
    % Calculate time-dependent parameters g, k, l and also first and second
    % order derivatives of g and k with respect to DMFA time points (tn).
    % This function calculates the values corresponding to the ith DMFA time point
    % at the requested time point t.  

    % The total number of DMFA time points
    n = length(tn);

    % Initialize
    g = 0;
    k = 0;
    l = 0;
    dg = zeros(1, n);
    dk = zeros(1, n);
    dg2 = zeros(1, n);
    dk2 = zeros(1, n);

    % Calculate parameter values for time point t
    if i > 1 && t >= tn(1) && t < tn(i-1)  % 1st row of table 1
        g = 0;
        k = 0;
        l = 0;
        dg = zeros(1, n);
        dk = zeros(1, n);
        dg2 = zeros(1, n);
        dk2 = zeros(1, n);

    % ... (rest of the code remains the same)
    
    else  % time point outside time domain
        errmsg = 'Requested time point is outside time domain.';
        disp(sprintf('%s: %s', mfilename, errmsg)); % display error message
    end
end
