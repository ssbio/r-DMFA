% Initialize figure number and EPSImage flag
figNo = 1;
EPSImage = false;

% Load color bar ID data
load('color_bar_id.mat')

% Extract relevant data from the output_sample structure
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

% Define a variable 'v' based on the value of 'i_d'
switch i_d
    case 3
        v = 'C';
    case 4
        v = 'D';
end

% Create a directory for saving plots if it doesn't exist
dir_name = [pwd '\Metab_fit\' v];
if ~(exist(dir_name,'dir') == 7)
    mkdir(dir_name);
end

% Calculate upper and lower bounds for measured data
meas_est_x_vals_ub = meas_est_x_vals + meas_est_x_vals_sd;
meas_est_x_vals_lb = meas_est_x_vals - meas_est_x_vals_sd;
meas_est_x_vals_lb(meas_est_x_vals_lb < 0) = 1e-6;

% Initialize variables for metabolite names and mean data
met_names = [];
Mean_data = [];

% Plot concentration profiles for each metabolite
for j = 1:numel(met_list_X)
    if EPSImage
        % Create higher resolution and editable EPS images
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

    % Create non-editable format PNG images
    figure(figNo + 1)
    plot(time_interp,  meas_est_x_vals(j, :), '-b', 'linewidth', 2.0)
    hold on
    scatter(time, conc_data_(j, :), 'r*')
    patch([time_interp fliplr(time_interp)], [meas_est_x_vals_ub(j, :) fliplr(meas_est_x_vals_lb(j, :))], 'g', 'FaceAlpha', 0.05)
    xline(2,'--g', 'linewidth', 3.0);
    xline(5,'--r', 'linewidth', 3.0);
    hold off
    legend({'predicted', 'rawdata', 'CI_{2.5 - 97.5%}'});
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

% Save metabolite mean data to a file (code commented out)

% Define folder paths for saving main reaction flux plots
folderPaths = {
    fullfile(pwd, 'MainRxnFluxes', 'NormFluxStats');
    fullfile(pwd, 'MainRxnFluxes')
};

% Create the necessary folders if they don't exist
for i = 1:numel(folderPaths)
    if exist(folderPaths{i}, 'dir') ~= 7
        mkdir(folderPaths{i});
    end
end

% Calculate time-dependent parameters gamma and kappa
if length(tn) > 1
    [g_t, k_t] = dynpars(tm, tn);
else
    g_t = tm - tm(1);
    k_t = ones(size(tm));
end

Nr = size(rxnIDs, 1);

% Normalize rates to the first reaction in the network
tmpCellRxnFluxes = fit(nDMFA).V * k_t;
tmpCellRxnFluxes_sd = fit(nDMFA).V_sd * k_t;

tmpCellRxnFluxes_ = tmpCellRxnFluxes;
tmpCellRxnFluxes_sd_ = tmpCellRxnFluxes_sd;

% Define UptakePool condition (code commented out)

% Extract normalized reaction fluxes
V_pred = tmpCellRxnFluxes_(1:Nr, :);
V_predNorm = V_pred ./ (V_pred(1, :));
V_pred_sd = tmpCellRxnFluxes_sd_(1:Nr, :);
V_pred_sd(abs(V_pred_sd) >= abs(V_pred)) = 0.25 * abs(V_pred(abs(V_pred_sd) >= abs(V_pred)));
V_predNorm_sd = V_pred_sd;

% Define reaction names
rxn_names = [];

% Plot and save main reaction fluxes
for j = 1:numel(rxnIDs)
    tmp_labels = rxnIDs(j);
    [~, tExp_idx_1] = min(abs(time_interp - 2));
    [~, tExp_idx_2] = min(abs(time_interp - 5));
    time_exphase = time_interp(tExp_idx_1:tExp_idx_2);
    tmp_mean_1 = mean(tmp_rxn_data(:, 3:6));
    tmp_mean_sd = mean(tmp_rxn_data_sd(:, 3:6));
    ss = spline(time, tmp_rxn_data, time_interp);
    ss_sd = spline(time, tmp_rxn_data_sd, time_interp);
    
    % Create and save plots
    figure(figNo + 1);
    plot(time_interp, ss, '-b', 'linewidth', 2.0);
    hold on
    patch([time_interp fliplr(time_interp)], [ss + ss_sd fliplr(ss - ss_sd)], 'g', 'FaceAlpha', 0.05);
    xline(2, '--g', 'linewidth', 3.0);
    xline(5, '--r', 'linewidth', 3.0);
    yline(0, '--k', 'linewidth', 3.0);
    hold off
    xlabel('Time [days]');
    ylabel('Normalized Flux_{|SPT_{rxn}|}');
    grid on
    title([char(string(strrep(tmp_labels, '_', '-'))) '  (V_{mean} = ' char(num2str(round(tmp_mean_1, 5))) ')']);
    fontSize = 12; 
    set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
    fig_name_plot = fullfile(pwd, 'MainRxnFluxes', 'NormFluxStats', v, [char(rxnIDs(j)) '.png']);
    folder_path = fileparts(fig_name_plot);
    if exist(folder_path, 'dir') ~= 7
        mkdir(folder_path);
    end
    saveas(gcf, fig_name_plot);
    
    % Create and save bar graphs
    figure(figNo + 2);
    bar(time, tmp_rxn_data, 'FaceColor', color_id_rxn{j});
    hold on
    er = errorbar(time, tmp_rxn_data, errlow, errhigh);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    xline(2 - 0.5, '--g', 'linewidth', 2.0);
    xline(5 + 0.5, '--r', 'linewidth', 2.0);
    hold off
    xlabel('Time [days]', 'FontSize', 20);
    ylabel('Normalized Flux_{|SPT_{rxn}|}', 'FontSize', 20);
    set(gca, 'box', 'off');
    title([char(string(strrep(tmp_labels, '_', '-'))) '  (V_{mean} = ' char(num2str(round(tmp_mean_1, 5))) ')']);
    fontSize = 12; 
    set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
    fig_name_plot = fullfile(pwd, 'MainRxnFluxes', 'NormFluxStats', v, 'Barplots', [char(rxnIDs(j)) '.png']);
    folder_path = fileparts(fig_name_plot);
    if exist(folder_path, 'dir') ~= 7
        mkdir(folder_path);
    end
    saveas(gcf, fig_name_plot);
    
    rxn_names = [rxn_names; rxnIDs(j)];
end

figNo = figNo + 2;

% Save mean reaction flux data to a file (code commented out)

% Save flux data to a MAT file (code commented out)

fprintf('\n Done.  \n');
MatFileName = [v '_DMFA_fit_N_', char(string(N)), '_', datestr(now, 'HH-MM-dd-mmm-yyyy'), '.mat'];
% save(MatFileName)
