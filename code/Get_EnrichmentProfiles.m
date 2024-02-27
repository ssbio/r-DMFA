clear; clc; close all;
N = 1000;
incorporation = 'N15';
combination_ = 'Fig_2C';
myLegend = {'A', 'B', 'C', 'D'};
kvec = 1:12;
color_combination = {'#1f77b4', '#9467bd', '#ff7f0e', '#e377c2'};
kkeycol = [color_combination, repelem(color_combination, 2)]';
kkeycol_2 = {'#1f77b4', '#e377c2'}';
files = dir(fullfile('../data/', ['N1*' 'Raw_data' '.xlsx']));
raw_DATA_num = {}; raw_DATA_den = {};
for ii = 1:numel(files)
    type_D = {'N15'; 'N14N15'; 'N14'};
    type = matchType(files(ii).name, type_D);
    input_filename = ['../data/' files(ii).name]; 
    input_filename_G = strrep(input_filename, 'Raw_data', 'Trend_data');
    [~, date] = version;
    file_sheets_names = loadSheetNames(input_filename, date);
    for k = 1:numel(file_sheets_names)
        tmp_sheet = file_sheets_names{k};
        if year(date) >= 2016
            tmp_tbl = readtable(input_filename, 'Sheet', tmp_sheet);
            tmp_tbl = tmp_tbl(1:38, 1:9);
            raw_DATA{k} = tmp_tbl;
            [met_list_X, meas_X_vals, conc_raw_data] = processTable(tmp_tbl);
            no_mets_measured = numel(met_list_X);
            time_points_meas = size(meas_X_vals, 2);
            time = 0:1:time_points_meas-1;
            conc_data = smoothAndCleanData(meas_X_vals, conc_raw_data);
            N_points = (time(end)*N)+1;
            time_interp = linspace(0, time(end), N_points);
            s = spline(time, conc_data, time_interp);
            s(s<0)= 0;
            My_time{k} = time;
            My_time_interp{k} = time_interp;
            meas_X_VALS{k} = meas_X_vals;
            CONC_DATA{k} = conc_data;
            CONC_DATA_2{k} = conc_raw_data;
            S2{k} = s;
        else
            fprintf('Code only works for versions of Matlab released 2016b or after \n');
        end
        v = [genvarname(tmp_sheet, who) '_' type];
        if strcmp(extractAfter(v, '_'), type_D{1}) && strcmp(incorporation, type_D{1})
            raw_DATA_num = [raw_DATA_num, conc_raw_data];
        elseif strcmp(extractAfter(v, '_'), type_D{3}) && strcmp(incorporation, type_D{3})
            raw_DATA_num = [raw_DATA_num, conc_raw_data];
        elseif strcmp(extractAfter(v, '_'), type_D{2})
            raw_DATA_den = [raw_DATA_den, conc_raw_data];
        end
    end      
end
% combination = [type '_' combination_];
dir_name = fullfile('../results/', 'CodeFigures', combination_);
if ~isfolder(dir_name)
    mkdir(dir_name);
end
num_data = length(raw_DATA_num);
R_CONC_DATA = cell(1, num_data);
R_CONC_DATA_2 = cell(1, num_data);
R_S2 = cell(1, num_data);
R_My_time = cell(1, num_data);
R_My_time_interp = cell(1, num_data);
for kk = 1:num_data
    DATA_ratio{kk} = raw_DATA_num{kk}./raw_DATA_den{kk};
    DATA_ratio{kk}(isinf(DATA_ratio{kk})) = NaN;
    conc_raw_data_r =  DATA_ratio{kk};
    conc_raw_data_r(:, 1) = strcmp(incorporation, type_D{1}) * 0 + strcmp(incorporation, type_D{3});
    conc_data_R = smoothAndCleanData(conc_raw_data_r, conc_raw_data_r);
    time_points_meas = size(conc_raw_data_r, 2);
    time = 0:time_points_meas-1;
    N_points = (time(end)*N)+1;
    time_interp = linspace(0, time(end), N_points);
    s_r = spline(time, conc_data_R , time_interp);
    s_r = max(min(s_r,1),0);
    s_r(:, 1) = conc_data_R(:, 1);
    R_CONC_DATA{kk} = conc_data_R;
    R_CONC_DATA_2{kk} = conc_raw_data_r;
    R_S2{kk} = s_r;
    R_My_time{kk} = time;
    R_My_time_interp{kk} = time_interp;
end
[~, selectedDays] = ismember(time, time_interp);
for j = 1:no_mets_measured
    figure(1)
    tmp_keep = [];
    for kk = 1:2
        tmp_keep = [tmp_keep; R_S2{kk}(j, :)];
    end
    N14_tmp_keep_mean = mean(tmp_keep, 1);
    pp(2) = plot(R_My_time_interp{2}, N14_tmp_keep_mean, 'color', kkeycol_2{1, 1}, 'linewidth', 3.0);
    hold on 
   tmp_keep = [];
    for kk = 3:4
        tmp_keep = [tmp_keep; R_S2{kk}(j, :)];
    end
    N15_tmp_keep_mean = mean(tmp_keep, 1);
    pp(1) = plot(R_My_time_interp{4}, N15_tmp_keep_mean, 'color', kkeycol_2{2, 1}, 'linewidth', 3.0);
    hold on 
    tmp_keep = [];
    for kk = 5:8
        tmp_keep = [tmp_keep; R_CONC_DATA{kk}(j, :)];
    end
    tmp_lb = N14_tmp_keep_mean(:, selectedDays) - std(tmp_keep, [], 1);
    tmp_ub = N14_tmp_keep_mean(:, selectedDays) + std(tmp_keep, [], 1);
    tmp_interp = linspace(0, time(end), (time(end)*100)+1);
    s_tmp_lb = spline(time, tmp_lb , tmp_interp);
    s_tmp_lb = max(min(s_tmp_lb,1),0);
    s_tmp_lb(:, 1) = tmp_lb(:, 1);
    s_tmp_ub = spline(time, tmp_ub , tmp_interp);
    s_tmp_ub = max(min(s_tmp_ub,1),0);
    s_tmp_ub(:, 1) = tmp_ub(:, 1);
    patch([tmp_interp fliplr(tmp_interp)], [s_tmp_ub fliplr(s_tmp_lb)], 'p', 'FaceColor', kkeycol_2{1, 1}, 'FaceAlpha',0.1);

    tmp_keep = [];
    for kk = 9:12
        tmp_keep = [tmp_keep; R_CONC_DATA{kk}(j, :)];
    end
    tmp_lb = N15_tmp_keep_mean(:, selectedDays) - std(tmp_keep, [], 1);
    tmp_ub = N15_tmp_keep_mean(:, selectedDays) + std(tmp_keep, [], 1);
    s_tmp_lb = spline(time, tmp_lb , tmp_interp);
    s_tmp_lb = max(min(s_tmp_lb,1),0);
    s_tmp_lb(:, 1) = tmp_lb(:, 1);
    s_tmp_ub = spline(time, tmp_ub , tmp_interp);
    s_tmp_ub = max(min(s_tmp_ub,1),0);
    s_tmp_ub(:, 1) = tmp_ub(:, 1);
    patch([tmp_interp fliplr(tmp_interp)], [s_tmp_ub fliplr(s_tmp_lb)], 'p', 'FaceColor',kkeycol_2{1, 1}, 'FaceAlpha',0.1);
    xline(2,'--g', 'linewidth', 2.5);
    xline(5,'--r', 'linewidth', 2.5);
    hold off
    legend(pp, '^{15}N-Labeling', '^{14}N-Labeling', 'Location', 'east', 'Orientation', 'horizontal', 'Box', 'off');
    xlabel('Time [days]');
    ylabel([incorporation ' / (N14 + N15)'])
    title(met_list_X{j});
    fontSize = 12; 
    set(gca, 'FontWeight', 'bold', 'FontSize', fontSize);
    fig_name_plot = fullfile(dir_name, [met_list_X{j}, '.png']);
    saveas(gcf, fig_name_plot);
end
function type = matchType(filename, type_D)
    if strncmp(filename, 'N14_N15', 7)
        type = type_D{2};
    elseif strncmp(filename, 'N15', 3)
        type = type_D{1};
    elseif strncmp(filename, 'N14', 3)
        type = type_D{3};
    end
end
function file_sheets_names = loadSheetNames(input_filename, date)
    if year(date)>=2020
        file_sheets_names = sheetnames(input_filename);
    elseif year(date)<2020 && year(date)>=2011
        file_sheets_names =xl_xlsfinfo((pwd+"\"+input_filename)) ;
    elseif year(date)<2011
        fprintf('Code only works for versions of Matlab released 2011b or after \n');
    end 
end
function [met_list_X, meas_X_vals, conc_raw_data] = processTable(tmp_tbl)
    met_list_X = table2cell(tmp_tbl(:, 2));
    tmp_id = find(not(cellfun('isempty', met_list_X)));
    met_list_X = met_list_X(tmp_id(1:end));
    meas_X = table2cell(tmp_tbl(:, 3:end));
    tmp_id = find(not(cellfun('isempty', meas_X)));
    meas_X_vals = cell2mat(meas_X(tmp_id));
    meas_X_vals  = reshape(meas_X_vals, size(meas_X));
    conc_raw_data = meas_X_vals;
end
function conc_data = smoothAndCleanData(meas_X_vals, conc_raw_data)
    kb = 2;
    kf = 2;
    meas_X_vals = smoothdata(meas_X_vals, 2, 'movmean', [kb kf]);
    conc_data = meas_X_vals;
    conc_data(isnan(conc_data))=0;
    conc_data(:, 1) = conc_raw_data(:, 1);  
end
