clear;
clc;
close all;

combination_ =  'Fig_2D';
myLegend = {'A', 'B', 'C', 'D'};
kvec = 1:12;
color_combination_5 = {'#1f77b4', '#9467bd', '#ff7f0e', '#e377c2'};
kkeycol = [color_combination_5, repelem(color_combination_5, 2)]';
N = 1000;
% files = dir(fullfile(pwd, ['N1*' 'Raw_data' '.xlsx']));
files = dir(fullfile('../data/', ['N1*' 'Raw_data' '.xlsx']));
for ii = 1:numel(files)
    type_D = {'N15'; 'N14N15'; 'N14' };
    if strncmp(files(ii).name, 'N14_N15', 7)
        type = type_D{2};
        y_label = ['^{14}N + ^{15}N'] ;
    elseif strncmp(files(ii).name, 'N15', 3)
        type = type_D{1};
        y_label = ['^{15}N'] ;
    elseif strncmp(files(ii).name, 'N14', 3)
        type = type_D{3};
        y_label = ['^{14}N'] ;
    end
    input_filename = fullfile('../data/', files(ii).name);
    [version_, date] = version;
    if year(date)>=2020
        file_sheets_names = sheetnames(input_filename);
    elseif year(date)<2020 && year(date)>=2011
        file_sheets_names =xl_xlsfinfo((pwd+"\"+input_filename)) ;
    elseif year(date)<2011
        fprintf('Code only works for versions of matlab released 2011b or after \n');
    end 
    for k = 1:numel(file_sheets_names)
        tmp_sheet = file_sheets_names{k};
        if year(date)>=2016
            tmp_tbl = readtable(input_filename , 'Sheet', tmp_sheet);
            tmp_tbl = tmp_tbl(1:38, 1:9);
        else
            fprintf('Code only works for versions of matlab released 2016b or after \n');
        end
        v = genvarname(tmp_sheet, who);
        assignin('base', v, tmp_tbl);
        raw_DATA{k} = tmp_tbl;
        fprintf('\n Cell Culture: %s  \t\n   ', v);
        met_list_X = table2cell(tmp_tbl(:, 2));
        tmp_id = find(not(cellfun('isempty', met_list_X)));
        met_list_X = met_list_X(tmp_id(1:end));
        no_mets_measured = numel(met_list_X);
        meas_X = table2cell(tmp_tbl(:, 3:end));
        tmp_id = find(not(cellfun('isempty', meas_X)));
        meas_X_vals = cell2mat(meas_X(tmp_id));
        meas_X_vals  = reshape(meas_X_vals, size(meas_X));
        conc_raw_data =  meas_X_vals;
        kb = 2;
        kf = 2;
        meas_X_vals = smoothdata(meas_X_vals, 2, 'movmean', [kb kf]);
        time_points_meas = size(meas_X_vals, 2);
        time = [0:1:time_points_meas-1]; 
        conc_data =  meas_X_vals;
        conc_data(isnan(conc_data))=0;
        conc_data(:, 1) = conc_raw_data(:, 1);  
        N_points = (time(end)*N)+1;
        time_interp = linspace(0, time(end), N_points);
        s = spline(time, conc_data, time_interp);
        s(s<0)= 0;
        My_time{k} =  time;
        My_time_interp{k} = time_interp;
        meas_X_VALS{k} = meas_X_vals;
        CONC_DATA{k} = conc_data;
        CONC_DATA_2{k} = conc_raw_data;
        S2{k} = s;
        eval(['clear ', v]);
    end
% combination = [type '_' combination_];
dir_name = fullfile('../results/', 'CodeFigures', combination_);
if ~isfolder(dir_name)
    mkdir(dir_name);
end
    for j = 1:no_mets_measured
        figure(1)
        for kk=1:numel(kvec(1:4))
            cc = kvec(kk);
            conc_data_use = CONC_DATA{cc};
            s2_use = S2{cc}; 
            tmp_mean(kk) = mean(conc_data_use(j,3:6));
            p(kk) = plot(My_time_interp{cc}, s2_use(j, :), 'color', kkeycol{kk, 1}, 'linewidth', 3.0);
            hold on        
        end 
        for kkk=kk+1:numel(kvec)
            cc = kvec(kkk);
            conc_data_use = CONC_DATA{cc};
            scatter(My_time{cc}, conc_data_use(j,:), 'Marker', '^', 'MarkerEdgeColor',  kkeycol{kkk, 1}, 'LineWidth', 1.5);
        end
        xline(2,'--g', 'linewidth', 2.5);
        xline(5,'--g', 'linewidth', 2.5);
        hold off
        xlim([0 5])
        fig = gcf;
        fig.Units = 'normalized';
%         fig.Color = 'none';
        fig.InvertHardcopy = 'off';
%         ax= gca ;
%         ax.Color = 'none';
        fontSize = 30; 
        set(gca, 'FontWeight', 'bold', 'FontSize', fontSize, 'FontName', 'Times New Roman');
        fig_name_plot_png = fullfile(dir_name, [char(met_list_X(j)), '.png']);
        saveas(gcf, fig_name_plot_png);
        
    end
end
