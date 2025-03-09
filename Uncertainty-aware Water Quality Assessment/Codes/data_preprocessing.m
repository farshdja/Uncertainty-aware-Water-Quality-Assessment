clc;
clear;
% Water Quality Parameter: 'TA', 'Ca', 'Cl⁻', 'DO', 'TH', 'Fe', 'Pb', 'Mg', 'Mn', 'TN', 'pH', 'TP', 'EC', 'SO₄²⁻', 'TDS'

%%
sta = 1:27;
pre_WQD = [];      % Preprocessed Water Quality Data

for i = 1:length(sta)
[x] = processWaterQuality(sta(i));
    if i == 1
    pre_WQD = x;
    else
    pre_WQD = cat(1,pre_WQD,x);
    end
end

for i = 5:size(pre_WQD,2)
[transdat(:,i-4),lambda(i-4)] = boxcox(pre_WQD(:,i));
end

Preprocessed_Data = zscore(transdat);



%% ************************************************************************
%  ***********                  Sub-functions                   ***********
%  ******                                                            ******
%  **                                                                    **
%% ************************************************************************

function [Matrix] = processWaterQuality(sta)
    % Initialize output matrix
    Matrix = [];
    a=1:15;
    e=length(sta);
        for q = 1:length(a)
            clear q_data Q_data 
            % Read data from the Excel file
            import = readtable(num2str(sta), 'Sheet', a(q));
            Ce = table2cell(import);
            Strng = string(Ce);
            q_data = double(Strng(:, 4));
            
            %% Outlier detection
            TF = isoutlier(q_data, 'quartiles');
            out(e, q) = (sum(TF) / length(q_data)) * 100;
            q_data(TF) = nanmean(q_data);
            
            %% Find missing data + monthly temporal scale
            date_serial = datetime(import{:, 1},'InputFormat', 'mm/dd/yyyy');
            y_m = datevec(import{:, 1});
            y_m = unique(y_m(:, 1:3), 'rows');
            year = unique(y_m(:, 1));
            year = sort(year,'ascend');
            month = 1:12;
            un = unique(y_m(:, 1:2), 'rows');
            un2 = unique(date_serial);
            d = 1;
            
            for i = 1:length(un2)
                dt = find(un2(i) == date_serial);
                if i > length(q_data)
                    q_data(i) = nan;
                else
                    q_data(i) = mean(q_data(dt));
                end
            end
            q_data(length(un2)+1:end) = [];
            date_serial = unique(date_serial);
            
            for i = 1:length(year)
                for j = 1:length(month)
                    b = find(un(:, 1) == year(i) & un(:, 2) == month(j), 1);
                    if isempty(b)
                        Miss(d, :) = [year(i), month(j), 15];
                        if j == 1
                            k = find(y_m(:, 1) == year(i-1) & y_m(:, 2) == j+11);
                        else
                            k = find(y_m(:, 1) == year(i) & y_m(:, 2) == month(j-1));
                        end
                        y_m = rowBetween(y_m, Miss(d, :), k(end));
                        q_data = rowBetween(q_data, nan, k(end));
                        date_serial = rowBetween(date_serial, datetime(Miss(d, :)), k(end));
                        d = d + 1;
                    end
                end
            end
            
            days = date_serial(1) - date_serial;
            q_data(1) = fillmissing(q_data(1), 'constant', nanmean(q_data));
            q_data(end) = fillmissing(q_data(end), 'constant', nanmean(q_data));
            
            %% Missing data percentage
            w(e, q) = sum(isnan(q_data));
            per(e, q) = (w(e, q) / length(q_data)) * 100;
            
            %% Cubic interpolation and averaging
            int = interp1(days, q_data, days, 'pchip')';
            int(int <= 0) = mean(int(int > 0));
            
            for i = 1:length(int)
                adrs = find(y_m(i, 1) == y_m(:, 1) & y_m(i, 2) == y_m(:, 2));
                Q_data(adrs(end)) = mean(int(adrs));
                Q_data(i:adrs(end)-1) = nan;
            end
            
            final_data = Q_data(~isnan(Q_data))';
            final_ym = y_m(~isnan(Q_data), :);
            
            for i = 1:length(year)
                max_data(i) = max(final_data(final_ym(:, 1) == year(i)));
                min_data(i) = min(final_data(final_ym(:, 1) == year(i)));
                mean_data(i) = mean(final_data(final_ym(:, 1) == year(i)));
            end
            
            %% Store in output matrix
            Matrix(:, 1:3, e) = final_ym;
            Matrix(:, 4, e) = e;
            Matrix(:, q+4, e) = final_data;
        end
        Matrix(:,[10 11 13]) = Matrix(:,[10 11 13])./1000;        % ug/l to mg/l
    end


function [ result ] = rowBetween( mat,row,n )

mat = mat';
row = row';

combine = [mat(:,1:n) row mat(:,n+1:end)];

result = combine';

end