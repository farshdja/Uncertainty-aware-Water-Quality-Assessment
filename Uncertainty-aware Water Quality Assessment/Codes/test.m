clc;
clear;


data = readtable('1');
data.Date = datetime(data{:,1}, 'InputFormat', 'MM/dd/yyyy'); % Adjust format as needed
[~, strings, ~] = xlsread(num2str(sta(e)), a(q));
            import = readtable(num2str(sta(e)), 'Sheet', a(q));
            Ce = table2cell(import);
            Strng = string(Ce);