function [ WQI, Q ] = WQI(W)

FileName   = 'mean_station.mat';
FolderName = 'C:\Users\Lenovo\Desktop\software';
File       = fullfile(FolderName, FileName);
load(File);   
%% Sub-indexing

Si(1)= 120; % TA
Si(2)= 75; % Ca
Si(3)= 250; % Cl
Si(4)= 5; % Do
Si(5)= 100; % TH
Si(6)= 0.3; % Fe
Si(7)= 0.01; % Pb
Si(8)= 30; % Mg
Si(9)= 0.1; % Mn
Si(10)= 50; % TN
Si(11)= 7.5; % pH
Si(12)= 0.5; % TP
Si(13)= 300; % EC
Si(14)= 250; % So4
Si(15)= 500; % TDS
m_data = meann;
Q=zeros(size(m_data,1),size(m_data,2));
 for i=1:size(m_data,1)
        for j=1:size(m_data,2)
            if j==11
                if m_data(i,j)<=6.5
                    Q(i,j)=(m_data(i,j)-7)/(6.5-7)*100;
                elseif m_data(i,j)>=8.5
                    Q(i,j)=(m_data(i,j)-7)/(8.5-7)*100;
                else
                    Q(i,j)=(m_data(i,j))/(Si(j))*100;
                end
            elseif j==4
                Q(i,j)=(m_data(i,j)-14.6)/(Si(j)-14.6)*100;
            else
                Q(i,j)=(m_data(i,j))/(Si(j))*100;
            end
        end
 end


%% WQI Calculation

for j = 1:size(Q,1)
    WQI(j) = sum(W'.*Q(j,:))/sum(W);
end
end