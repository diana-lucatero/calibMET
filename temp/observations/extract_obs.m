%% Extract data observations

tic

clear all; clc;

inf = xlsread('F:\ECMWF_Seasonal_data\BC_lt_T\ensem_info.xlsx','info_T');
inf2 = xlsread('F:\ECMWF_Seasonal_data\2-layer-filer\grid_obs.xlsx','grid_obs3');

%% Load observations from dfs2 file

NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;
dfs2b = DfsFileFactory.Dfs2FileOpen('F:\ECMWF_Seasonal_data\2-layer-filer\20kmgrid_Ta_1990-2014.dfs2');

for i = 1:1:dfs2b.FileInfo.TimeAxis.NumberOfTimeSteps
    data1 = double(dfs2b.ReadItemTimeStep(1,i-1).To2DArray()); 
    data(:,:,i) = flip(data1,2);
end

data(  data <= -1.0e-30 & data >= -2.0e-30 )=NaN; %% Replaces negative values with NaN
%data(  data >= -1.0e-30 | data <= -2.0e-30 )=NaN;
% data(  data >= -2.0e-30 | data <= -1.0e-30 )=NaN;
% data(  abs(data) <= 2.0e-5 )=NaN;

for j = 1:1:inf2(end,1)
    datafi = [];
    for fc = 0:1:293-1 % Until start 2014-06
        if fc == 0
            dataf = [fc+1 NaN squeeze(data(inf2(j,4),inf2(j,5),inf(fc+1,29):inf(fc+1,29)+213))'];
            datafi = [datafi; dataf];
        else
        dataf = [fc+1 squeeze(data(inf2(j,4),inf2(j,5),inf(fc+1,29):inf(fc+1,29)+214))'];
        datafi = [datafi; dataf];
     end
    end
    
  save(['obs_station_',num2str(inf2(j,1)),'.txt'], 'datafi', '-ASCII');
   disp(['obs ',num2str(j)]);
end


toc