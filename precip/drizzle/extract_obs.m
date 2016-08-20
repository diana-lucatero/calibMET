%% Bias correction and downscaling seasonal EPS - observation

clear all; clc;

inf = xlsread('ensem_info.xlsx','info_obs');
inf2 = xlsread('grid_obs.xlsx','grid_obs3');
%% Load observations from dfs2 file

NET.addAssembly('DHI.Generic.MikeZero.DFS');
import DHI.Generic.MikeZero.DFS.*;
dfs2b = DfsFileFactory.Dfs2FileOpen('10kmgrid_DynCorrPrecip_1990-2014.dfs2');

for i = 1:1:dfs2b.FileInfo.TimeAxis.NumberOfTimeSteps
    data1 = double(dfs2b.ReadItemTimeStep(1,i-1).To2DArray()); 
    data(:,:,i) = flip(data1,2);
end

data(data < 0) = NaN; %% Replaces negative values with NaN

years = 1990:2015;
days = yeardays(years);
finish = cumsum(days);
start = [1,finish+1];

dataf = NaN(1,367);
for j = 1:1:inf2(end,1)
    datafi = NaN(25,367);
    for fc = 0:1:25-1
        dataf = [years(fc+1) squeeze(data(inf2(j,4),inf2(j,5),start(fc+1):finish(fc+1)))'];
        if days(fc+1) == 365
            dataf = [dataf,NaN];
        else
            dataf = dataf;
        end
        datafi(fc+1,:) = dataf;
     end
    
  save(['observation/obs_station_',num2str(inf2(j,1)),'.txt'], 'datafi', '-ASCII');
  disp(['obs ',num2str(j)]);
end

