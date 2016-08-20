%% Read dfs2 and compute yearly minimum

clear all; clc;
addpath('C:\Users\DianaL\Documents\matlab\DHIMatlabToolbox_20130222\mbin');
%% Minimum flow
name_file = ['C:\Users\DianaL\Desktop\POSTER_AGU\poster\observations\Q_data\HOBE_Q-obs_usedAHL1990-2014.dfs0'];
dfs0  = dfsTSO(name_file);
year = [1990:1:2014]; ndays = yeardays(year);nd = cumsum(ndays);
station = [20,21,82]
%corrF=[];corrA=[];corrJ=[];corrFm

for st = 1:1:3
dataQ = double(dfs0(st+1));
data0 = NaN(length(year),366);

for ny = 0:1:length(year)-1 
    if ny==0 
        st0 = 1; 
    else st0 = nd(ny)+1; 
    end
    data0(ny+1,1:ndays(ny+1)) = dataQ(st0:nd(ny+1));
end
 file = ['obs_station_',num2str(station(st)),'.txt'];   
 save(file,'data0','-ascii');
end



