clc;clear;

addpath 'C:\Users\XUSHENG\software\croco_tools-forcast\Preprocessing_tools\';
addpath 'D:\application\SCS';
stafile = 'sta_out.nc';

% 读取每个观测站的位置
lon = ncread(stafile,'lon_rho');
lat = ncread(stafile,'lat_rho');
% station1[114,12](南海中心) station2[108,20](北部湾北)

csr = ncread('sta_out.nc','Cs_r');
h = ncread('sta_out.nc','h'); % 观测点水深
u = ncread(stafile,'u_eastward');
v = ncread(stafile,'v_northward');
w = ncread(stafile,'w');
temp = ncread(stafile,'temp');
salt = ncread(stafile,'salt');
rho = ncread(stafile,'rho');
zeta = ncread(stafile,'zeta'); %观测点水位

% sta1:110,18
% sta2:110.5,17.5
% sta3:111,17
max_time_length = size(zeta);
max_time_length = max_time_length(2);

%% 更换站点后此段重新运行
station_num = 1;
layer = 38; % 注意40是表层

endtime = max_time_length
starttime = 1
targettime = 5000

x = linspace(starttime,endtime,endtime-starttime+1);
y = linspace(0,0,endtime-starttime+1);
x = x';
y = y';
z = zlevs(h(station_num),zeta(station_num,targettime),4.5,1.5,5,40,'r',2);    
zz = (h(station_num)+zeta(station_num,targettime))*(linspace(1,1,40)+csr');
%% 流速
figure(1)
quiver(x,y,squeeze(u(layer,station_num,starttime:endtime)),squeeze(v(layer,station_num,starttime:endtime)));

%% 温盐
figure(2) % time-series in single layer
plot(x,squeeze(temp(layer,station_num,starttime:endtime)))
title(["temperature in station",num2str(station_num)])

figure(3) 
plot(x,squeeze(salt(layer,station_num,starttime:endtime)))
title(["salinity in station",num2str(station_num)])

%% 海表面高程
figure(4)
plot(x,squeeze(zeta(station_num,starttime:endtime)))
title(["ssh in station",num2str(station_num)])


%% vertical profiles
figure(5)
subplot(1,2,1)
plot(salt(:,station_num,targettime),z)
xlabel('盐度')
ylabel('深度')
subplot(1,2,2)
plot(temp(:,station_num,targettime),z)
xlabel('温度')
ylabel('深度')


