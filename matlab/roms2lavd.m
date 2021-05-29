clc;clear;

%% 读取roms中的表层流场
addpath d:/application/SCS/

% nc文件使用如下语句生成:
%  !!!大乌龙!!!
%  -0.9875是sigma底层!!! 应该用-0.0125!!!
%cdo -select,name=u_eastward,v_northward,level=-0.9875 $(ls his*) uv_surface.nc4
%cdo -setmisstoc,0 uv_surface.nc4 uv_surface_0.nc4
ncfile = 'uv_surface_0.nc4';

u = squeeze(ncread(ncfile,'u_eastward')); %lon,lat,time
v = squeeze(ncread(ncfile,'v_northward'));

ocean_time = ncread(ncfile,'ocean_time');
lon_rho = ncread(ncfile,'lon_rho');
lat_rho = ncread(ncfile,'lat_rho');

%% 调整格式

ocean_time = ocean_time(1:6:end);
time = ocean_time/86400;
disp("time变量已转换完成");

u = permute(u,[2,1,3]); %lat,lon,time is required by LAVD tool
v = permute(v,[2,1,3]);
UT = u(:,:,1:6:end); 
VT = v(:,:,1:6:end);
disp("速度分量已转换完成");

lon = lon_rho(:,1);
lat = lat_rho(1,:);
lat = lat';
disp("坐标变量已转换完成");

%% 目标区域 不能包含陆地

disp("正在保存变量到: roms_surface_flowfield_for_LAVD.mat");
save('roms_surface_flowfield_for_LAVD.mat','time','lon','lat','UT','VT')
disp("保存完成");