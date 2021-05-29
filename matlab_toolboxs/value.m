% function value(lon,lat,depth,time,array)
% lon: degrees
% lat: degrees
% depth: meters
% time: timesteps in the target 4darray
% array: 4darray, like water_u(time,depth,lat,lon)
clc;clear;
addpath(['C:\Users\XUSHENG\software\croco_tools-forcast\Preprocessing_tools\'])

file = 'C:\Users\XUSHENG\Desktop\roms\hycom.nc'
time = ncread(file,'time');
lev = ncread(file,'LEV');
lat = ncread(file,'latitude');
lon = ncread(file,'longitude');
u = ncread(file,'water_u'); % size = time,depth,lat,lon
v = ncread(file,'water_v');

% 目标位置的mask是1,  才开展后续的插值计算, 否则直接给个nan;
% 在上述条件满足的情况下:
% 如果目标位置src中有值,则正常插值,如果目标位置src中无值,采用上下左右四点的值.

theta_s = 4.5
theta_b = 1.5
hc = 1.
N = 40

% 读取grid基本信息
grid = 'C:\Users\XUSHENG\Desktop\roms\Grd.nc';
mask_rho = ncread(grid,'mask_rho');
lon_rho = ncread(grid,'lon_rho');
lat_rho = ncread(grid,'lat_rho');
mask_u = ncread(grid,'mask_u');
lon_u = ncread(grid,'lon_u');
lat_u = ncread(grid,'lat_u');
mask_v = ncread(grid,'mask_v');
lon_v = ncread(grid,'lon_v');
lat_v = ncread(grid,'lat_v');
h = ncread(grid,'h'); %全正


% 网格大小
[xi_rho,eta_rho] = size(mask_rho);
[xi_u,eta_u] = size(mask_u);
[xi_v,eta_v] = size(mask_v);

% 插出u网格和v网格上的水深
[x_rho,y_rho] = meshgrid(1:xi_rho,1:eta_rho);
x_rho=x_rho';y_rho = y_rho';
[x_u,y_u] = meshgrid(1:xi_u,1:eta_u);
x_u = x_u';y_u = y_u';
x_u = x_u+0.5;
[x_v,y_v] = meshgrid(1:xi_v,1:eta_v);
x_v = x_v';y_v = y_v';
y_v = y_v+0.5;
h_u = griddata(x_rho,y_rho,h,x_u,y_u,'linear');
h_v = griddata(x_rho,y_rho,h,x_v,y_v,'linear'); %griddata支持不规则待提取区域.


% z = zlevs(h,zeta,theta_s,theta_b,hc,N,'r',2);  
% 坐标变换数组
z_rho = -zlevs(h,0,theta_s,theta_b,hc,N,'r',2);  
z_rho = permute(z_rho,[2,3,1]);   
z_u = -zlevs(h_u,0,theta_s,theta_b,hc,N,'r',2);  
z_u = permute(z_u,[2,3,1]);   
z_v = -zlevs(h_v,0,theta_s,theta_b,hc,N,'r',2);  
z_v = permute(z_v,[2,3,1]);   

% hycom源数据坐标数组
[src_lon,src_lat,src_lev] = meshgrid(lon,lat,lev);
[lon_u_mesh,lat_u_mesh,wastes] = meshgrid(lon_u(:,1),lat_u(1,:)',[1:N]');
lon_u_mesh = permute(lon_u_mesh,[2,1,3]);
lat_u_mesh = permute(lat_u_mesh,[2,1,3]);

i = 120;
j = 120;
k = 1;
lon_u_mesh(i,j,k)
lat_u_mesh(i,j,k)
z_u(i,j,k)
h_u(i,j)

138 137
137 138

% a= lon_u_mesh(:,:,2)
% b = src_lon(:,:,2);
% 
% i = 480;
% j = 301;

% 插值
uu = griddata(src_lon,src_lat,src_lev,squeeze(u(:,:,:,1)),lon_u_mesh,lat_u_mesh,z_u,'linear');

uu = griddata(src_lon,src_lat,src_lev,squeeze(u(:,:,:,1)),114.9583,14.9167,4183.5,'linear');


% inperp1(srcloc,srcvalue,targetloc,method) -> targetvalue
u_ini = zeros(size(z));
v_ini = zeros(size(z));



