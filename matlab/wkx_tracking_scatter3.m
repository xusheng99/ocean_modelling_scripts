clear;clc;

addpath 'D:\lg'

x=ncread('rz_out_lag3d.nc','x');
y=ncread('rz_out_lag3d.nc','y');
z=ncread('rz_out_lag3d.nc','z');

addpath 'C:\Users\XUSHENG\Desktop\';



scatter3(x(:,1000),y(:,1000),z(:,1000),'filled');
axis([510000,630000,3820000,3950000,-30,-1])

