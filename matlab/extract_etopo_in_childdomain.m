cd '.'/read_roms_result/
%% 
lon = ncread('etopo1.nc','lon');
lat = ncread('etopo1.nc','lat');
h = ncread('etopo1.nc','topo');
%%
a = lon>110;
b = lon<115.7;
c = a.*b;
lons = find(c);

d = lat>20.5;
e = lat<23.2;
f = d.*e;
lats = find(f);

hh = h(lons,lats);

x = lon(lons);
y = lat(lats);
%%
savetxt('h.dat',x,y,hh)

%% 验证
clc;clear;
hh = load('h.dat');
% hh(hh>0) = NaN;
x = reshape(hh(:,1),[342,162]);
hhh = hh(:,3);
hhh(hhh<=0) = 0.1;
scatter(hh(:,1),hh(:,2),hhh);
contourf(hh');axis equal tight;colorbar;