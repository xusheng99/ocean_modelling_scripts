clc
clear

cd 'd:\application\SCS\'

ff = 'flt_track.nc';
start = [1,1]
count = [Inf,216]
lon = ncread(ff,'lon',start,count);
lat = ncread(ff,'lat',start,count);
% lon(lon == 0) = NaN;
% lat(lat == 0) = NaN;
depth = ncread(ff,'depth');
thesize = size(lon);
N = thesize(1)
T = thesize(2)
lonnew = [];
latnew = [];
for i = 1:N % 对每一个粒子
    if lon(i,T) ~= 0 % 如果最后时刻没有进入开边界, 则保留.
        lonnew = [lonnew,lon(i,:)];
        latnew = [latnew,lat(i,:)];
    end
end


dlmwrite('initial.dat',[lon(:,2),lat(:,2)]);
dlmwrite('final.dat',[lon(:,end),lat(:,end)]);
            

ftle_t0 = 100;
ftle_dt = ftle_t0+48;
lon_start = lon(:,ftle_t0);
lon_end = lon(:,ftle_dt);
lat_start = lat(:,ftle_t0);
lat_end = lat(:,ftle_dt);
ftle_start = [lon_start,lat_start];
ftle_end = [lon_end,lat_end];
ftle_all = [ftle_start,ftle_end];
dlmwrite('all.dat',ftle_all);
dlmwrite('.dat',ftle_start)
dlmwrite('2.dat',ftle_end)

%% 指定粒子的三维轨迹
flt_num = 404;  % 1~N任取一个粒子
x = lon(flt_num,:);
y = lat(flt_num,:);
z = depth(flt_num,:);

figure(1)
scatter3(x,y,z)
axis equal

%% 全体粒子的水平轨迹
figure(3)
for i = 1:10:N
    x2 = lon(i,:);
    y2 = lat(i,:);
    plot(x2,y2)
    hold on
end
axis equal

figure(4)
scatter(lon(:,end),lat(:,end))
axis equal
axis tight



his_plot(1,'his_0001.nc');