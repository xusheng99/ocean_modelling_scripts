clear;clc;

cd 'C:\Users\XUSHENG\Desktop\particle_tracking_results\LTRANS_output'

%% 从csv加载数据
file = dir('*.csv');
for i = 1:length(file)
    filecontent = load(file(i,1).name);
    lon_ = filecontent(:,3);
    lat_ = filecontent(:,4);
    height_ = filecontent(:,1);
    lon(:,i) = lon_;
    lat(:,i) = lat_;
    depth(:,i) = height_;
end
lonsize = size(lon);
N = lonsize(1)
T = lonsize(2)



%% 从nc加载数据
lon = ncread('./output.nc','lon');
lat = ncread('./output.nc','lat');
depth = ncread('./output.nc','depth');      
age = ncread('./output.nc','age');
lonsize = size(lon);
N = lonsize(1)
T = lonsize(2)


maxx =  72
for time = 1:2:maxx
%% 岸线
load('C:\Users\XUSHENG\documents\scripts\matlab\all_china_sea_i.mat');  %连续岸线
plot(ncst(:,1),ncst(:,2));    
axis equal
hold on

draw_coastline
%% 
% time = 200;
scatter(lon(:,time),lat(:,time),0.15,'black');

r = 0.2
lon_stride = max(lon(:,maxx)) - min(lon(:,maxx));
lat_stride = max(lat(:,maxx)) - min(lat(:,maxx));

axis([min(lon(:,maxx))-r*lon_stride,max(lon(:,maxx))+r*lon_stride,...
    min(lat(:,maxx))-r*lat_stride,max(lat(:,maxx))+r*lat_stride]) % 合适的显示范围

fig = gcf  
print(fig, '-dpng', '-r600', [num2str(time,'%04d'),'.png'])   
hold off
clf    
close(fig)
end

%% 二维轨迹示意图
%对于每一个点, 做连线
for i  = 1:100:N
    plot(lon(i,:),lat(i,:))
    hold on
end
axis equal


%% 三维动图
% particle_num = 1
% for time = 1:T
%     scatter3(lon(particle_num,time),lat(particle_num,time),depth(particle_num,time),'filled')
%     hold on 
%     %{
%     fig = gcf  
%     print(fig, '-dpng', '-r600', [num2str(time),'.png'])   
%     hold off
%     clf    
%     close(fig)
%     %}
% end
% box on



