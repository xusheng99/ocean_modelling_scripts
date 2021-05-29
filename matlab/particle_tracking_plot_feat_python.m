clc;clear;

cd c:\users\xusheng\desktop\particle_tracking_results\plant
addpath C:\Users\XUSHENG\Documents\scripts\matlab_toolboxs\npy-matlab
% his_plot(1,'d:/application/SCS/his_0001.nc')
% hold on
% file=dir('*.dat');
%% 从npy文件读取
a = readNPY('result.npy');
[N,d,t] = size(a);
a(a == 19981231) = NaN;
for i = 1:N
    if a(i,1,1) == a(i,1,t)
        a(i,:,:) = NaN;
    end
end


%% 轨迹
for i = 1:20:N
    disp(i);
    plot(squeeze(a(i,1,:)),squeeze(a(i,2,:)));
    hold on
end
axis equal

%% 截图
for time = 1:2:t
    draw_coastline()
    scatter(a(:,1,time),a(:,2,time),0.15,'black');
    axis equal tight
    fig = gcf  
    print(fig, '-dpng', '-r600',[num2str(time,'%04d'),'.png'])   
    hold off
    clf    
    close(fig) 
end
%% 轨迹
for i = 1:100:length(file)
    disp(i)
    ttt = load(file(i,1).name);
    ttt(ttt == 19981231) = NaN;
%     if ttt(1,1) == ttt(end,1) & ttt(1,2) == ttt(end,2)
%         ttt(ttt>0)=NaN;
%     end
    plot(ttt(1:1:end,1),ttt(1:1:end,2))   
    hold on
end
axis equal

%% 终点
for i = 1:10:length(file)
    ttt = load(file(i,1).name);
    scatter(ttt(end,1),ttt(end,2));
    hold on
end

%% 释放位置标识
scatter(ttt(1,1),ttt(1,2),'filled','red');

%% 保存
fig = gcf  
print(fig, '-dpng', '-r600','trace_python.png')   
hold off
clf    
close(fig) 

%axis([min(ttt(:,1)zoom,max(ttt(:,1))zoom,min(ttt(:,2))zoom,max(ttt(:,2))*zoom])
