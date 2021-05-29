clear;clc;

%% 读取数据
addpath 'C:\Users\XUSHENG\Desktop\roms相关脚本\matlab'
cd 'C:\Users\XUSHENG\Desktop\particle_tracking_results\TDN\180-30'

[X1,Y1] = textread('initial.dat','%f %f');
[X2,Y2] = textread('final.dat','%f %f');

%% 参数
ii=500;jj=500;
title_ = 'FTLE_t0=20160629_dt=30days'    

%% 计算FTLE矩阵
qx=zeros(jj,ii);
qy=zeros(jj,ii);
hx=zeros(jj,ii);
hy=zeros(jj,ii);
FTLE=zeros(jj,ii);
for k=1:ii
     js=(k-1)*jj+1;
     je=k*jj;
    % SS(1:end,k)=ZH(hs,js:je)';
    qx(1:end,k)=X1(js:je);
    qy(1:end,k)=Y1(js:je);
    hx(1:end,k)=X2(js:je);
    hy(1:end,k)=Y2(js:je);
end

for i=2:ii-1
  for j=2:jj-1  
     b1=(hx(j,i+1)-hx(j,i-1))/(qx(j,i+1)-qx(j,i-1));
     b2=(hx(j+1,i)-hx(j-1,i))/(qy(j+1,i)-qy(j-1,i));
     b3=(hy(j,i+1)-hy(j,i-1))/(qx(j,i+1)-qx(j,i-1));
     b4=(hy(j+1,i)-hy(j-1,i))/(qy(j+1,i)-qy(j-1,i));
     B=[b1,b2;b3,b4];
     C=B'*B;
     tezheng=eig(C);
     lama=max(tezheng);
     FTLE(j,i)=log(lama^0.5)/240;
    end
end

%% 处理异常值
mask_rho = ncread('his_0001.nc','mask_rho');
FTLE(isinf(FTLE))=NaN; % 无穷点
for i = 1:ii
    for j = 1:jj
        if X2((i-1)*ii + j) > 1000 % 离开开边界的粒子标记为NaN
%             ^ ^
%             | | ^    3列
%             | | |    2列
%             | | | .. 1列
%             1 2 3 行 是这样子排列的
            FTLE(ii,jj) = NaN;
        end
    end    
end

histogram(FTLE)
disp("根据直方图结果设置合适的阈值")
disp("press any key to continue")
pause();

thresh = 0.03 % 阈值
FTLE(FTLE>thresh) = NaN;
FTLE(FTLE == 0) = NaN;
%% 保存FTLE 用于其他程序作图
% lon = linspace(106,124,120);
% lat = linspace(6,24,120);
% savetxt('FTLE.dat',lon,lat,FTLE');

%% 绘图
figure(1);
min = min(min(FTLE))
max = max(max(FTLE))
minX = 105.5
maxX = 124.5
minY = 5.5
maxY = 24.5
xx = linspace(minX,maxX,ii);  
yy = linspace(minY,maxY,jj);
[xxx,yyy] = meshgrid(xx,yy);

% 使用基于三角的双线性插值 
Z=griddata(qx,qy,FTLE,xxx,yyy,'cubic');
[c,hh]=contourf(xxx,yyy,Z,100); %画分色图 

% 擦掉等值线 
set(hh,'linestyle','none');
caxis([min,max]);
colorbar ;
hold on

% 绘制岸线
coastlinefile = 'all_china_sea_i.mat'    %指定岸线文件
load(coastlinefile);  %岸线
plot(ncst(:,1),ncst(:,2));   
hold on
axis equal;
axis ([105,125,5,25])

% 增加标题
title(title_);
