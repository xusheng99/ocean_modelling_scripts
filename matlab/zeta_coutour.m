clc
clear

lon1 = 105
lon2 = 124
lat1 = 5
lat2 = 24
lonD = lon2 - lon1
latD = lat2 - lat1
xi_rho = 229
eta_rho = 238
xi_u = xi_rho - 1
xi_v = xi_rho
eta_u = eta_rho 
eta_v = eta_rho -1
u_interval = 4
v_interval = 4

 %地理坐标系
for i = 1:xi_rho;
    x(i)=lon1+(lonD/(xi_rho-1))*(i-1)
end;
for j = 1:eta_rho;
    y(j)=lat1+(latD/(eta_rho-1))*(j-1)
end;
[xx,yy]=meshgrid(x,y)
xx=xx'
yy=yy'


%水高程
zeta = load('60_zeta.txt')
levels = [-2.9,-2,-1,0,1,2]
[c,zeta]=contourf(xx,yy,zeta,levels)
clabel(c,zeta,'FontSize',5,'Color','blue')
hold on

 %岸线
load('coastline_f.mat')
% ncst(isnan(ncst))=0;
plot(ncst(:,1),ncst(:,2),'color','red');    %两个同长度一维数组绘制连线图 scatter为散点，plot则将这些点顺序连接
axis equal
axis([lon1,lon2,lat1,lat2]);    %岸线数据范围是全球，我们只画目标区域的


    