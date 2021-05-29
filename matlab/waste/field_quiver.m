clear;
clc;

%-----------填写关键参数---------
lon1 = 105
lon2 = 124
lat1 = 5
lat2 = 24
lonD = lon2 - lon1
latD = lat2 - lat1

xi_rho = 121
eta_rho = 85
xi_u = xi_rho - 1
xi_v = xi_rho
eta_u = eta_rho 
eta_v = eta_rho -1

u_interval = 3
v_interval = 3

%-----------读取流场数据-------------------
%u=load('u_his_field.txt');
%v=load('v_his_field.txt');
u = ncread('roms_ini.nc','u')
v = ncread('roms_ini.nc','v')
u = u(:,:,20)
v = v(:,:,20)

    %稍移位使格点匹配
ux=u(:,1:eta_v);
vx=v(1:xi_u,:);

%----------异常剔除--------------
for lat = 1:eta_v;
    for lon = 1:xi_u;
        if abs(ux(lon,lat))>10;
            ux(lon,lat)=0;
            vx(lon,lat)=0;
        end;
        if abs(vx(lon,lat))>10;
            ux(lon,lat)=0;
            vx(lon,lat)=0;
        end;
    end;
end;

%----------为成品图添加经纬度信息--------------
for i = 1:xi_u;
    x(i)=lon1+(lonD/(xi_rho-1))*(i-1)
end;
for j = 1:eta_v;
    y(j)=lat1+(latD/(eta_rho-1))*(j-1)
end;


%--------------网格化---------------
[xx,yy]=meshgrid(x,y)
xx=xx'
yy=yy'



%{
for ii = 1:200;
    for jj = 1:224;
        quiver(xx(ii,jj),yy(ii,jj),u(ii,jj),v(ii,jj));
        hold on;
    end;
end;
%}

%--------------绘图--------------
    %读取岸线数据
load('coastline_f.mat')
quiver(xx(1:u_interval:xi_u,1:v_interval:eta_v),yy(1:u_interval:xi_u,1:v_interval:eta_v),ux(1:u_interval:xi_u,1:v_interval:eta_v),vx(1:u_interval:xi_u,1:v_interval:eta_v),2);
axis equal 
% ncst(isnan(ncst))=0;
hold on
plot(ncst(:,1),ncst(:,2));
axis equal
axis([lon1,lon2,lat1,lat2]);

for i = 1:80
for j = 1:159
if abs(u(i,j)) > 10
u(i,j) = 0
v(i,j) = 0
end
end
end

quiver(xx,yy,u,v)
axis equal