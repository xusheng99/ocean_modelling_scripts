clear;
clc;

addpath(['C:\Users\XUSHENG\software\m_map'])
%sd_list = [2,3,5,8,10,20,30,50,100,150,200,300,500,1000]    %把想要作图的深度写在这里，单位是m，前提是数据已经放在“作图用文件”里面了
%sd_list=[10,50,100,200,300,400,500]
%time_list = [720,4344]    %hours since 2019-01-01-00:00
sd_list = [10]
time_list = [720]

%{
sd = 10
t = 720
%}
for i = sd_list
    for j = time_list
        myfunction(i,j)
    end
end

function [] = myfunction(sd,t)
%----组合出图片名字----
sd = num2str(sd);
ttt = t
t = num2str(t);
uname = ['u_',sd,'meters',t,'hours_flowfield.txt'];
vname = ['v_',sd,'meters',t,'hours_flowfield.txt'];
Vname = ['Velocity_',sd,'meters',t,'hours_flowfield.txt'];
saltname = ['salt_',sd,'meters',t,'hours_field.txt'];
tempname = ['temp_',sd,'meters',t,'hours_field.txt'];
if ttt == 720
    date = '2019-01-31-00-00'
end
if ttt == 4344
    date = '2019-07-01-00:00'
end
uvfigname = ['flowfield_',sd,'meters_',date,'.png'];    %uv的名字
Vfigname = ['V',sd,'_meters',t,'hours_img.png'];
saltfigname = ['salt',sd,'_meters',t,'hours_img.png']
tempfigname = ['temp',sd,'_meters',t,'hours_img.png']

ux = load(uname);
vx = load(vname);
V = load(Vname);
salt = load(saltname);
temp = load(tempname);
%----参数指定----
lon1 = 105;
lon2 = 124;
lat1 = 5;
lat2 = 24;
lonD = lon2 - lon1;
latD = lat2 - lat1;

xi_rho = 229;
eta_rho = 238;
xi_u = xi_rho - 1;
xi_v = xi_rho;
eta_u = eta_rho ;
eta_v = eta_rho -1;

u_interval = 3;
v_interval = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %右上偏移
%{
uux = zeros(228,237);
vvx = uux;
Vx = uux;
for i = 1:228
    for j = 1:237
        uux(i,j) = ux(i+1,j);
        vvx(i,j) = vx(i+1,j);
        Vx(i,j) = V(i+1,j);
    end
end
u2 = [ux(1:228,1),uux];
temp = [ux(229,1),ux(229,2:238)];
ux = [u2;temp];
v2 = [vx(1:228,1),vvx];
temp = [vx(229,1),vx(229,2:238)];
vx = [v2;temp];
V2 = [V(1:228,1),Vx];
temp = [V(229,1),V(229,2:238)];
V = [V2;temp];
%向流场向右上角平移了一段，勉强和岸线对上了
%造成数据整体移动的原因还没搞明白
%}
%{
    %上偏移
uxx = zeros(223,238);
vxx = uxx;
Vx = uxx;
for i = 1:223
    for j = 1:238
        uxx(i,j) = ux(i+6,j);
        vxx(i,j) = vx(i+6,j);
        Vx(i,j) = V(i+6,j);
    end
end
ux = [uxx;ux(224:229,:)];
vx = [vxx;vx(224:229,:)];
V = [Vx;V(224:229,:)];
%}
    
%----加入地理坐标网格-----
for i = 1:xi_rho
    x(i)=lon1+(lonD/(xi_rho-0))*(i-1);
end
for j = 1:eta_rho
    y(j)=lat1+(latD/(eta_rho-0))*(j-1);
end
[xx,yy]=meshgrid(x,y);
xx=xx';
yy=yy';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------制图部分--------


%%%%%%%%第一张图:流场箭头图%%%%%%%%
    
%{
%岸线
load('coastline_f.mat');
% ncst(isnan(ncst))=0;
plot(ncst(:,1),ncst(:,2));    %两个同长度一维数组绘制连线图 scatter为散点，plot则将这些点顺序连接
axis equal
axis([lon1,lon2,lat1,lat2]);    %岸线数据范围是全球，我们只画目标区域的
hold on
%}
%{
    %等高线
m = linspace(1,570,229);
n = linspace(1,570,238);
[mm,nn]=meshgrid(m,n);
mm = mm';
nn = nn';
h = load('etoposcs.txt');
X = 1:570;
Y = 1:570;
[XX,YY] = meshgrid(X,Y);
H = interp2(XX,YY,h,mm,nn,'linear');
H = -H;
levels = [10,500,2500]%,20,50,200,500,1000,1500,2000,2500,2800]
[c,h]=contour(xx,yy,H,levels);
clabel(c,h,'FontSize',5,'Color','red')
axis equal
hold on    %hold on的意思是之后的所有作图均叠加在hold on所指定的前一幅图上
           %只需要写一次hold on就够了
%}           
    %流场
quiver(xx(1:u_interval:xi_u,1:v_interval:eta_v),yy(1:u_interval:xi_u,1:v_interval:eta_v),ux(1:u_interval:xi_u,1:v_interval:eta_v),vx(1:u_interval:xi_u,1:v_interval:eta_v),2);
axis equal 
hold on

    %保存图片
fig = gcf     %获取图片文件句柄
print(fig, '-dpng', '-r600', uvfigname)    %格式png，dpi1200，名称由uvfigname自动生成
hold off
clf    %关闭figure 1中的绘图，以防下一张图和前一张图重叠
close(fig)
%%%%%%%%第一张图绘制完毕%%%%%%%%


%%%%%%%%绘制第二张图：速率等值线图%%%%%%%%
%{
%岸线
load('coastline_f.mat')
% ncst(isnan(ncst))=0;
plot(ncst(:,1),ncst(:,2));    %两个同长度一维数组绘制连线图 scatter为散点，plot则将这些点顺序连接
axis equal
axis([lon1,lon2,lat1,lat2]);    %岸线数据范围是全球，我们只画目标区域的
hold on
%}
    %速率等值线
levels = linspace(0,2,10)
contourf(xx,yy,V,levels)
colorbar;
axis equal
hold on
    %保存图片
fig = gcf     %获取图片文件句柄
print(fig, '-dpng', '-r600', Vfigname)    %格式png，dpi1200，名称由Vfigname自动生成
hold off
clf    
close(fig)
%%%%%%%%第二张图绘制完成%%%%%%%%%

%%%%%%%%绘制第三张图：温度场%%%%%%%%
contourf(xx,yy,temp)
colorbar;
axis equal
hold on
fig = gcf     %获取图片文件句柄
print(fig, '-dpng', '-r600', tempfigname)    %格式png，dpi1200，名称由Vfigname自动生成
hold off
clf    
close(fig) 
%%%%%%%%绘制完毕%%%%%%%%

%%%%%%%%绘制第四张图：盐度图'
levels = linspace(32,35,16)
contourf(xx,yy,salt,levels)
colorbar;
axis equal
colorbar;
hold on
fig = gcf     %获取图片文件句柄
print(fig, '-dpng', '-r600', saltfigname)    %格式png，dpi1200，名称由Vfigname自动生成
hold off
clf    
close(fig) 
%%%%%%%%绘制完毕%%%%%%%%

end


