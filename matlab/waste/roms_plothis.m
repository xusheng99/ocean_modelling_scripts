function roms_plot(time,ncfilename)
% 函数文件和普通子函数的区别就是函数文件把function写在第一行, 且一个文件只定义一个函数

targettime =  time 
hisname = ncfilename

addpath(['C:\Users\XUSHENG\software\croco_tools-forcast\Preprocessing_tools\'])
addpath(['C:\Users\XUSHENG\software\croco_tools-forcast\Visualization_tools\'])

%% 可能需要指定的参数
%gridname = 'D:\bhd_sea_model\bhd3_refine_undone\roms_grid.nc'
%hisname = 'D:\bhd_sea_model\bhd3_refine_undone\out_his.nc'
%gridname = 'D:\bhd_sea_model\DATA_LARGE\roms_grid.nc'
%hisname = 'D:\bhd_sea_model\out_his.nc'
%gridname = 'D:\scs_hainan\DATA_LARGE\roms_grid.nc';
%hisname = 'D:\scs_model\scs_model_only_uv_in_rho_grid\scs_his1.nc';

ifsave = 0;  %是否自动保存生成的图片

coastlinefile = 'all_china_sea_i.mat';    %指定岸线文件

field_interval = 2   %流场图两点之间的间隔，单位：grid

%zarray = [1;2;3;4;5;6;7;8;9;10;12;15;20;30;40;50;60;70;80;90;100;110;120;130;140;150;160;170;180;190;200;250;300;400;500;600;700;800;900;1000;1200;1500;1800;2000;2500];
zarray = [-1];    %########关心的z分层深度（meters）写在这里########## 正值代表meters, -1代表海表
zarray = flip(zarray);    %翻转，使之顺序和zz以及value一致：从海底到海面


plot_switch = [0,0,0,1,0]    %作图开关:   [温度、盐度、流场、流速、水位]
%选0不作图
%对于温盐,1作contourf，2作pcolor
%对于流场箭头图,1做真实岸线,2做网格岸线
%第三个开关"流场"将速率映射到箭头长度
%第四个开关"流速"将速率映射到颜色

%% 读取原始数据
%######## DOMAIN ########
timelen = length(ncread(hisname,'ocean_time'));
if timelen < time
    error('nc文件中不存在指定的时间步长! 终止运行!')  %手动抛出异常
end
lon_rho = ncread(hisname,'lon_rho');  %范围，经纬度
lat_rho = ncread(hisname,'lat_rho');
lon1 = min(min(lon_rho));  
lon2 = max(max(lon_rho));    
lat1 = min(min(lat_rho));
lat2 = max(max(lat_rho));
%########垂向变换系数##########
theta_s = ncread(hisname,'theta_s')
theta_b = ncread(hisname,'theta_b')
hc = ncread(hisname,'hc')
N = length(ncread(hisname,'s_rho'))
%######## 温盐流水位等数据 ########
depth = ncread(hisname,'h');    %水深数据，时不变数据
start = [1,1,1,targettime];    %lon-lat-sigmalevel-timeslice
count = [Inf,Inf,Inf,1];    %这里不是end，是个数计数
stride= [1,1,1,1];    %stride is one
u = ncread(hisname,'u_eastward',start,count);
v = ncread(hisname,'v_northward',start,count);
temp = ncread(hisname,'temp',start,count);    %注意，从温度容易看出，1层是底层，40层是表层，不要搞反
salt = ncread(hisname,'salt',start,count);
zeta = ncread(hisname,'zeta',[1,1,targettime],[Inf,Inf,1]);   %自由高程，时变数据
maskr = ncread(hisname,'mask_rho',[1,1],[Inf,Inf]);
%######## zlevs函数 ########
z = zlevs(depth,zeta,theta_s,theta_b,hc,N,'r',2);    %%%%%！！！！关键的变换数组！！！！%%%%
%置换维度，符合个人记忆习惯
z = permute(z,[2,3,1]);   %数字表示原来的维度，例如原来的第一维放到了后来的第三维
z = -z;

%% 创建空数组用于存储
[lenxi,leneta,N] = size(u);  %xi方向网格数目，eta方向网格数目，sigma层层数
lenz = length(zarray);    %待插值的z分层层数
tempz = zeros(lenxi,leneta,lenz);    %创建同大小空数组，存放插值后数据
saltz = zeros(lenxi,leneta,lenz);
uz = zeros(lenxi,leneta,lenz);
vz = zeros(lenxi,leneta,lenz);
Vz = zeros(lenxi,leneta,lenz);
uz_nor = zeros(lenxi,leneta,lenz);
vz_nor = zeros(lenxi,leneta,lenz);

%% 插值运算
if zarray > 0  %关注z分层某深度
    for i = 1:lenxi
        for j = 1:leneta 
            if maskr(i,j) == 1    % 对于陆点，zeta全为NaN，会导致插值出来的z在该点全为NaN，导致后续错误
                % 准备location_vector和value_vector
                zz = z(i,j,:);
                zz = squeeze(zz);
                tempp = temp(i,j,:);    % 某一点的垂向剖面数据，临时变量
                tempp = squeeze(tempp);    % 删除多余的长度为1的维度
                saltp = salt(i,j,:);
                saltp = squeeze(saltp);
                up = u(i,j,:);
                vp = v(i,j,:);
                up = squeeze(up);
                vp = squeeze(vp);
                % Interpolation
                tempz(i,j,:) = interp1(zz,tempp,zarray,'linear');
                saltz(i,j,:) = interp1(zz,saltp,zarray,'linear');
                uz(i,j,:) = interp1(zz,up,zarray,'linear');
                vz(i,j,:) = interp1(zz,vp,zarray,'linear');
                Vz(i,j,:) = sqrt(uz(i,j,:).*uz(i,j,:)+vz(i,j,:).*vz(i,j,:));     
            else    % 陆点参与插值运算会报错
                tempz(i,j,:) = NaN;
                saltz(i,j,:) = NaN;
                uz(i,j,:) = NaN;
                vz(i,j,:) = NaN;
                Vz(i,j,:) = NaN;
            end
        end
    end
else  % 关注sigma表面
    sigma = N + 1 + zarray(1);
    tempz = temp(:,:,sigma);   
    tempz = squeeze(tempz); 
    saltz = salt(:,:,sigma);
    saltz = squeeze(saltz);
    uz = u(:,:,sigma);
    vz = v(:,:,sigma);
    uz = squeeze(uz);
    vz = squeeze(vz);
    Vz = sqrt(uz.*uz+vz.*vz); 
end
 
if plot_switch(4) ~= 0
    % 正规化的速率，方向不变，大小缩放到1
    for i = 1:lenxi
        for j = 1:leneta
            for k = 1:lenz
                uz_nor(i,j,k) = uz(i,j,k)/Vz(i,j,k);
                vz_nor(i,j,k) = vz(i,j,k)/Vz(i,j,k);
            end
        end
    end
end

%% 作图
len = length(zarray);
for dkdk = 1:len  %对每一个指定要作图的水深...
    ux = uz(:,:,dkdk);  %读取该层数据
    vx = vz(:,:,dkdk);
    tempzz = tempz(:,:,dkdk);
    saltzz = saltz(:,:,dkdk);
    Vzz = Vz(:,:,dkdk);
    ux_nor = uz_nor(:,:,dkdk);
    vx_nor = vz_nor(:,:,dkdk);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_interval = field_interval;   
    v_interval = field_interval;
    
    xi_rho = lenxi;
    eta_rho = leneta; %底下有不少用这个名字的地方，不再统一改成leneta了
    x = linspace(lon1,lon2,lenxi);
    y = linspace(lat1,lat2,leneta);
    
    %{ 
    典型错误方法
    for i = 1:xi_rho
    x(i)=lon1+(lonD/(xi_rho-0))*(i-1);
    end
    for j = 1:eta_rho
        y(j)=lat1+(latD/(eta_rho-0))*(j-1);
    end
    %}
    
    [xx,yy]=meshgrid(x,y);  %meshgrid for lon and lat
    xx=xx';
    yy=yy';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%%
    if plot_switch(3) == 1 || plot_switch(3) == 2   %流场
        if zarray > 0
            uvfigname = [num2str(targettime,'%04d'),'timestep_',num2str(zarray(dkdk)),'meters_flowfield.png']
        else
            uvfigname = [num2str(targettime,'%04d'),'timestep_','SeaSurfaceFlowfield.png']
        end
        
        if plot_switch(3) == 1
            figure(3)
            load(coastlinefile);  %连续岸线
            % ncst(isnan(ncst))=0;
            plot(ncst(:,1),ncst(:,2));    %两个同长度一维数组绘制连线图 scatter为散点，plot则将这些点顺序连接
            hold on
        end
        if plot_switch(3) == 2
            figure(3);
            ffp = pcolor(xx,yy,maskr);  %网格岸线
            ffp.EdgeColor = 'none';
            axis equal 
            colormap(gray(2))
            hold on
        end
        q = quiver(xx(1:u_interval:xi_rho,1:v_interval:eta_rho),yy(1:u_interval:xi_rho,1:v_interval:eta_rho),ux(1:u_interval:xi_rho,1:v_interval:eta_rho),vx(1:u_interval:xi_rho,1:v_interval:eta_rho),2);
        q.ShowArrowHead = 'on';
        q.LineWidth = 0.1  %线太粗会盖住尖锐的箭头头部
        axis tight
        axis equal
        axis([lon1,lon2,lat1,lat2]);    %岸线数据范围是全球，我们只画目标区域的
        title([num2str(zarray(dkdk)),'meters-flowfield'])
        
        if ifsave == 1
            fig = gcf    
            print(fig, '-dpng', '-r600', uvfigname)   
            hold off
            clf   
            close(fig) 
        end
    end
        
        
 %%  温度
    if plot_switch(1) == 1    %温度contourf
        if zarray > 0
            tempfigname = [num2str(targettime,'%04d'),'timestep_',num2str(zarray(dkdk)),'meters_temp.png']
        else
            tempfigname = [num2str(targettime,'%04d'),'timestep_','SeaSurfaceTemperature.png']
        end
        figure(1);
        tempmin = min(min(tempzz));
        tempmax = max(max(tempzz));
        templevels = linspace(tempmin,tempmax,16);
        [C,h] = contourf(xx,yy,tempzz,templevels)
        clabel(C,h,'FontSize',5) %显示等高线上的数字
        colorbar;
        axis tight
        axis equal
        title([num2str(zarray(dkdk)),'meters-temperature'])
        hold on
         if ifsave == 1
            fig = gcf    
            print(fig, '-dpng', '-r600', uvfigname)   
            hold off
            clf   
            close(fig) 
        end
    end
    if plot_switch(1) == 2    %温度网格图
        if zarray > 0
            tempfigname = [num2str(targettime,'%04d'),'timestep_',num2str(zarray(dkdk)),'meters_temp.png']
        else
            tempfigname = [num2str(targettime,'%04d'),'timestep_','SeaSurfaceTemperature.png']
        end
        figure(1);
        tempmin = min(min(tempzz));
        tempmax = max(max(tempzz));
        ttp = pcolor(xx,yy,tempzz)
        ttp.EdgeColor = 'none';
        ttp.FaceColor = 'interp';
        caxis([tempmin,tempmax]);   
        axis tight
        axis equal
        colorbar
        box on
        title([num2str(zarray(dkdk)),'meters-temperature'])
        hold on
        if ifsave == 1
            fig = gcf    
            print(fig, '-dpng', '-r600', uvfigname)   
            hold off
            clf   
            close(fig) 
        end
    end
%%  盐度
    if plot_switch(2) == 1    %盐度contourf
        if zarray > 0
            saltfigname = [num2str(targettime,'%04d'),'timestep_',num2str(zarray(dkdk)),'meters_salt.png']
        else
            saltfigname = [num2str(targettime,'%04d'),'timestep_','SeaSurfaceSalinity.png']
        end
        figure(2);
        levels = linspace(26,35,32)
        contourf(xx,yy,saltzz,levels)
        colorbar;
        axis tight
        axis equal
        colorbar;
        title([num2str(zarray(dkdk)),'meters-salinity'])
        hold on
        if ifsave == 1
            fig = gcf    
            print(fig, '-dpng', '-r600', uvfigname)   
            hold off
            clf   
            close(fig) 
        end
    end
    if plot_switch(2) == 2    %盐度网格图
        if zarray > 0
            saltfigname = [num2str(targettime,'%04d'),'timestep_',num2str(zarray(dkdk)),'meters_salt.png']
        else
            saltfigname = [num2str(targettime,'%04d'),'timestep_','SeaSurfaceSalinity.png']
        end
        figure(2);
        stp = pcolor(xx,yy,saltzz)
        stp.EdgeColor = 'none';
        stp.FaceColor = 'interp';
        caxis([20,36]);   
        axis tight
        axis equal
        colorbar;
        title([num2str(zarray(dkdk)),'meters-salinity'])
        box on
        hold on
        if ifsave == 1
            fig = gcf    
            print(fig, '-dpng', '-r600', uvfigname)   
            hold off
            clf   
            close(fig) 
        end
    end
 %%  速度场
    if plot_switch(4) == 1 || plot_switch(4) == 2
        if zarray > 0
            Vfigname = [num2str(targettime,'%04d'),'timestep_',num2str(zarray(dkdk)),'meters_current.png']
            title_ = [num2str(zarray(dkdk)),'meters_current']
        else
            Vfigname = [num2str(targettime,'%04d'),'timestep_','SeaSurfaceCurrent.png']
            title_ = ['Sea Surface Current']
        end
        if plot_switch(4) == 2
            figure(4);
            vtp = pcolor(xx,yy,Vzz);
            vtp.EdgeColor = 'none';
            vtp.FaceColor = 'interp';
            caxis([min(min(Vzz)),max(max(Vzz))]);
            axis tight
            axis equal
            colorbar;
            load('VelocityColorMap.mat')
            colormap(VelocityColorMap)
            hold on
            title(title_)
            box on
            hold on
            quiversize = 0.5  %箭头的大小，默认0.4
            q = quiver(xx(1:u_interval:xi_rho,1:v_interval:eta_rho),yy(1:u_interval:xi_rho,1:v_interval:eta_rho),ux_nor(1:u_interval:xi_rho,1:v_interval:eta_rho),vx_nor(1:u_interval:xi_rho,1:v_interval:eta_rho),quiversize);
            q.Color = 'black';
            q.ShowArrowHead = 'on';
            q.LineWidth = 0.1  %线太粗会盖住尖锐的箭头头部
        end
        if plot_switch(4) == 1
            figure(3)
            load(coastlinefile);  %连续岸线
            plot(ncst(:,1),ncst(:,2));    %两个同长度一维数组绘制连线图 scatter为散点，plot则将这些点顺序连接
            hold on
            title(title_)
            box on
            hold on
            quiversize = 0.5  %箭头的大小，默认0.4
            q = quiver(xx(1:u_interval:xi_rho,1:v_interval:eta_rho),yy(1:u_interval:xi_rho,1:v_interval:eta_rho),ux_nor(1:u_interval:xi_rho,1:v_interval:eta_rho),vx_nor(1:u_interval:xi_rho,1:v_interval:eta_rho),quiversize);            q.ShowArrowHead = 'on';
            q.LineWidth = 0.1  %线太粗会盖住尖锐的箭头头部
            axis equal
            axis([lon1,lon2,lat1,lat2]);
        end
       
        
        if ifsave == 1
            fig = gcf    
            print(fig, '-dpng', '-r600', uvfigname)   
            hold off
            clf   
            close(fig) 
        end
        
    end
%% 水位
    if plot_switch(5) == 1    %水位pcolor
        figure(5);
        s = pcolor(xx,yy,zeta);
        s.EdgeColor = 'none';
       % s.LineWidth = 0;
        s.FaceColor = 'interp';
        caxis([-2,2]);   %指定pcolor绘图的colorbar范围，便于多张图进行比较
        axis equal
        colorbar;
        title('水位')
    end
end
end


