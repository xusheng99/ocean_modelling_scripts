%多张水位图合成gif用，临时代码
clc;clear;


hisname = 'D:\bhd_sea_model\wrongobc\scs_his.nc'
times = linspace(4345,4416,72);
times


    lon1 = 117;
    lon2 = 132;
    lat1 = 24;
    lat2 = 41;
    lonD = lon2 - lon1;
    latD = lat2 - lat1;

    xi_rho = 91;
    eta_rho = 123;
    xi_u = xi_rho - 1;
    xi_v = xi_rho;
    eta_u = eta_rho ;
    eta_v = eta_rho -1;

    u_interval = 3;    %流场图两点之间的间隔，单位：grid
    v_interval = 3;

    for dd = 1:xi_rho
        x(dd)=lon1+(lonD/(xi_rho-0))*(dd-1);
    end
    for ee = 1:eta_rho
        y(ee)=lat1+(latD/(eta_rho-0))*(ee-1);
    end
    [xx,yy]=meshgrid(x,y);
    xx=xx';
    yy=yy';
    
for time = times
    zeta = ncread(hisname,'zeta',[1,1,time],[Inf,Inf,1]);   %自由高程，时变数据
    s = pcolor(xx,yy,zeta)
    axis equal
    s.EdgeColor = 'none';
    %s.LineWidth = 0;
    s.FaceColor = 'interp';
    caxis([-2,2]);   %指定pcolor绘图的colorbar范围，便于多张图进行比较
    
    colorbar;
    title('水位')
    
    name = [num2str(time),'.png']
    
    fig = gcf     %获取图片文件句柄
    print(fig, '-dpng', '-r300',name)    %格式png，dpi1200，名称由uvfigname自动生成
    hold off
    clf    %关闭figure 1中的绘图，以防下一张图和前一张图重叠
    close(fig)
end
