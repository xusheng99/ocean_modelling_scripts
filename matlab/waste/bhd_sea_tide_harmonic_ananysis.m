clc;
clear;

addpath('C:\Users\XUSHENG\software\t_tide\');
addpath('C:\Users\XUSHENG\Desktop\roms相关脚本\matlab\');

%zetafile = 'D:\bhd_sea_model\bhd3_加细_未跑完\bhd_zeta.nc';
%gridfile = 'D:\bhd_sea_model\bhd3_加细_未跑完\roms_grid.nc';
%storage_file_name = 'bhd_sea_tide_harmonic_ananysis.txt';
zetafile = 'D:\scs_hainan_coarse_M2_alone\scs_his.nc';
gridfile = 'D:\scs_hainan_coarse_M2_alone\DATA_LARGE\roms_grid.nc';
storage_file_name = 'scs_hainan_coarse_M2_alone_harmonic_ananysis.txt';


zeta = ncread(zetafile,'zeta');
lon = ncread(zetafile,'lon_rho');
lat = ncread(zetafile,'lat_rho');
mask_rho = ncread(gridfile,'mask_rho');
rrr = size(zeta);
len_xi = rrr(1)
len_eta = rrr(2)
len_time = rrr(3)

% 创建空数组用于存储
tide = zeros(len_xi,len_eta,4);
phase = zeros(len_xi,len_eta,4);

%创建文件用于存储
outfile = fopen(storage_file_name,'w+');
for i = 1:len_xi
    for j = 1:len_eta
        
         if mask_rho(i,j) == 1
             [NAME,F,TIDE,XOUT]=t_tide(...
                 zeta(i,j,:),...
                 'interval',1,... %单位: 小时
                 'latitude',lat(i,j),...  %纬度控制科氏力
                 'start time',[2019 1 1 0 0 0]);  %第一个时刻
                 %'rayleigh',['M2';'S2';'K1';'O1']...  %仅分析四大分潮，所有字符串长度必须和最长的一个对齐
                 %);
                 
                for k=1:length(F)
                    % 使用不同的频率来区分不同的分潮
                    if abs(F(k)-0.0387306544)<1.0e-4          %o1 
                        tide(i,j,4)=TIDE(k,1);phase(i,j,4)=TIDE(k,3);
                    elseif abs(F(k)-0.0417807462)<1.0e-4      %k1
                        tide(i,j,3)=TIDE(k,1);phase(i,j,3)=TIDE(k,3); 
                    elseif abs(F(k)-0.0805114007)<1.0e-4   %m2
                        tide(i,j,1)=TIDE(k,1);phase(i,j,1)=TIDE(k,3);
                    elseif abs(F(k)-0.0833333333)<1.0e-4   %s2
                        tide(i,j,2)=TIDE(k,1);phase(i,j,2)=TIDE(k,3);
                    end 
                 end
     
         else  %陆地点
             tide(i,j,:) = NaN;phase(i,j,:) = NaN;
             
         end
         
            fprintf(...
                outfile,...
                '%15.10f %15.10f %20.15f %20.15f %20.15f %20.15f  %20.15f %20.15f  %20.15f %20.15f',...
                lon(i,j),lat(i,j),tide(i,j,1),phase(i,j,1),tide(i,j,2),phase(i,j,2),tide(i,j,3),phase(i,j,3),tide(i,j,4),phase(i,j,4)...
                );
            %fprintf(outfile,' %20.15f %20.15f  %20.15f %20.15f  %20.15f %20.15f '
            fprintf(outfile,'\n');
    end
end
fclose(outfile);

%% 创建绘图用网格
lon1 = 105.5
lon2 = 115
lat1 = 15
lat2 = 22

xi_rho = len_xi
eta_rho = len_eta

lonD = lon2 - lon1;
latD = lat2 - lat1;

for dd = 1:xi_rho
    x(dd)=lon1+(lonD/(xi_rho-0))*(dd-1);
end
for ee = 1:eta_rho
    y(ee)=lat1+(latD/(eta_rho-0))*(ee-1);
end

[xx,yy]=meshgrid(x,y);  %meshgrid for lon and lat
xx=xx';
yy=yy';
fontsize = 4
%% M2
tidevar = tide(:,:,1);
phasevar = phase(:,:,1);
titlecontent = 'M2分潮'

tmax = max(max(tidevar));
tmin = min(min(tidevar));
ddd = pcolor(xx,yy,tidevar);
ddd.EdgeColor = 'none';
ddd.FaceColor = 'interp';
caxis([tmin,tmax]);
colorbar
axis equal;
axis tight;
title(titlecontent);
hold on
[C,h] = contour(xx,yy,phasevar,linspace(0,330,11))
clabel(C,h,'FontSize',fontsize) %显示等高线上的数字

fig = gcf;     
print(fig, '-dpng', '-r600',[titlecontent,'.png'])   
hold off
clf    
close(fig)
