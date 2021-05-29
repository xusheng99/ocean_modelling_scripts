clc;
clear;

addpath('C:\Users\XUSHENG\software\t_tide\');
addpath('C:\Users\XUSHENG\Desktop\roms相关脚本\matlab\');

zetafile = 'D:\application\SCS\2016half\zetafile.nc4';
grid = 'D:\application\SCS\2016half\his_0001.nc';

zeta = ncread(zetafile,'zeta');
lon = ncread(zetafile,'lon_rho');
lat = ncread(zetafile,'lat_rho');
mask_rho = ncread(grid,'mask_rho');
rrr = size(zeta);
len_xi = rrr(1)
len_eta = rrr(2)
len_time = rrr(3)
lon1 = min(min(lon))
lon2 = max(max(lon))
lat1 = min(min(lat))
lat2 = max(max(lat))

%% 从保存的文件中读取
cd d:/application/SCS/2016half/
load('scs_tide_amplitude')
load('scs_tide_phase')

%% 从zeta文件中计算
% 创建空数组用于存储
tide = zeros(len_xi,len_eta,4);
phase = zeros(len_xi,len_eta,4);

%创建文件用于存储
outfile = fopen('scs_tide_harmonic_analysis.txt','w+');

for i = 1:len_xi
    for j = 1:len_eta
        
         if mask_rho(i,j) == 1
             [NAME,F,TIDE,XOUT]=t_tide(...
                 zeta(i,j,:),...
                 'interval',1,...
                 'latitude',lat(i,j),...
                 'start time',[2019 1 1 0 0 0]);
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
fontsize = 7

levels = [30,60,120,180,240,300,330]

%% 模式输出同潮图
figure(1)
% M2
subplot(2,2,1)
tidevar = tide(:,:,1);
phasevar = phase(:,:,1);
phasevar(phasevar>350 | phasevar<10) = NaN; % 避免在360和0之间巨大的梯度
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
[C,h] = contour(xx,yy,phasevar,levels)
clabel(C,h,'FontSize',fontsize) %显示等高线上的数字

% S2
subplot(2,2,2)
tidevar = tide(:,:,2);
phasevar = phase(:,:,2);
phasevar(phasevar>350 | phasevar<10) = NaN;
titlecontent = 'S2分潮'

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
[C,h] = contour(xx,yy,phasevar,levels)
clabel(C,h,'FontSize',fontsize) %显示等高线上的数字

% K1
subplot(2,2,3)
tidevar = tide(:,:,3);
phasevar = phase(:,:,3);
phasevar(phasevar>350 | phasevar<10) = NaN;
titlecontent = 'K1分潮'

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
[C,h] = contour(xx,yy,phasevar,levels)
clabel(C,h,'FontSize',fontsize) %显示等高线上的数字

% O1
subplot(2,2,4)
tidevar = tide(:,:,4);
phasevar = phase(:,:,4);
phasevar(phasevar>350 | phasevar<10) = NaN;
titlecontent = 'O1分潮'

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
[C,h] = contour(xx,yy,phasevar,levels)
clabel(C,h,'FontSize',fontsize) %显示等高线上的数字


