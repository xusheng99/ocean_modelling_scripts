clear;clc;


addpath 'C:\Users\XUSHENG\Desktop\roms相关脚本\matlab'
addpath 'C:\Users\XUSHENG\Desktop\particle_tracking_results\'
addpath 'd:\application\SCS\'

[X1,Y1] = textread('initial.dat','%f %f');
[X2,Y2] = textread('final.dat','%f %f');
% X1(X1 == 19981231) = NaN;
% Y1(Y1 == 19981231) = NaN;
% X2(X2 == 19981231) = NaN;
% Y2(Y2 == 19981231) = NaN;

% %%
% N = length(X1);  %原始释放粒子数目
% % x1 = find(X1);
% x2 = find(X2);  % 判断粒子是否撞开边界
% 
% len = length(x2);  %x2长度一定短于x1,因为粒子只会损失不会增添
% XX1 = linspace(0,0,len);
% YY1 = linspace(0,0,len);
% XX2 = linspace(0,0,len);
% YY2 = linspace(0,0,len);
% XX1 = XX1';
% YY1 = YY1';
% XX2 = XX2';
% YY2 = YY2';
% ddd = []
% 
% 
% count = 1  %不要把count写在循环里面, 否则每一次循环都会将count初始化为1.
% for i = 1:N
%     if ismember(i,x2)
%         XX1(count) = X1(i);
%         YY1(count) = Y1(i);
%         XX2(count) = X2(i);
%         YY2(count) = Y2(i);
%         count = count + 1;
%     else
%         ddd = [ddd,i];
%     end   
% end
% X1 = XX1;X2 = XX2;Y1 = YY1;Y2 = YY2;

%%


 ii=120;jj=120;%lzw
            
%  ii=2081;jj=2131;%bh
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
%      SS(j,i)~=0&&
%      if(qx(j,i+1)~=0&&qx(j,i-1)~=0&&qy(j+1,i)~=0&&qy(j-1,i)~=0&&qx(j,i+1)-qx(j,i-1)~=0&&qy(j+1,i)-qy(j-1,i)~=0)
%       if(qx(j,i+1)-qx(j,i-1)~=0&&qy(j+1,i)-qy(j-1,i)~=0)
     b1=(hx(j,i+1)-hx(j,i-1))/(qx(j,i+1)-qx(j,i-1));
     b2=(hx(j+1,i)-hx(j-1,i))/(qy(j+1,i)-qy(j-1,i));
     b3=(hy(j,i+1)-hy(j,i-1))/(qx(j,i+1)-qx(j,i-1));
     b4=(hy(j+1,i)-hy(j-1,i))/(qy(j+1,i)-qy(j-1,i));
     B=[b1,b2;b3,b4];
     C=B'*B;
     tezheng=eig(C);
     lama=max(tezheng);
     FTLE(j,i)=log(lama^0.5)/240;
%      else
%          FTLS(j,i)=0;  
%      end
  
    end
end

FTLE(isinf(FTLE))=0;
% FTLE2=reshape(FTLE,250000,1);
% 
% max=max(max(FTLE));
% min=min(min(FTLE));
% xishu=max-0.2*(max-min);
% [a b]=find(FTLE>xishu);
% L=zeros(jj,ii);
% L(a,b)=1;

% xx = 25000:100:145000;  %bhw
% yy = 4215000:100:4355000;
figure
minX = 106
maxX = 124
minY = 6
maxY = 24
xx = linspace(minX,maxX,ii);   %lzw
yy = linspace(minY,maxY,jj);
% xx = 26000:100:442000;
% yy = 4110000:100:4536000;

[xxx,yyy] = meshgrid(xx,yy);

Z=griddata(qx,qy,FTLE,xxx,yyy,'cubic');%使用基于三角的双线性插值 
[c,hh]=contourf(xxx,yyy,Z,100); %画分色图 

% Z=griddata(qx,qy,L,xxx,yyy,'linear');%使用基于三角的双线性插值 
% [c,hh]=contourf(xxx,yyy,Z,100); %画分色图 


set(hh,'linestyle','none');%擦掉等值线 
      caxis([0.000,0.02]);
% %    caxis([-0.002,0.006]);
 colorbar ;
%   load jiegou;
%   colormap(jiegou);


% hold on;
% X3=ZX(hs,1:end);Y3=ZY(hs,1:end);
% X4=ZX(he,1:end);Y4=ZY(he,1:end);
% Qu=(X4-X3)/900;Qv=(Y4-Y3)/900;
% 
% uu = griddata(X3,Y3,Qu,xxx,yyy);
% vv = griddata(X3,Y3,Qv,xxx,yyy);
% density=6;
% verts=streamslice(xxx,yyy,uu,vv,density);
% 
% 
axis equal tight;
hold on

coastlinefile = 'all_china_sea_i.mat'    %指定岸线文件

load(coastlinefile);  %岸线
plot(ncst(:,1),ncst(:,2));   
axis equal
axis([105,119,14,22]);   
hold on

title('FTLE test 20190101-15days')

% fig = gcf;     
% print(fig, '-dpng', '-r600','ftletest.png')   
% hold off
% clf    
% close(fig)
YU

%{
%%
hold on;
A=importdata('I:\张燕伟\output2018\coast.bln');
Bx=A(:,1);
By=A(:,2);
for j=1:length(Bx)
if ((By(j)==0)|(By(j)==1))
bx=Bx(j+1:j+Bx(j));by=By(j+1:j+Bx(j)); %注释掉首行
plot(bx,by,'color','b'); %(绘制隐形岸线)
patch(bx,by,[0.7,0.7,0.7]);
hold on;
end
end
% hold on;
% A=importdata('岸线\island.bln');
% Bx=A(:,1);
% By=A(:,2);
% for j=1:length(Bx)
% if ((By(j)==0)|(By(j)==1))
% bx=Bx(j+1:j+Bx(j));by=By(j+1:j+Bx(j)); %注释掉首行
% plot(bx,by,'color','w'); %(绘制隐形岸线)
% patch(bx,by,[0.7,0.7,0.7]);
% hold on;
% end
% end


 axis equal;
%    axis([26000 443000 4100000 4545000]);
%     axis([138000 300000 4110000 4250000]);%渤海湾
    axis([100000 350000 4100000 4250000]);%莱州湾
%  set(gca,'ytick',[4100750,4211800,4322850,4433900,4544950]);
% set(gca, 'YTickLabel', {'37°N', '38°N', '39°N', '40°N', '41°N'});
% set(gca,'xtick',[66800,143400,240000,326600,413200]);
% set(gca, 'XTickLabel', {'118°E', '119°E', '120°E', '121°E', '122°E'});
  
 box on;
   saveas (gcf,['1','.tif']);
%%
%%存储FTLE
[m,n]=size(qx);
num=m*n;

xxx=reshape(qx,num,1);
yyy=reshape(qy,num,1);
f=reshape(FTLE,num,1);


fid=fopen('2018ftle.txt','w');
fprintf(fid,'%12.2f %12.2f %12.5f\n',[xxx,yyy,f]');
fclose(fid);
%}
