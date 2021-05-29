    %% --------Particle Track based on the FVCOM Flow Filed----------
%%%%%%%%%%%%%%%%%%%%%%采用四阶龙格库塔算法%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Yi+1 = Yi +1/6*K1 +1/3*K2+1/3*K3+1/6*K4;                      %%% 
%%%       K1 = hf(Xi,Yi);                                               %%%
%%%       K2 = hf(Xi+h/2, Yi+K1/2);                                     %%%
%%%       K3 = hf(Xi+h/2, Yi+K2/2);                                     %%%
%%%       K4 = hf(Xi+h, Yi+K3);                                         %%%
%%%       追踪模型引入了随机扩散                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   文件格式说明   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       主要的输入文件有：释放粒子的原始位置文件‘lag_ini.dat’，         %%%
%%%       与fvcom粒子追踪初始粒子文件格式相同,                             %%%

%%%       数据读取和赋值时已经将无用数值删除                               %%%
%%%       fcvom 流场数据文件，直接读取系列文件；                           %%%
%%%       分区各点坐标，格式为surfer可划分区的bln文件；                    %%%
%%%       主要的输出文件有：释放粒子的到达位置文件‘lag_out.dat’，         %%%
%%%       格式为每个粒子到达位置的时间循环结束后，循环粒子个数；             %%%
%%%       粒子各时间到达的分区、到达数量、源汇比例等                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
fprintf('Particle Track Model based on FVCOM Start Now!\n');
%% ---------------------------------------------------------------
%输入文件
% Particle_initial_file='2018岸线结果\粒子数量敏感分析\bohai_lag_ini_jy200m.dat';%fvcom运行追踪的输入文件
Particle_initial_file='初始位置2500.dat';%fvcom运行追踪的输入文件
% Particle_initial_file='2018岸线结果无风\粒子数量敏感分析\释放粒子路径文件\lag_nondif001_243.dat';%fvcom运行追踪的输入文件%需要热启动
% UV_folder_name='jingliu_UV\';%旧版fvcom输出结果，主要利用uv
UV_folder_name='2018UV365day\';%新版fvcom输出结果，主要利用uv
SubRegion='bohaifenqu.bln';%可在surfer中画出分区的数据文件
A = importdata('2018bhx.dat');
xv=A(:,1);
yv=A(:,2);
%输出文件
% Particle_Output_file=fopen('lag_results\lag_dif_out.dat','w');
%% ---------------------------------------------------------------
%参数设置
SubRegion_Number=751;%分区数量，如果下面两个数据的系数分区之间不同，可以用文件给出
Particle_SubRegion_Number=2500.*ones(SubRegion_Number,1);%初始时每个分区粒子数量
SubRegion_Point_Number=5.*ones(SubRegion_Number,1);%一个分区所需端点数量，首尾一致
TimeStep=3600;%时间步长，单位为s
TimeStep_h=TimeStep/3600;%时间步长，单位为h，用来提供时间循环步长
% StartTime=165;%释放开始时间，单位为h
StartTime=2879;%+90*24-1;%释放开始时间，单位为h,这里加是因为5天结果已经保存，要以这个数据为起点继续追踪，类似热启动
EndTime  =3240;%追踪结束时间，单位为h
OpenDiffusion=1;%扩散计算控制，1为开启，0为关闭

%% ---------------------------------------------------------------
%粒子初始位置赋值
Particle_initial=textread(Particle_initial_file);
Particle_initial_x=Particle_initial(1:end,1);
Particle_initial_y=Particle_initial(1:end,2);
% Particle_initial_x=Particle_initial(:,3);%热起动时
% Particle_initial_y=Particle_initial(:,4);
Particle_sum_number=length(Particle_initial_x);
Particle_initial_SubRegion=zeros(Particle_sum_number,1);
Particle_number=1;
for Number=1:SubRegion_Number;
     Line0=Particle_number;
     Line1=Particle_SubRegion_Number(Number,1);
     Line2=Particle_number+Line1-1;
     Particle_number=Line2+1;
     Particle_initial_SubRegion(Line0:Line2,1)=Number;
end
%速度场文件系列名称
% UV_file=dir([UV_folder_name,'SMSUV*.XY']);%旧版fvcom
UV_file=dir([UV_folder_name,'bohai_*.nc']);%新版fvcom
%各分区范围赋值
SubRegion_Range=textread(SubRegion);
SubRegion_Range_x=SubRegion_Range(:,1);
SubRegion_Range_y=SubRegion_Range(:,2);
SubRegion_Range_x(1:(SubRegion_Point_Number+1):end)=[];
SubRegion_Range_y(1:(SubRegion_Point_Number+1):end)=[];


%% 每天输出一个结果
% for filenum=1:6
    filenum=5;       %月份
    file=sprintf('%.3d',filenum)
    file_day=0;
% Particle_Output_file=fopen(['test\175lag_nondif.dat'],'w');   %质点测试所用
%% ---------------------------------------------------------------
%%%%%%    Particle Track Model      ---------------------------------------
for IsStep=StartTime:TimeStep_h:EndTime;
    IsStep
    ISStep=IsStep-StartTime+1;
%对流项
%  [X0,Y0,U,V,F1,F2,F3,F4]=textread([UV_folder_name,UV_file(IsStep).name], ...
%                       '%f %f %f %f %f %f %f %f','headerlines',2);%旧版fvcom
                  
 X0=double(ncread([UV_folder_name,UV_file(IsStep).name],'xc'));
 Y0=double(ncread([UV_folder_name,UV_file(IsStep).name],'yc'));
 uu=double(ncread([UV_folder_name,UV_file(IsStep).name],'u'));
 vv=double(ncread([UV_folder_name,UV_file(IsStep).name],'v'));
 U=uu(:,1);
 V=vv(:,1);%新版fvcom
 
 
 Particle_Now_U=griddata(X0,Y0,U,Particle_initial_x,Particle_initial_y);%粒子当前步的u
 Particle_Now_V=griddata(X0,Y0,V,Particle_initial_x,Particle_initial_y);%粒子当前步的v
 
 K1_Particle_X=Particle_Now_U.*TimeStep;%K1步走的距离
 K1_Particle_Y=Particle_Now_V.*TimeStep;%K1步走的距离
 
 Particle_K1_X=Particle_initial_x+0.5.*K1_Particle_X;%K2计算所需位置
 Particle_K1_Y=Particle_initial_y+0.5.*K1_Particle_Y;%K2计算所需位置
 
 Particle_K1_U=griddata(X0,Y0,U,Particle_K1_X,Particle_K1_Y);%粒子K2步计算所需u
 Particle_K1_V=griddata(X0,Y0,V,Particle_K1_X,Particle_K1_Y);%粒子K2步计算所需v 
 
 K2_Particle_X=Particle_K1_U.*TimeStep;%K2步走的距离
 K2_Particle_Y=Particle_K1_V.*TimeStep;%K2步走的距离
 
 Particle_K2_X=Particle_initial_x+0.5.*K2_Particle_X;%K3步计算所需位置
 Particle_K2_Y=Particle_initial_y+0.5.*K2_Particle_Y;%K3步计算所需位置
 
 Particle_K2_U=griddata(X0,Y0,U,Particle_K2_X,Particle_K2_Y);%粒子K3步计算所需u
 Particle_K2_V=griddata(X0,Y0,V,Particle_K2_X,Particle_K2_Y);%粒子K3步计算所需v 
 
 K3_Particle_X=Particle_K2_U.*TimeStep;%K3步走的距离
 K3_Particle_Y=Particle_K2_V.*TimeStep;%K3步走的距离
 
 Particle_K3_X=Particle_initial_x+0.5.*K3_Particle_X;%K4步计算所需位置
 Particle_K3_Y=Particle_initial_y+0.5.*K3_Particle_Y;%K4步计算所需位置
 
 Particle_K3_U=griddata(X0,Y0,U,Particle_K3_X,Particle_K3_Y);%粒子K4步计算所需u
 Particle_K3_V=griddata(X0,Y0,V,Particle_K3_X,Particle_K3_Y);%粒子K4步计算所需v 
 
 K4_Particle_X=Particle_K3_U.*TimeStep;%K4步走的距离
 K4_Particle_Y=Particle_K3_V.*TimeStep;%K4步走的距离
 
%扩散项
if(OpenDiffusion==1)
    dx=0.2;%扩散系数
    dy=0.2;%扩散系数
 Diff_Random_X=(rand(Particle_sum_number,1)-0.5).*2;
 Diff_Random_Y=(rand(Particle_sum_number,1)-0.5).*2 ;
 Diff_X=Diff_Random_X.*sqrt(dx.*2.*TimeStep); 
 Diff_Y=Diff_Random_Y.*sqrt(dy.*2.*TimeStep);
else
 Diff_X=0;
 Diff_Y=0;  
end
%总位移 
Particle_X=Particle_initial_x+(K1_Particle_X+2.*K2_Particle_X+2.*K3_Particle_X+K4_Particle_X)./6+Diff_X;
Particle_Y=Particle_initial_y+(K1_Particle_Y+2.*K2_Particle_Y+2.*K3_Particle_Y+K4_Particle_Y)./6+Diff_Y;

%[in,on]= inpolygon(Particle_X,Particle_Y,xv,yv); %岸反射
%index=find(in==1);
%Particle_initial_x(index)=Particle_X(index);
%Particle_initial_y(index)=Particle_Y(index);

Particle_initial_x=Particle_X;
Particle_initial_y=Particle_Y;

%结果输出，可用作后期检查使用
if(mod(ISStep,24)==0) %每隔一天输出一次粒子位置
    file_day=file_day+1;
    file_Day=sprintf('%.3d',file_day);
   Particle_Output_file=fopen(['5月\',file,'_',file_Day,'.dat'],'w');
  for out_particle=1:Particle_sum_number;
    fprintf(Particle_Output_file,'%10.0f %10.0f %20.7f %20.7f\n', ...
    ISStep,Particle_initial_SubRegion(out_particle,1), ...
     Particle_X(out_particle,1),Particle_Y(out_particle,1));
  end
  fclose(Particle_Output_file);
end
 end

% end


fprintf('Particle Track Model based on FVCOM end now!\n');
%%---------------------------------------------------------------

