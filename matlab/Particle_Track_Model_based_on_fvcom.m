    %% --------Particle Track based on the FVCOM Flow Filed----------
%%%%%%%%%%%%%%%%%%%%%%�����Ľ���������㷨%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       Yi+1 = Yi +1/6*K1 +1/3*K2+1/3*K3+1/6*K4;                      %%% 
%%%       K1 = hf(Xi,Yi);                                               %%%
%%%       K2 = hf(Xi+h/2, Yi+K1/2);                                     %%%
%%%       K3 = hf(Xi+h/2, Yi+K2/2);                                     %%%
%%%       K4 = hf(Xi+h, Yi+K3);                                         %%%
%%%       ׷��ģ�������������ɢ                                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   �ļ���ʽ˵��   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       ��Ҫ�������ļ��У��ͷ����ӵ�ԭʼλ���ļ���lag_ini.dat����         %%%
%%%       ��fvcom����׷�ٳ�ʼ�����ļ���ʽ��ͬ,                             %%%

%%%       ���ݶ�ȡ�͸�ֵʱ�Ѿ���������ֵɾ��                               %%%
%%%       fcvom ���������ļ���ֱ�Ӷ�ȡϵ���ļ���                           %%%
%%%       �����������꣬��ʽΪsurfer�ɻ�������bln�ļ���                    %%%
%%%       ��Ҫ������ļ��У��ͷ����ӵĵ���λ���ļ���lag_out.dat����         %%%
%%%       ��ʽΪÿ�����ӵ���λ�õ�ʱ��ѭ��������ѭ�����Ӹ�����             %%%
%%%       ���Ӹ�ʱ�䵽��ķ���������������Դ�������                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;
fprintf('Particle Track Model based on FVCOM Start Now!\n');
%% ---------------------------------------------------------------
%�����ļ�
% Particle_initial_file='2018���߽��\�����������з���\bohai_lag_ini_jy200m.dat';%fvcom����׷�ٵ������ļ�
Particle_initial_file='��ʼλ��2500.dat';%fvcom����׷�ٵ������ļ�
% Particle_initial_file='2018���߽���޷�\�����������з���\�ͷ�����·���ļ�\lag_nondif001_243.dat';%fvcom����׷�ٵ������ļ�%��Ҫ������
% UV_folder_name='jingliu_UV\';%�ɰ�fvcom����������Ҫ����uv
UV_folder_name='2018UV365day\';%�°�fvcom����������Ҫ����uv
SubRegion='bohaifenqu.bln';%����surfer�л��������������ļ�
A = importdata('2018bhx.dat');
xv=A(:,1);
yv=A(:,2);
%����ļ�
% Particle_Output_file=fopen('lag_results\lag_dif_out.dat','w');
%% ---------------------------------------------------------------
%��������
SubRegion_Number=751;%������������������������ݵ�ϵ������֮�䲻ͬ���������ļ�����
Particle_SubRegion_Number=2500.*ones(SubRegion_Number,1);%��ʼʱÿ��������������
SubRegion_Point_Number=5.*ones(SubRegion_Number,1);%һ����������˵���������βһ��
TimeStep=3600;%ʱ�䲽������λΪs
TimeStep_h=TimeStep/3600;%ʱ�䲽������λΪh�������ṩʱ��ѭ������
% StartTime=165;%�ͷſ�ʼʱ�䣬��λΪh
StartTime=2879;%+90*24-1;%�ͷſ�ʼʱ�䣬��λΪh,���������Ϊ5�����Ѿ����棬Ҫ���������Ϊ������׷�٣�����������
EndTime  =3240;%׷�ٽ���ʱ�䣬��λΪh
OpenDiffusion=1;%��ɢ������ƣ�1Ϊ������0Ϊ�ر�

%% ---------------------------------------------------------------
%���ӳ�ʼλ�ø�ֵ
Particle_initial=textread(Particle_initial_file);
Particle_initial_x=Particle_initial(1:end,1);
Particle_initial_y=Particle_initial(1:end,2);
% Particle_initial_x=Particle_initial(:,3);%����ʱ
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
%�ٶȳ��ļ�ϵ������
% UV_file=dir([UV_folder_name,'SMSUV*.XY']);%�ɰ�fvcom
UV_file=dir([UV_folder_name,'bohai_*.nc']);%�°�fvcom
%��������Χ��ֵ
SubRegion_Range=textread(SubRegion);
SubRegion_Range_x=SubRegion_Range(:,1);
SubRegion_Range_y=SubRegion_Range(:,2);
SubRegion_Range_x(1:(SubRegion_Point_Number+1):end)=[];
SubRegion_Range_y(1:(SubRegion_Point_Number+1):end)=[];


%% ÿ�����һ�����
% for filenum=1:6
    filenum=5;       %�·�
    file=sprintf('%.3d',filenum)
    file_day=0;
% Particle_Output_file=fopen(['test\175lag_nondif.dat'],'w');   %�ʵ��������
%% ---------------------------------------------------------------
%%%%%%    Particle Track Model      ---------------------------------------
for IsStep=StartTime:TimeStep_h:EndTime;
    IsStep
    ISStep=IsStep-StartTime+1;
%������
%  [X0,Y0,U,V,F1,F2,F3,F4]=textread([UV_folder_name,UV_file(IsStep).name], ...
%                       '%f %f %f %f %f %f %f %f','headerlines',2);%�ɰ�fvcom
                  
 X0=double(ncread([UV_folder_name,UV_file(IsStep).name],'xc'));
 Y0=double(ncread([UV_folder_name,UV_file(IsStep).name],'yc'));
 uu=double(ncread([UV_folder_name,UV_file(IsStep).name],'u'));
 vv=double(ncread([UV_folder_name,UV_file(IsStep).name],'v'));
 U=uu(:,1);
 V=vv(:,1);%�°�fvcom
 
 
 Particle_Now_U=griddata(X0,Y0,U,Particle_initial_x,Particle_initial_y);%���ӵ�ǰ����u
 Particle_Now_V=griddata(X0,Y0,V,Particle_initial_x,Particle_initial_y);%���ӵ�ǰ����v
 
 K1_Particle_X=Particle_Now_U.*TimeStep;%K1���ߵľ���
 K1_Particle_Y=Particle_Now_V.*TimeStep;%K1���ߵľ���
 
 Particle_K1_X=Particle_initial_x+0.5.*K1_Particle_X;%K2��������λ��
 Particle_K1_Y=Particle_initial_y+0.5.*K1_Particle_Y;%K2��������λ��
 
 Particle_K1_U=griddata(X0,Y0,U,Particle_K1_X,Particle_K1_Y);%����K2����������u
 Particle_K1_V=griddata(X0,Y0,V,Particle_K1_X,Particle_K1_Y);%����K2����������v 
 
 K2_Particle_X=Particle_K1_U.*TimeStep;%K2���ߵľ���
 K2_Particle_Y=Particle_K1_V.*TimeStep;%K2���ߵľ���
 
 Particle_K2_X=Particle_initial_x+0.5.*K2_Particle_X;%K3����������λ��
 Particle_K2_Y=Particle_initial_y+0.5.*K2_Particle_Y;%K3����������λ��
 
 Particle_K2_U=griddata(X0,Y0,U,Particle_K2_X,Particle_K2_Y);%����K3����������u
 Particle_K2_V=griddata(X0,Y0,V,Particle_K2_X,Particle_K2_Y);%����K3����������v 
 
 K3_Particle_X=Particle_K2_U.*TimeStep;%K3���ߵľ���
 K3_Particle_Y=Particle_K2_V.*TimeStep;%K3���ߵľ���
 
 Particle_K3_X=Particle_initial_x+0.5.*K3_Particle_X;%K4����������λ��
 Particle_K3_Y=Particle_initial_y+0.5.*K3_Particle_Y;%K4����������λ��
 
 Particle_K3_U=griddata(X0,Y0,U,Particle_K3_X,Particle_K3_Y);%����K4����������u
 Particle_K3_V=griddata(X0,Y0,V,Particle_K3_X,Particle_K3_Y);%����K4����������v 
 
 K4_Particle_X=Particle_K3_U.*TimeStep;%K4���ߵľ���
 K4_Particle_Y=Particle_K3_V.*TimeStep;%K4���ߵľ���
 
%��ɢ��
if(OpenDiffusion==1)
    dx=0.2;%��ɢϵ��
    dy=0.2;%��ɢϵ��
 Diff_Random_X=(rand(Particle_sum_number,1)-0.5).*2;
 Diff_Random_Y=(rand(Particle_sum_number,1)-0.5).*2 ;
 Diff_X=Diff_Random_X.*sqrt(dx.*2.*TimeStep); 
 Diff_Y=Diff_Random_Y.*sqrt(dy.*2.*TimeStep);
else
 Diff_X=0;
 Diff_Y=0;  
end
%��λ�� 
Particle_X=Particle_initial_x+(K1_Particle_X+2.*K2_Particle_X+2.*K3_Particle_X+K4_Particle_X)./6+Diff_X;
Particle_Y=Particle_initial_y+(K1_Particle_Y+2.*K2_Particle_Y+2.*K3_Particle_Y+K4_Particle_Y)./6+Diff_Y;

%[in,on]= inpolygon(Particle_X,Particle_Y,xv,yv); %������
%index=find(in==1);
%Particle_initial_x(index)=Particle_X(index);
%Particle_initial_y(index)=Particle_Y(index);

Particle_initial_x=Particle_X;
Particle_initial_y=Particle_Y;

%�����������������ڼ��ʹ��
if(mod(ISStep,24)==0) %ÿ��һ�����һ������λ��
    file_day=file_day+1;
    file_Day=sprintf('%.3d',file_day);
   Particle_Output_file=fopen(['5��\',file,'_',file_Day,'.dat'],'w');
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

