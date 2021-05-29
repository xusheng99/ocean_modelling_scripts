clear
clc

%----指定参数----
%请指定timeinterval和endtime_hour，注意endtime_hour不能超过timelength
%(endtime_hour-starttime_hour)必须要能整除timeinterval
timeinterval =1
starttime_hour = 1
endtime_hour = 300
%%%%%%%%%%%%%%%%%%%


addpath C:\Users\XUSHENG\Desktop\作图用文件

addpath d:\application\SCS\
hisname = 'his_0003.nc'


u=load('[36,168]_u_his_point.txt');
v=load('[36,168]_v_his_point.txt');
size = size(u);
timelength = size(1);
if endtime_hour <= timelength
    u22=u(starttime_hour+1:timeinterval:endtime_hour,1);
    v22=v(starttime_hour+1:timeinterval:endtime_hour,1);
else
    disp('endtime_hour设定过大，请检查！')
end  
temp = (endtime_hour-starttime_hour)/timeinterval;
for j = 1:temp
    x(j)=starttime_hour+timeinterval*j;
    if v22(j,1) > 0
        y(j) = 0;
    else 
        y(j) = 0;
    end
end
x=x';
y=y';
quiver(x,y,u22,v22)
axis equal
ploth = floor((endtime_hour-starttime_hour)/20);
axis([starttime_hour,endtime_hour,-ploth*2,ploth*2]);