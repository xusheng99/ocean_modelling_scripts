clc
clear

%数据源
ftx = load('cod3.txt');


%作图点间隔，可选0.05，0.1
interval = 0.05


load('bohai_coastline_i.mat')


lon3 = ftx(:,1);
lat3 = ftx(:,2);
value3 = ftx(:,3);




x = 117.2:interval:122.5;
y = 37:interval:41;
[xx,yy] = meshgrid(x,y);
z3 = []
z8 = []



plot(ncst(:,1),ncst(:,2))
axis([117.2,122.5,37,41])
axis equal
hold on


F3 = scatteredInterpolant(lon3,lat3,value3);


for loni = 117.2:interval:122.5
    for latj = 37:interval:41
        in = inpolygon(loni,latj,ncst(:,1),ncst(:,2));
        i = ((loni-117.2)/interval)+1
        i = int8(i)
        j = ((latj-37)/interval)+1
        j = int8(j)
        if in == 1
           z3(i,j) = NaN;
        end
        if in == 0
            z3(i,j) = F3(loni,latj);
        end
    end
end
zz  = z3
zz = zz'

%此处-2 2 10分别代表作图时colorbar上想要显示的最小值、最大值、和分割份数 
%最大最小值可以先手动看一下zz这个数组，看看里面的数字大概都是什么数量范围的
levels = linspace(-2,2,16);
contourf(xx,yy,zz,levels)
axis equal
colorbar