clc
clear

%数据源
cod3 = load('cod3.txt');
cod8 = load('cod8.txt');

%作图点间隔，可选0.05，0.1
interval = 0.05

%尾缀 l i f zi'zuo
load('bohai_coastline_f.mat')


lon3 = cod3(:,1);
lat3 = cod3(:,2);
value3 = cod3(:,3);


lon8 = cod8(:,1);
lat8 = cod8(:,2);
value8 = cod8(:,3);

x = 117.2:interval:122.5;
y = 37:interval:41;
[xx,yy] = meshgrid(x,y);
z3 = [];
z8 = [];



plot(ncst(:,1),ncst(:,2))
axis([117.2,122.5,37,41])
axis equal
hold on


F3 = scatteredInterpolant(lon3,lat3,value3);
F8 = scatteredInterpolant(lon8,lat8,value8);

for loni = 117.2:interval:122.5
    for latj = 37:interval:41
        in = inpolygon(loni,latj,ncst(:,1),ncst(:,2));
        i = ((loni-117.2)/interval)+1;
        i = int8(i);
        j = ((latj-37)/interval)+1;
        j = int8(j);
        if in == 1
           z3(i,j) = NaN;
           z8(i,j) = NaN;
        end
        if in == 0
            z3(i,j) = F3(loni,latj);
            z8(i,j) = F8(loni,latj);

        end
    end
end
zz = z8 - z3;
zz = zz';

        
levels = linspace(-2,2,10);
contourf(xx,yy,zz,levels)
axis equal
colorbar
