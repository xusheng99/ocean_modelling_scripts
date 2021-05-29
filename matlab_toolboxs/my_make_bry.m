clc; clear;

% 1. 做出给定位置(lon,lat,depth)返回值的插值函数.
% 2. 根据zlevs函数给出边界上二维平面的所有网格点的lon_rho(lon_u),lat_rho(lat_u)以及sigma层--zeta+h-->depth
% 3. 根据2中给出的网格点的坐标, 给出当地的数值, 并写到合适大小的数组的正确的位置上.
% 4. 使用ncwrite将3中生成的数组装入nc文件这一容器中.nc文件使用"ncgen -b bry.cdl"生成.

%his_plot中用的zlevs插值方法是正确的,和roms_gui生成的1000m深处的温盐图进行了对比,各种细节纹理都是一致的.







