import numpy as np
import netCDF4 as nc
import matplotlib.pylab as plt
import sys
import os

workdir = sys.argv[1]
# workdir = 'C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\plant_tdn\\20'
os.chdir(workdir)

# N = sys.argv[1]

## 读取网格大小信息
file_domain = 'D:/application/PLANT/his_0001.nc'
content = nc.Dataset(file_domain)
mask_rho = content.variables['mask_rho'][:].data
lon_rho = content.variables['lon_rho'][:].data
lat_rho = content.variables['lat_rho'][:].data

# ## 看看
# plt.pcolor(mask_rho)
# plt.show()

## 提取海点 1/12 degrees --> 1/4 degrees
N = 20
length1 = (mask_rho.shape[0]-1)//N+1
length2 = (mask_rho.shape[1]-1)//N+1
maskr = mask_rho[::N,::N].reshape([length1*length2,1])
lonr = lon_rho[::N,::N].reshape([length1*length2,1])
latr = lat_rho[::N,::N].reshape([length1*length2,1])
points_inocean = [[lonr[i][0],latr[i][0]] for i in range(len(maskr)) if maskr[i] == 1]

with open('grid_center_inocean.dat','w',encoding='utf-8') as f:
    count = 0
    for line in points_inocean:
        f.write(str(line[0])+'    '+str(line[1])+'    '+str(count)+'\n')
        count += 1
f.close()
























# 顺带统计下有多少个点在海里
# count = 0
# for i in np.linspace(105.5,124.5,120):
#     for j in np.linspace(5.5,24.5,120):
#         if inocean([i,j]):
#             count += 1
