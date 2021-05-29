import netCDF4 as nc
import os
import numpy as np 
import math
import shutil
import matplotlib.path as mplPath

import platform

system = platform.system()
if system == 'Windows':
    file_domain = 'D:/application/SCS/DATA/Grd.nc'
    # coast = 'C:\\Users\\XUSHENG\\Documents\\scripts\\paticle_tracking_2dim_for_ROMS\\llbounds.bln'
    coast = 'D:\\application\\PLANT\\llbounds.bln'
    # flowfield_prefix = 'c:\\Users\\XUSHENG\\Desktop\\roms\\'
# if system == 'Linux':
#     file_domain = '/mnt/d/application/SCS/DATA/Grd.nc'
#     coast = '/mnt/c/Users/XUSHENG/Documents/scripts/particle_tracking_2dim_for_ROMS/llbounds.bln'

## 读取网格大小信息
# file_domain = 'D:/roms_ini.nc'
content = nc.Dataset(file_domain)
__mask_rho = content.variables['mask_rho'][:].data
lon_rho = content.variables['lon_rho'][:].data
lat_rho = content.variables['lat_rho'][:].data
lonmin, lonmax, latmin, latmax = np.min(lon_rho), np.max(lon_rho), np.min(lat_rho), np.max(lat_rho)
lonarray,latarray = lon_rho[0,:],lat_rho[:,0]
len_xi,len_eta = len(lonarray),len(latarray) 

# pnts = [[106,24],[116,15],[110,19]] # 分别是大陆, 海, 台湾中心
# pnts = [[110,19]]

a = np.loadtxt(coast,delimiter=',')
# 找出分割符所在行
b = a==1
# 读取其下标
c = np.where(b == True)[0]
# 找到每个多边形的边数
dis = np.array([c[i+1]-c[i] for i in range(len(c)-1)]) - 1
d = np.zeros([len(c)-1,max(dis),2]) # 多边形编号, 节点编号, 经纬
# 提取每个多边形的节点, 存储在不同的行, 不够的位数补零
for i in range(d.shape[0]):
    d[i,0:dis[i],:] = a[c[i]+1:c[i+1],:]
# 创建多个多边形, 用于补齐的零是不计入的
polys = []
for i in range(d.shape[0]): # 创建多个多边形
    # i= 1
    crd = d[i,0:dis[i]-1,:]
    poly = mplPath.Path(crd)
    polys.append(poly)

# pnt.shape = (2,2)  && poly.shape=(n,2)
def inpoly(pnt,poly):
    poly = mplPath.Path(poly)
    r = 0
    return (poly.contains_point(pnt,radius=r) or poly.contains_point(pnt,radius=-r))

# 使用LTRANS构建的surfer多边形来判定的方法, 更准确
def inocean(pnt):
    isIns = []
    r = 0
    for i in range(d.shape[0]): # 对每个多边形判定
        isIn = [polys[i].contains_point(pnt,radius=r) or polys[i].contains_point(pnt,radius=-r)]
        isIns.append(isIn)
    isIns = np.squeeze(np.array(isIns)) # 删除多余维度

    pi = False
    for element in isIns[1:]: # 所有小多边形, 必须全False才返回False
        pi = pi + element
    if isIns[0] and (not pi): # 必须要在大多边形内, 且不在其余的任何一个小多边形内, 才为海点
        return True
    else:
        return False

# 使用掩码来判定海陆的方法
def inocean2(cordinate):
    lon = cordinate[0]
    lat = cordinate[1]
    xinterval = (len_xi - 1)/(lonarray[-1] - lonarray[0])
    yinterval = (len_eta - 1)/(latarray[-1] - latarray[0])
    x = (lon - lonarray[0])*xinterval
    y = (lat - latarray[0])*yinterval
    def is_edge_or_boundary(x1,y1):
        return [(__mask_rho[y1,x1] == 1 ), (__mask_rho[y1,x1+1] == 1), (__mask_rho[y1+1,x1+1] == 1), (__mask_rho[y1+1,x1] == 1)]
    x = int(x//1) 
    y = int(y//1)
    inocean = is_edge_or_boundary(x,y)
    return inocean[0]+inocean[1]+inocean[2]+inocean[3]

def ncread(ncfile,varname):
    dataset = nc.Dataset(ncfile)
    a = dataset.variables[varname][:]
    return a

def ncwrite(matrix,varname,ncfile):
    content = nc.Dataset(ncfile,'r+') #务必r+
    content.variables[varname][:] = matrix
    print(ncfile+'\'s variable:',varname,'\n','has been written!')
    content.close()