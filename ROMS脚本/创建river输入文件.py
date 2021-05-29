from netCDF4 import Dataset
import numpy as np 
import os 

def ncread(ncfile,varname):
    dataset = Dataset(ncfile)
    a = dataset.variables[varname][:]
    return a

#默认river文件
riverfile = 'river.nc'
def ncwrite(matrix,varname,ncfile=riverfile):
    content = Dataset(ncfile,'r+') #务必r+
    content.variables[varname][:] = matrix
    print(ncfile+'\'s variable:',varname,'\n','has been written!')
    content.close()

os.chdir('c:\\Users\XUSHENG\\Desktop\\river_input\\')

time = [
15.2188, 
45.6563, 
76.0938,
106.5313, 
136.9688, 
167.4063, 
197.8438, 
228.2813, 
258.7188, 
289.1563, 
319.5938, 
350.0313
]

ncwrite(np.array([0.0,1.0]),'river') #河流标识符
ncwrite(time,'river_time') #河流时间维度

ncwrite([101,15],'river_Xposition') #格点位置
ncwrite([214,53],'river_Eposition')

#径流
runoff = [[2501.46163682864,
2683.41317126574,
3638.02736570646,
7180.29411764706,
13331.6988308367,
19037.1849104859,
18080.5498355864,
16656.65991962,
11118.7714285714,
6559.19682539683,
4634.75555555556,
3027.73333333333],
[2793.12009634623,
2149.37459359518,
1807.73249333079,
1786.67235151544,
2922.55966149755,
8942.78921545557,
17380.7837208208,
27390.9081089985,
27577.3344785392,
16108.3613045934,
8031.46569677146,
4265.9013008777]]
runoff_ = np.zeros((12,2))
runoff_[:,0] = runoff[0]
runoff_[:,1] = runoff[1]
ncwrite(runoff_,'river_transport')
a = ncread(riverfile,'river_transport')

river_Vshape = np.zeros((40,2)) #垂向流量分布, 直接弄成简单的均匀分布
river_Vshape[:,0] = np.linspace(0.025,0.025,40)
river_Vshape[:,1] = np.linspace(0.025,0.025,40)
ncwrite(river_Vshape,'river_Vshape')

river_salt = np.ones((12,40,2)) #盐度设成1
ncwrite(river_salt,'river_salt')
