from netCDF4 import Dataset
import numpy as np
import os

workdir = 'C:\\Users\\XUSHENG\\Desktop\\read_roms_result\\'    #将某些较小的nc文件直接放到这里
plotdir = 'C:\\Users\\XUSHENG\\Desktop\\作图用文件\\'    #读取的数据移动到这里供matlab使用
os.chdir(workdir)

lon1 = 110
lat1 = 20.5
lon2 = 115.7
lat2 = 23.2


hl = int(((lon1+180)/360)*10800)
hr = int(((lon2+180)/360)*10800)
hd = int(((lat1+90)/180)*5400)
hu = int(((lat2+90)/180)*5400)

h = Dataset('etopo1.nc')
h = h.variables['topo'][:,:]
hthumb = h[0:5400:10,0:10800:10]
h_scs = h[hd:hu,hl:hr]
np.savetxt('etopothumbnail.txt',hthumb)
# np.savetxt('etopo.txt',h)
np.savetxt('etoposcs.txt',h_scs)


