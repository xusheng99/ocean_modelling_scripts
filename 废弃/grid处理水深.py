dimport netCDF4 as nc4
import numpy as np 
import os 


#指定grid文件所在路径
A='C:\\Users\\XUSHENG\\Desktop\\'
gridfilename='roms_grid.nc'

os.chdir(A)
#不要忘记以r+模式打开，否则会写入失败而且不报错
gridfile=nc4.Dataset(gridfilename,'r+')
h=gridfile.variables['h'][:]

i,j = h.shape


'''
print(i,j,'is the shape of h')
for ii in range(i):
    for jj in range(j):
        if h[i,j] < 10:
            print('find',i,j,h[i,j])
            h[i,j]=10
'''

h[h<10]=10
print(h.max(),h.min())

gridfile.variables['h'][:]=h[:]
gridfile.close()

# check 
gridfile=nc4.Dataset(gridfilename)
h=gridfile.variables['h'][:]
if h.min() > 0:
    print('replacement is done')
else:
    print('something went wrong')




