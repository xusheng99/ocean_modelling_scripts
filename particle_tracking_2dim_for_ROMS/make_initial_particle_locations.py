import numpy as np
import os
from myfunc.nc_post_process import inocean
from myfunc.nc_post_process import inocean2


os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\')

minX = 105.5
maxX = 124.5

minY = 5.5
maxY = 24.5


x = np.linspace(minX,maxX,100)
y = np.linspace(minY,maxY,100)

inocean([124.3080808080808,13.56060606060606])

# for LTRANS, with comma
with open('Initial_particle_locations.csv','w',encoding='utf-8') as f:
    for eachlon in x:
        for eachlat in y:
            if inocean([eachlon,eachlat]) and inocean2([eachlon,eachlat]):
                f.write(str(eachlon))
                f.write(',')
                f.write(str(eachlat))
                f.write(',')
                f.write('-3.0')  #深度
                f.write(',')
                f.write('3600')  #date of birth
                f.write('\n')


# # for my own python script, without comma
# with open('particles_initial.txt','w',encoding='utf-8') as f:
#     for eachlon in x:
#         for eachlat in y:
#             f.write(str(eachlon)+' ')
#             f.write(str(eachlat)+' ')
#             f.write('3600'+' ')  #date of birth
#             f.write('\n')