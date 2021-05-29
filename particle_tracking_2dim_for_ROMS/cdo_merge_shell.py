# cdo生成目标时间段的surface_uv.nc4的命令
import os
os.chdir('d:/application/PLANT/')

def hisname(i):
    return 'his_'+"{:0>4d}".format(i)+'.nc'

start = 20
stride = 90

with open('cdo_merge.sh','w',encoding='utf-8') as f:
    f.write('# /usr/bin/bash\n')
    f.write('# extract surface flowfield\n')
    f.write('cdo -select,name=u_eastward,v_northward,level=-0.03125 ')
    for i in range(start,start+stride):
        f.write(str(hisname(i))+' ')
    f.write('surface_uv.nc4')
f.close()


