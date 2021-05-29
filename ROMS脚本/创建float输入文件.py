import os
import numpy as np
import re
import shutil

os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\')

def move(sourcefile,targetdir):    #移动sourcefile至targetdir,如果已存在同名文件则覆盖
    pwd = os.getcwd()
    os.chdir(targetdir)
    filelist = os.listdir(targetdir)
    if sourcefile in filelist:
        os.remove(sourcefile)
    os.chdir(pwd)
    shutil.move(sourcefile,targetdir)

with open('float_roms_prefix.txt','r',encoding='utf-8') as f:
    prefix = f.read()
f.close()

with open('float_roms_suffix.txt','r',encoding='utf-8') as f:
    suffix = f.read()
f.close()

minX = 106
maxX = 124
minY = 6
maxY = 24
nn = 200
x = np.linspace(minX,maxX,nn)
y = np.linspace(minY,maxY,nn)

info = []
ft0 = '0.0d0'
z0 = '-1.d0'
for i in range(nn):
    for j in range(nn):
        info.append(['1','1','1','1',ft0,str(x[i]),str(y[j]),z0,'0.d0','0.d0','0.d0','0.d0'])
# .replace('.','.d')

with open('floats.in','w',encoding='utf-8') as f:
    f.write(prefix)
    f.write('\n')
    for eachline in info:
        f.write('\n')
        for eachelement in eachline:
            f.write(eachelement)
            f.write('    ')
    f.write('\n')
    f.write('\n')
    f.write(suffix)

move('floats.in','D:\\application\\SCS\\')