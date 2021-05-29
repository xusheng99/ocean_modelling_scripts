import os
import numpy as np

os.chdir('c:/users/xusheng/documents/scripts/flow_networks/')

path =[1479, 1525, 1572, 1665, 1669, 1763, 1910, 2007, 2004, 1950, 1848, 1799, 1902]
p = np.loadtxt('grid_center_inocean.dat')

for i in range(len(p)):
    p[i,2] = int(p[i,2])

def cor(num):
    return [p[num,0],p[num,1]]

with open('path.dat','w',encoding='utf-8') as f:
    for node in path:
        f.write(str(cor(node)[0])+'    '+str(cor(node)[1])+'\n')
f.close()


