# 构建几个若干个区域的中心点
# 根据长方形边长原则给出每个点在哪个区域的函数
# 对始末状态, 由经纬度转化成区域
# 统计初始状态每个区域的数目, 归一化处理
# 构建邻接矩阵, 维度很小


import numpy as np
import matplotlib.path as mplPath
import os

def inpoly(pnt,poly):
    poly = mplPath.Path(poly)
    r = 0
    return (poly.contains_point(pnt,radius=r) or poly.contains_point(pnt,radius=-r))

poly1 = np.array([[111,20.5],[111,23],[112,23],[112,20.7]])
poly2 = np.array([[112,20.7],[112,23],[113,23],[113,20.9]])
poly3= np.array([[113,20.9],[113,23],[114,23],[114,21.1]])
poly4 = np.array([[114,21.1],[114,23],[115,23],[115,21.3]])
poly5 = np.array([[115,21.3],[115,23],[116,23],[116,21.5]])
polys = np.zeros([poly1.shape[0],poly1.shape[1],5])


def which_poly(pnt):
    if inpoly(pnt,poly1):
        return int(1)
    elif inpoly(pnt,poly2):
        return int(2)
    elif inpoly(pnt,poly3):
        return int(3)
    elif inpoly(pnt,poly4):
        return int(4)
    elif inpoly(pnt,poly5):
        return int(5)
    else:
        return 0
    
os.chdir('C:\\users\\xusheng\\desktop\\roms\\')

# 载入初始时刻的位置和结束时刻的位置, 经纬度格式表示
ini = np.loadtxt('initial.dat',delimiter=' ')  # 不同的粒子追踪程序输出的格式不同, 使用comma或space作为间隔.
final = np.loadtxt('final.dat',delimiter=' ')

## 清洗数据: 删除离开了开边界的
ini = [ini[i] for i in range(ini.shape[0]) if 0.1<final[i,0]<1000]
final = [final[i] for i in range(final.shape[0]) if 0.1<final[i,0]<1000]

## 清洗数据: 删除在陆地上的点
ini2 = [[ini[i][0],ini[i][1]] for i in range(len(ini)) if ini[i][0] != final[i][0] and ini[i][1] != final[i][1]]
final2 = [[final[i][0],final[i][1]] for i in range(len(final)) if ini[i][0] != final[i][0] and ini[i][1] != final[i][1]]

ini_grid = np.zeros(len(ini2),dtype='int')
final_grid = np.zeros(len(final2),dtype='int')
for i in range(len(ini2)):
    ini_grid[i] = which_poly(ini2[i])
    final_grid[i] = which_poly(final2[i])

np.savetxt('')

ini_grid2 = [ini_grid[i] for i in range(len(ini_grid)) if 1<=ini_grid[i]<=5]
final_grid2 = [final_grid[i] for i in range(len(ini_grid)) if 1<=ini_grid[i]<=5]

length = 6
A = np.zeros([length,length])
for i in range(len(ini_grid2)):
    A[ini_grid2[i],final_grid2[i]] += 1
A