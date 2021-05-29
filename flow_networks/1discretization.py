import os
import numpy as np
import matplotlib.pyplot as plt
import sys

workdir = sys.argv[1]
os.chdir(workdir)
# os.chdir('C:\\Users\\XUSHENG\\Desktop\\scripts\\flow_networks\\')
p = np.loadtxt('grid_center_inocean.dat',dtype=type(9.9))

N = 20
edge = N/60 # 网格大小: 1/3度

# 输入经纬度, 返回格点坐标, 如果输入的点在海里, 那么返回nan
def grid_num(cordinate):
    lon = cordinate[0]
    lat = cordinate[1]
    for i in range(len(p)):
        # 此处的左闭右开判定方法导致了出度图奇怪的"格栅化"
        if (p[i,0]-edge*0.5<=lon<p[i,0]+edge*0.5) and (p[i,1] -edge*0.5<=lat<p[i,1]+edge*0.5): 
            return i
    return np.NaN

## 载入初始时刻的位置和结束时刻的位置, 经纬度格式表示
ini = np.loadtxt('initial.dat',delimiter=' ')  # 不同的粒子追踪程序输出的格式不同, 使用comma或space作为间隔.
final = np.loadtxt('final.dat',delimiter=' ')

## 清洗数据: 删除离开了开边界的
ini = [ini[i] for i in range(ini.shape[0]) if 0.1<final[i,0]<1000]
final = [final[i] for i in range(final.shape[0]) if 0.1<final[i,0]<1000]

## 清洗数据: 删除在陆地上的点, 注意, 一朝为陆, 终生为陆, ini和final都是陆地的
ini2 = [[ini[i][0],ini[i][1]] for i in range(len(ini)) if ini[i][0] != final[i][0] and ini[i][1] != final[i][1]]
final2 = [[final[i][0],final[i][1]] for i in range(len(final)) if ini[i][0] != final[i][0] and ini[i][1] != final[i][1]]

## 保存剩余的有效的点的始末位置
np.savetxt('initial_valid.dat',ini2)
np.savetxt('final_valid.dat',final2)

## 对剩余的有效点, 执行网格化操作
ini_grid = []
final_grid = []
for each_particle in ini2:
    ini_grid.append(grid_num(each_particle))
for each_particle in final2:
    final_grid.append(grid_num(each_particle))
ini_gridd = [ini_grid[i] for i in range(len(ini_grid)) if type(ini_grid[i]) == type(9) and type(final_grid[i]) == type(9)]
final_gridd = [final_grid[i] for i in range(len(final_grid)) if type(ini_grid[i]) == type(9) and type(final_grid[i]) == type(9)]

## 按两列写入数据
with open('Ini2Final.dat','w',encoding='utf-8') as f:
    for i in range(len(ini_gridd)):
        f.write(str(ini_gridd[i])+' '+str(final_gridd[i])+'\n')
f.close()   


