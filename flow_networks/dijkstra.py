'''二维Dijkstra算法'''

import os
import numpy as np
import math

# os.chdir('c:/users/xusheng/desktop/scripts/flow_networks/')
os.chdir('c:/users/xusheng/documents/scripts/flow_networks/')

# for i in range(len(mgraph)):
#     print(sum(mgraph[i,:]))

lon1 = 115
lat1 = 15
lon2 = 115
lat2 = 18

p = np.loadtxt('grid_center_inocean.dat',dtype=type(9.9))

N = 4
edge = N/12 # 网格大小: 0.5度

inf = 10086      # 给无限大定义一个值

mgraph = np.load('Adjacency_Matrix.npy')
mgraph = mgraph / 307.78

# dijstra 不能处理负边权问题, 务必选取合适的底数,使得取到的对数均为正数
for i in range(len(mgraph)):
    for j in range(len(mgraph)):
        if mgraph[i,j] != 0:
            mgraph[i,j] = -1*math.log(mgraph[i,j])

mgraph[mgraph==0] = inf

# 输入经纬度, 返回格点坐标, 如果输入的点在海里, 那么返回nan
def grid_num(cordinate):
    lon = cordinate[0]
    lat = cordinate[1]
    for i in range(len(p)):
        if (p[i,0]-edge*0.5<=lon<p[i,0]+edge*0.5) and (p[i,1] -edge*0.5<=lat<p[i,1]+edge*0.5): 
            return i
    return np.NaN

def dijkstra(start,end):
    # # 始末位置网格编号（从1开始计）

    # mgraph = [[0, 1, 12, inf, inf, inf],
    #           [inf, 0, 9, 3, inf, inf],
    #           [inf, inf, 0, inf, 5, inf],
    #           [inf, inf, 4, 0, 13, 15],
    #           [inf, inf, inf, inf, 0, 4],
    #           [inf, inf, inf, inf, inf, 0]]
    # mgraph = np.array(mgraph)

    passed = [start-1]
    nopass = [x for x in range(len(mgraph)) if x != start-1]
    dis = mgraph[start-1]

    path = [[start-1] for _ in range(len(mgraph))]
    for n in range(len(mgraph)):
        if dis[n] < inf:
            path[n].append(n)
    # print(path)

    while len(nopass) != 0:
        # 找到nopass中距离最短的点
        idx = nopass[0]
        for i in nopass:
            if dis[i] < dis[idx]:
                idx = i
        # 找到后把这个点移到passed中
        nopass.remove(idx)
        passed.append(idx)
        # 根据找到的这个点，更新最短路径
        for i in nopass:
            if dis[idx] + mgraph[idx][i] < dis[i]:
                dis[i] = dis[idx] + mgraph[idx][i]   # 更新距离
                path[i] = path[idx] + [i]            # 更新路径n
    # print(dis)
    # print(path)
    return path[end-1]

print(dijkstra(grid_num([lon1,lat1]),grid_num([lon2,lat2])))
# print('{} to {}:\ndis={}\npath={}'.format(start, end, dis[end-1], path[end-1]))
