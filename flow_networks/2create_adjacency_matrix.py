#### 说明 ####
# 将始末位置的格点编号转化为邻接矩阵
##############

import os
import numpy as np
import matplotlib.pyplot as plt

import sys

workdir = sys.argv[1]
os.chdir(workdir)

# os.chdir('C:\\Users\\XUSHENG\\Desktop\\scripts\\flow_networks\\')

## 读取邻接矩阵长度
grid_center = np.loadtxt('grid_center_inocean.dat')
length = len(grid_center)

## 导入对应的始末格点编号
data = np.loadtxt('Ini2Final.dat',dtype=type(21))
ini = data[:,0]
final = data[:,1]

## 创建合适大小的邻接矩阵, 遍历所有粒子, 并向对应的位置加值
A = np.zeros([length,length])
for i in range(len(ini)):
    A[ini[i],final[i]] += 1

## 保存
np.save('Adjacency_Matrix.npy',A) #npy版本
np.savetxt('Adjacency_Matrix.dat',A) #dat版本

## 简单可视化
# plt.contour(A)
# plt.show()