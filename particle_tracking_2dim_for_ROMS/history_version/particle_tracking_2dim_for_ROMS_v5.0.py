import os
import shutil
import sys
import numpy as np
from netCDF4 import Dataset
import math
import random
import scipy
import time
import re
from multiprocessing import Process
os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\')

# 5.0 更新说明
# 使用四阶龙格库塔算法替代原来的欧拉法, 增加精度
# 运行方式: 命令行运行 python ./particle_tracking_2dim_for_ROMS_v4.0.py 

################# 初始释放粒子范围和数量控制 ##############
# minX = 112
# maxX = 118
# minY = 15
# maxY = 21

minX = 114
maxX = 114
minY = 12
maxY = 12

x = np.linspace(minX,maxX,2)
y = np.linspace(minY,maxY,2)

# for my own python script, without comma
with open('particles_initial.txt','w',encoding='utf-8') as f:
    for eachlon in x:
        for eachlat in y:
            f.write(str(eachlon)+' ')
            f.write(str(eachlat)+' ')
            f.write('3600'+' ')  #date of birth, unit:seconds
            f.write('\n')
f.close()
ini_particle_matrix = np.loadtxt('particles_initial.txt')

# 控制参数
diffusion = 'off' # 是否开启扩散? 'on'/'off'
dt = 360  #粒子运动步长, 单位秒
random_speed_sigma = 0.06 #随机速度的标准差（正态分布模式） / 波动范围的四分之一（几何概型）
total_time = 86400*90 #总的模拟时间, 86400秒为一天
timestep_numbers = int(total_time/dt) # 时间步长数目, 后面要用这个变量

flowfield_file = 'D:\\application\\SCS\\fortrack.nc4'  #南海
#flowfield_file = 'D:\\bhd_sea_model\\bhd3_refine_undone\\out_his.nc'  #东海








#####################################################################
content = Dataset(flowfield_file)
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
lon_rho = content.variables['lon_rho'][:].data
lat_rho = content.variables['lat_rho'][:].data
#timeseries = content.variables['ocean_time'][:].data  #0/3600/7200... 秒
# 读取网格大小信息
lonmin = np.min(lon_rho)
lonmax = np.max(lon_rho)
latmin = np.min(lat_rho)
latmax = np.max(lat_rho)
lonarray = lon_rho[0,:]
latarray = lat_rho[:,0]
len_xi = len(lonarray)
len_eta = len(latarray)  #grid size
length = ini_particle_matrix.shape 
length = length[0]
N = length  # 读取质点总数
#####################################################################

# 返回任意位置，任意时间的速度;  单位：度，度，秒； 时间从0开始;
def V(lon,lat,time): #
    # 将经纬度转化为0~229,0~239之间的一个小数坐标
    xinterval = lonarray[1] - lonarray[0]
    yinterval = latarray[1] - latarray[0]
    x = 114514
    y = 114514
    for i in range(1,len_xi): #排除边缘网格 从1开始而不是0
        if lonarray[i] == lon:
            x = i
            break
        else:
            if i == len_xi - 2:
                pass
                break #return后的break是无意义的
            if lonarray[i] < lon < lonarray[i+1]:
                x = i + ((lon - lonarray[i])/(xinterval))
                break
    for j in range(1,len_eta): #排除边缘网格 
        if latarray[j] == lat:
            y = j
            break
        else:
            if j == len_eta - 2:  #到达边缘就停止计算，边缘插值误差大
                break                #break
            if latarray[j] < lat < latarray[j+1]:
                y = j + ((lat - latarray[j])/(yinterval))  #这里的j写成了i...查了好久才查出来
                break

    if x == 114514 or y == 114514: #未能正确赋值, 说明碰到了开边界, 黏附
        print(lon,lat,'now arrived at edge')
        return ['obc','obc'] # return完了就跳出这个函数

    time = time/3600  #刚好对应0秒-时间索引0；3600秒-时间索引1...
    #print(x,y,time)

    time1 = int(time//1)
    time2 = time1 + 1
    x1 = int(x//1)  #整数部分  
    x2 = x1 + 1
    y1 = int(y//1)
    y2 = y1 + 1

    #print(x,y,time)
    #新算法，基于空间中最近8点的反距离权重插值(IDW)
    u = []
    u.append(content.variables['u_eastward'][time1,-1,y1,x1])
    u.append(content.variables['u_eastward'][time1,-1,y1,x2])
    u.append(content.variables['u_eastward'][time1,-1,y2,x2])
    u.append(content.variables['u_eastward'][time1,-1,y2,x1])
    u.append(content.variables['u_eastward'][time2,-1,y1,x1])
    u.append(content.variables['u_eastward'][time2,-1,y1,x2])
    u.append(content.variables['u_eastward'][time2,-1,y2,x2])
    u.append(content.variables['u_eastward'][time2,-1,y2,x1])

    v = []
    v.append(content.variables['v_northward'][time1,-1,y1,x1])
    v.append(content.variables['v_northward'][time1,-1,y1,x2])
    v.append(content.variables['v_northward'][time1,-1,y2,x2])
    v.append(content.variables['v_northward'][time1,-1,y2,x1])
    v.append(content.variables['v_northward'][time2,-1,y1,x1])
    v.append(content.variables['v_northward'][time2,-1,y1,x2])
    v.append(content.variables['v_northward'][time2,-1,y2,x2])
    v.append(content.variables['v_northward'][time2,-1,y2,x1])

    for point in range(4):  #只需检查第一个时刻的周边四点, 提高效率
        if u[point].mask:  #如果周边有任意一个点是陆点, 则黏附, 返回零速率
            return ['land','land']

    def distance(m1,n1,k1,m2,n2,k2):  # 欧几里得距离
        return math.sqrt((m1-m2)*(m1-m2)+(n1-n2)*(n1-n2)+(k1-k2)*(k1-k2))

    d = []
    d.append(distance(x,y,time,x1,y1,time1))
    d.append(distance(x,y,time,x2,y1,time1))
    d.append(distance(x,y,time,x2,y2,time1))
    d.append(distance(x,y,time,x1,y2,time1))
    d.append(distance(x,y,time,x1,y1,time2))
    d.append(distance(x,y,time,x2,y1,time2))
    d.append(distance(x,y,time,x2,y2,time2))
    d.append(distance(x,y,time,x1,y2,time2))

    for point in range(len(d)):  # 正卡在了整数网格上，手动处理分母为零的情况
        if abs(d[point]) < 1e-15:
            return [1*u[point],1*v[point]]
        
    def func_d(distance_list):  # 每点的权重因子, 目前使用反比关系, 或可选平方反比关系
        fd = []
        for point in range(len(distance_list)):
            fd.append(1/distance_list[point])
        return fd
    fd_list = func_d(d)

    def weight(fd_list):  # 每点的权
        weight = []
        for point in range(len(fd_list)):
            weight.append((fd_list[point])/(sum(fd_list)))
        return weight
    weight = weight(fd_list)

    for point in range(8): 
        u[point] = 1*u[point].data #乘1改变类型
        v[point] = 1*v[point].data
    
    # 加权平均, 即 权向量u/v 和 速率向量weight 点乘
    u_cross_weight = [i*j for i,j in zip(u,weight)]
    v_cross_weight = [i*j for i,j in zip(v,weight)]
    u = sum(u_cross_weight)
    v = sum(v_cross_weight)

    # 改单位, 改成degree/second
    u = u*(360/(6371000*2*3.14159*math.cos(lat*(math.pi/180))))
    v = v*(360/(6371000*2*3.14159))
    return [u,v]

# 某位置，某时刻，返回一个时间步长后的状态 
# 四阶龙格库塔算法  巨慢...
def jump(position):
    lon_0 = position[0]
    lat_0 = position[1]
    time_0 = position[2]
    [u,v] = V(lon_0,lat_0,time_0)
    if u != 0 and v != 0:
        V1 = V( lon_0, lat_0, time_0 )
        u1,v1 = V1[0],V1[1]
        V2 = V( lon_0+dt*0.5*u1, lat_0+dt*0.5*v1, time_0+dt*0.5 )
        u2,v2 = V2[0],V2[1]
        V3 = V( lon_0+dt*0.5*u2, lat_0+dt*0.5*v2, time_0+dt*0.5 )
        u3,v3 = V3[0],V3[1]
        V4 = V( lon_0+dt*u3, lat_0+dt*v3, time_0+dt)
        u4,v4 = V4[0],V4[1]

        u = (1/6)*(u1+2*u2+2*u3+u4)
        v = (1/6)*(v1+2*v2+2*v3+v4)
        
        if diffusion == 'on':  #扩散项
            u = u + np.random.normal(0,random_speed_sigma) #正态随机扩散，0为速度均值， 0.05为标准差
            v = v + np.random.normal(0,random_speed_sigma)
        #u = u + random.uniform(-random_speed_sigma*2,random_speed_sigma*2) #随机扩散，0为速度均值， 0.05为标准差
        #v = v + random.uniform(-random_speed_sigma*2,random_speed_sigma*2)

        lon = lon_0 + dt*u
        lat = lat_0 + dt*v
    else:
        lon = lon_0
        lat = lat_0
    return [lon,lat,time_0+dt]

# 定义每个质点的编号
def storage_filename(x):
    return 'particle_'+str(x+1)+'.dat'

# 释放指定编号范围的粒子
def release(basenum,NN):
    for particle_num in range(basenum,basenum+NN):  #对每一个粒子
        t1 = time.perf_counter()
        with open(storage_filename(particle_num),'w',encoding='utf-8') as f:  #开启对应编号的文件
            print('now tracking: ','particle_'+str(particle_num+1))
            condition = list(ini_particle_matrix[particle_num][:])  #读取对应编号粒子的初始状态
            f.write(str(condition[0])+' '+str(condition[1])+'\n')
            count = 0
            # 开始运动
            while condition[2] < total_time + ini_particle_matrix[0][2]: #超过截至时间，停止运行
                temp = condition
                condition = jump(condition)
                count = count + 1
                f.write(str(condition[0])+' '+str(condition[1])+'\n')
                if temp[0] == condition[0] and temp[1] == condition[1]:  #如果停止不动，说明撞岸，停止运行节约资源
                    #print('hit the land/edge, end!')
                    break
        f.close()
        t2 = time.perf_counter()
        # print('该质点运动了',count,'步/',int((count*dt)/(3600*24)),'天, 耗时',(t2-t1),'秒','\n')

# 将按粒子编号分别存储的文件, 重新整合成按时刻存储的文件
def reconstruct():  
    print('正在重新写结果')
    interval = 900  # 每隔interval个时间步长(dt)写入一次全场的所有粒子的位置到同一个文件里

    # 读取数据, 如果interval过小, 可能会爆内存
    a = np.zeros([(timestep_numbers//interval),N,2]) #维度顺序: 时间,粒子编号,经或纬
    for i in range(N):
        each_particle_data = np.loadtxt('particle_'+str(i+1)+'.dat')
        a[:,i,:] = each_particle_data[0:-1:interval,:]
    print('all data loaded in,now reconstructing...')

    for each_timestep in range(timestep_numbers//interval):
        print(each_timestep)
        with open(str((each_timestep+1)*interval)+'.dat','w',encoding='utf-8') as f:
            for each_particle in range(N):
                f.write(str(a[each_timestep,each_particle,0])+' '+str(a[each_timestep,each_particle,1])+'\n')
        f.close()
    print('重写完成')


if __name__ == '__main__':
    # 按粒子, 把任务分配到八个进程计算
    T1 = time.perf_counter()
    
    interval = N//8
    p1 = Process(target=release,args=(0,interval))
    p2 = Process(target=release,args=(interval,interval))
    p3 = Process(target=release,args=(interval*2,interval))
    p4 = Process(target=release,args=(interval*3,interval))
    p5 = Process(target=release,args=(interval*4,interval))
    p6 = Process(target=release,args=(interval*5,interval))
    p7 = Process(target=release,args=(interval*6,interval))
    p8 = Process(target=release,args=(interval*7,N-7*interval))
    print('Child process start')
    p1.start()
    p2.start()
    p3.start()
    p4.start()
    p5.start()
    p6.start()
    p7.start()
    p8.start()
    p1.join()
    p2.join()
    p3.join()
    p4.join()
    p5.join()
    p6.join()
    p7.join()
    p8.join()
    print('Child process end')

    reconstruct() 

    T2 = time.perf_counter()
    print()
    print('共计耗时: ',T2-T1,'秒')



