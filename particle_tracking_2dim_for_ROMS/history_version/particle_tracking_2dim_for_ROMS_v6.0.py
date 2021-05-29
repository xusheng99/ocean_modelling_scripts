# 6.0 更新说明
# 支持离散的多个his文件
# 增加了对陆地和开边界的处理
# 进一步降低性能~ [😀]
# 运行方式: 命令行运行 python ./particle_tracking_2dim_for_ROMS_v4.0.py 

import os
import shutil
import sys
import numpy as np
from netCDF4 import Dataset
import math
import random
import time
# import re
from multiprocessing import Process
os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\')


############################ 控制参数 ##############################
## 流场文件参数
flowfield_prefix = 'D:\\application\\SCS\\'  # 南海
his_timestep = 3600 # history文件中, 两个时刻之间的时间差 单位秒
his_writestep = 24 # 非初始的history文件中的时刻数目, 单位个; 注意, his_0001.nc比后续的正常的文件多一个时刻. 

## 粒子行为控制
diffusion = 'off' # 是否开启扩散? 'on'/'off'
random_speed_sigma = 0.06 # 随机速度的标准差（正态分布模式） / 波动范围的四分之一（几何概型）
stick_stop = 'off'  # 碰到陆地或者是开边界后是否停止计算改粒子

## 积分步长和时间控制
dt = 900  # 粒子运动步长, 单位秒
total_time = 86400*30 # 总的模拟时间, 86400秒为一天
release_time = 10800 # 释放的初始时刻, 单位秒
timestep_numbers = int(total_time/dt) # 总共运动多少个时间步长数目, 后面要用这个变量

## 初始释放粒子范围和数量控制 
minX = 106
maxX = 124
minY = 6
maxY = 24
N_edge = 51 # 在上述矩形中均匀释放N_edge*N_edge个粒子

####################################################################





















# 第一个问题: 如何判定粒子有没有跑到陆地上去?   每一步都搞一个多边形判定?
# 第二个问题: 怎么实现镜面反射?
# 第三个问题: 边界的V插值如何去做? 
# 目前使用的方法中, 任意连续位置的速度是根据对角四个网格的值插出来的
# 靠近陆地的点是的周边是缺失值的, 如何处理?  缺失的值采用剩余有效点的平均值

x = np.linspace(minX,maxX,N_edge)
y = np.linspace(minY,maxY,N_edge)
# for my own python script, without comma
with open('particles_initial.txt','w',encoding='utf-8') as f:
    for eachlon in x:
        for eachlat in y:
            f.write(str(eachlon)+' ')
            f.write(str(eachlat)+' ')
            f.write(str(release_time)+' ')  #date of birth, unit:seconds
            f.write('\n')
f.close()
ini_particle_matrix = np.loadtxt('particles_initial.txt')

def whichfile(time):
        gulugulu = his_timestep*his_writestep
        if time == 0:
            return 'his_0001.nc'
        else:
            if time%(gulugulu) == 0:
                return 'his_'+"{:0>4d}".format(int(time//gulugulu))+'.nc'
            else:
                return 'his_'+"{:0>4d}".format(int(time//gulugulu + 1))+'.nc'

#####################################################################
content = Dataset(flowfield_prefix + whichfile(1))
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
lon_rho = content.variables['lon_rho'][:].data
lat_rho = content.variables['lat_rho'][:].data
content.close()
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
def V(lon,lat,time_): 
    t1 = time.perf_counter()
    flag = 0  # 返回状态flag, 1表示
    # lon = 110.4
    # lat = 18.62
    # time_ = 4000
    # 将经纬度转化为0~229,0~239之间的一个小数坐标
    lonstart = lonarray[0]
    latstart = latarray[0]
    xinterval = (len_xi - 1)/(lonarray[-1] - lonstart)
    yinterval = (len_eta - 1)/(latarray[-1] - latstart)
    x = (lon - lonstart)*xinterval
    y = (lat - latstart)*yinterval
    xmin,ymin = 0,0
    xmax,ymax = (lonarray[-1] - lonstart)*xinterval,(latarray[-1] - latstart)*yinterval
 
    if x > xmax or x < xmin or y > ymax or y < ymin: #离开开边界判定
        flag = 2
        return [0.0,0.0,flag]

    # if x == 114514 or y == 114514: #未能正确赋值, 说明碰到了开边界, 返回开边界标识符
    #     print(lon,lat,'now arrived at edge')
    #     return ['obc','obc'] # return完了就跳出这个函数

    time_ = time_/3600  #刚好对应0秒-时间索引0；3600秒-时间索引1...
    time1 = int(time_//1)
    time2 = time1 + 1
    x1 = int(x//1)  #整数部分  
    x2 = x1 + 1
    y1 = int(y//1)
    y2 = y1 + 1

    t2 = time.perf_counter()
    print("完成坐标转换的时间:",t2-t1)

    def is_edge_or_boundary(x1,y1):
        return [(mask_rho[y1,x1] == 1 ), (mask_rho[y1,x1+1] == 1), (mask_rho[y1+1,x1+1] == 1), (mask_rho[y1+1,x1] == 1)]
    def append_0(u,v):
        u.append(0.0)
        v.append(0.0)
    
    u = []
    v = []
    inocean = is_edge_or_boundary(x1,y1) # 周边四个点是否在海里的布尔值数组
    # 陆地点判定
    if  not inocean[0]+inocean[1]+inocean[2]+inocean[3]:
    # if not inocean[0]*inocean[1]*inocean[2]*inocean[3]:
        flag = 1
        return [0.0,0.0,flag]

    content = Dataset(flowfield_prefix + whichfile(time1*his_timestep))
    if inocean[0]:
        u.append(content.variables['u_eastward'][time1%his_writestep,-1,y1,x1].data)
        v.append(content.variables['v_northward'][time1%his_writestep,-1,y1,x1].data)
    else:
        append_0(u,v)
    if inocean[1]:
        u.append(content.variables['u_eastward'][time1%his_writestep,-1,y1,x2].data)
        v.append(content.variables['v_northward'][time1%his_writestep,-1,y1,x2].data)
    else:
        append_0(u,v)
    if inocean[2]:
        u.append(content.variables['u_eastward'][time1%his_writestep,-1,y2,x2].data)
        v.append(content.variables['v_northward'][time1%his_writestep,-1,y2,x2].data)
    else:
        append_0(u,v)
    if inocean[3]:
        u.append(content.variables['u_eastward'][time1%his_writestep,-1,y2,x1].data)
        v.append(content.variables['v_northward'][time1%his_writestep,-1,y2,x1].data)
    else:
        append_0(u,v)

    # 也许time1和time2分布在两个不同的文件里, 需要对content做出更新
    if whichfile(time1) != whichfile(time2):
        content = Dataset(flowfield_prefix + whichfile(time2*his_timestep))

    if inocean[0]:
        u.append(content.variables['u_eastward'][time2%his_writestep,-1,y1,x1].data)
        v.append(content.variables['v_northward'][time2%his_writestep,-1,y1,x1].data)
    else:
        append_0(u,v)
    if inocean[1]:
        u.append(content.variables['u_eastward'][time2%his_writestep,-1,y1,x2].data)
        v.append(content.variables['v_northward'][time2%his_writestep,-1,y1,x2].data)
    else:
        append_0(u,v)
    if inocean[2]:
        u.append(content.variables['u_eastward'][time2%his_writestep,-1,y2,x2].data)
        v.append(content.variables['v_northward'][time2%his_writestep,-1,y2,x2].data)
    else:
        append_0(u,v)
    if inocean[3]:
        u.append(content.variables['u_eastward'][time2%his_writestep,-1,y2,x1].data)
        v.append(content.variables['v_northward'][time2%his_writestep,-1,y2,x1].data)
    else:
        append_0(u,v)
    content.close()
    t3 = time.perf_counter()
    print("完成周边八点值读取耗时:",t3-t2)

    # # 这个地方有问题, 无法实现反射, 而且大大扩展了陆地的边界
    # for point in range(4):  #只需检查第一个时刻的周边四点, 提高效率
    #     if u[point].mask:  #如果周边有任意一个点是陆点, 则黏附, 返回陆地标识符
    #         return ['land','land']

    def distance(m1,n1,k1,m2,n2,k2):  # 欧几里得距离
        return math.sqrt((m1-m2)*(m1-m2)+(n1-n2)*(n1-n2)+(k1-k2)*(k1-k2))

    d = []
    d.append(distance(x,y,time_,x1,y1,time1))
    d.append(distance(x,y,time_,x2,y1,time1))
    d.append(distance(x,y,time_,x2,y2,time1))
    d.append(distance(x,y,time_,x1,y2,time1))
    d.append(distance(x,y,time_,x1,y1,time2))
    d.append(distance(x,y,time_,x2,y1,time2))
    d.append(distance(x,y,time_,x2,y2,time2))
    d.append(distance(x,y,time_,x1,y2,time2))

    for point in range(len(d)):  # 正巧卡在了整数编号的网格点上，手动处理分母为零的情况
        if abs(d[point]) < 1e-15:
            return [1*u[point]*(360/(6371000*2*3.14159*math.cos(lat*(math.pi/180)))),1*v[point]*(360/(6371000*2*3.14159)),flag]
        
    def func_d(distance_list):  # 每点的权重因子, 平方反比关系
        fd = []
        for point in range(len(distance_list)):
            fd.append(1/(distance_list[point]*distance_list[point]))
        return fd
    fd_list = func_d(d)

    def weight(fd_list):  # 每点的权
        weight = []
        for point in range(len(fd_list)):
            weight.append((fd_list[point])/(sum(fd_list)))
        return weight
    weight = weight(fd_list)

    # for point in range(8): 
    #     u[point] = 1*u[point].data #乘1改变类型
    #     v[point] = 1*v[point].data
    
    # 加权平均, 即 权向量u/v 和 速率向量weight 点乘
    u_cross_weight = [i*j for i,j in zip(u,weight)]
    v_cross_weight = [i*j for i,j in zip(v,weight)]
    u = sum(u_cross_weight)
    v = sum(v_cross_weight)

    if diffusion == 'on':  #扩散项
        u = u + np.random.normal(0,random_speed_sigma) #正态随机扩散，0为速度均值， 0.05为标准差
        v = v + np.random.normal(0,random_speed_sigma)
        #u = u + random.uniform(-random_speed_sigma*2,random_speed_sigma*2) 
        #v = v + random.uniform(-random_speed_sigma*2,random_speed_sigma*2)
    # print(u,v)

    # 改单位, 改成degrees/second
    u = u*(360/(6371000*2*3.14159*math.cos(lat*(math.pi/180))))
    v = v*(360/(6371000*2*3.14159))
    t4 = time.perf_counter()
    print("空间插值耗时:",t4-t3)
    return [u,v,flag]

# V(113.90131752706657,12.109197143383554,123123)

# 某位置，某时刻，返回一个时间步长后的状态 
# 四阶龙格库塔算法  
def jump(position):
    # position = [114,12,4000]
    lon_0 = position[0]
    lat_0 = position[1]
    time_0 = position[2]
    [u,v,flag] = V(lon_0,lat_0,time_0) 
    if flag == 0:
        [u1,v1,f1] = V( lon_0, lat_0, time_0 )
        [u2,v2,f2] = V( lon_0+dt*0.5*u1, lat_0+dt*0.5*v1, time_0+dt*0.5 )
        [u3,v3,f3] = V( lon_0+dt*0.5*u2, lat_0+dt*0.5*v2, time_0+dt*0.5 )
        [u4,v4,f4] = V( lon_0+dt*u3, lat_0+dt*v3, time_0+dt)

        u,v = (1/6)*(u1+2*u2+2*u3+u4), (1/6)*(v1+2*v2+2*v3+v4)
        
        lon,lat = lon_0 + dt*u, lat_0 + dt*v
      
        newflag = V(lon,lat,time_0+dt)[2] # 检查新位置是否上陆或出界
        if newflag == 1: # 某一步跳上陆地返回上一步, 回到海里
            lon,lat = lon_0,lat_0
        if newflag == 2: # 某一步跳出边界,则抛弃粒子,设置成NaN
            lon,lat = 19981231,19981231
    if flag == 1: # 初始时刻就释放在陆地上的粒子, 从头到尾都别动
        lon,lat = lon_0,lat_0
    if flag == 2: # 初始时刻就在外面的粒子, 抛弃
        lon,lat = 19981231,19981231
    return [lon,lat,time_0+dt]

# jump([114,12,4000])

# 定义每个质点的编号
def storage_filename(x):
    return 'particle_'+str(x+1)+'.dat'

def release_single(initialcondition):
    # initialcondition =[114,12,4000]
    with open('1.dat','w',encoding='utf-8') as f:  
        condition = initialcondition 
        f.write(str(condition[0])+' '+str(condition[1])+'\n') #写入初始状态到第一行
        count = 0
        while condition[2] < total_time + initialcondition[2]:
            temp = condition
            condition = jump(condition)
            print(condition)
            count = count + 1
            f.write(str(condition[0])+' '+str(condition[1])+'\n')
            if stick_stop == 'on':
                if temp[0] == condition[0] and temp[1] == condition[1]: 
                    #print('hit the land/edge, end!')
                    break
    f.close()

# # release_single([114,12,4000])
# release_single([106.3,14.25,86080]) #陆地
# release_single([111.05,19.58,86080])

# test_position = [111.012, 19.59, 3110080]
# release_single(test_position)
# jump(test_position)
# jump(jump(test_position))
# V(test_position[0],test_position[1],test_position[2])

# jump([112.0,12.0,3600])
# V(112.0,12.0,3600)

# 释放指定编号范围的粒子
def release(basenum,NN):
    for particle_num in range(basenum,basenum+NN):  #对每一个粒子
        # t1 = time.perf_counter()
        with open(storage_filename(particle_num),'w',encoding='utf-8') as f:  #开启对应编号的文件
            print('now tracking: ','particle_'+str(particle_num+1))
            condition = list(ini_particle_matrix[particle_num][:])  #读取对应编号粒子的初始状态
            f.write(str(condition[0])+' '+str(condition[1])+'\n') #写入初始状态到第一行
            count = 0
            # 开始运动
            while condition[2] < total_time + ini_particle_matrix[0][2]: #超过截至时间，停止运行
                temp = condition
                condition = jump(condition) #iteration
                count = count + 1
                f.write(str(condition[0])+' '+str(condition[1])+'\n')
                if stick_stop == 'on':
                    if temp[0] == condition[0] and temp[1] == condition[1]:  #如果停止不动，停止运行节约资源
                        #print('hit the land/edge, end!')
                        break
        f.close()
        # t2 = time.perf_counter()
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

    # release_single([114,12,4000])


