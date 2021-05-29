# 9.0 更新说明
# 使用了新编写的inocean函数


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
from multiprocessing import cpu_count
import functools
import platform
import linecache
from myfunc.nc_operate import inocean

################### 命令行参数  ##################
# workdir = sys.argv[1] 
# release_time = (int(sys.argv[2])-1)*86400 + 0 # 释放的初始时刻, 单位秒
workdir = 'C:\\Users\\xusheng\\desktop\\particle_tracking_results\\test2\\'
release_time = 0

############################ 控制参数 ##############################
## 工作目录
os.chdir(workdir)

## 流场文件参数
system = platform.system()
if system == 'Windows':
    flowfield_prefix = 'D:\\application\\SCS\\'  
if system == 'Linux':
    flowfield_prefix = '/mnt/d/application/SCS/'
his_timestep = 3600 # history文件中, 两个时刻之间的时间差 单位秒
his_writestep = 24 # 非初始的history文件中的时刻数目, 单位个; 注意, 初始的his_0001.nc比后续的正常的文件多一个时刻. 

## 粒子行为控制
diffusion = 'off' # 是否开启扩散? 'on'/'off'
random_speed_sigma = 0.06 # 随机速度的标准差（正态分布模式） / 波动范围的四分之一（几何概型）
stick_stop = 'off'  # 碰到陆地或者是开边界后是否停止计算改粒子

## 积分步长和时间控制
dt = 1200  # 粒子运动步长, 单位秒
total_time = 86400*3 # 总的模拟时间, 86400秒为一天
timestep_numbers = int(total_time/dt) # 总共运动多少个时间步长数目, 后面要用这个变量

## 初始释放粒子范围和数量控制 
minX = 113.3
maxX = 114.3
minY = 21.5
maxY = 22.5
N_edge = 100 # 在上述矩形中均匀释放N_edge*N_edge个粒子
####################################################################



















## 创建初始粒子文件
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
#####################################################################
content = Dataset(flowfield_prefix + 'surface_uv.nc4')
content2 = Dataset(flowfield_prefix+'DATA/'+'Grd.nc')
mask_rho = content2.variables['mask_rho'][:].data
lon_rho = content2.variables['lon_rho'][:].data
lat_rho = content2.variables['lat_rho'][:].data
time_rho = content.variables['ocean_time'][:].data
print('now loading u & v in rho points, this might take several minutes\n')
u_rho = np.squeeze(content.variables['u_eastward'][:].data)
v_rho = np.squeeze(content.variables['v_northward'][:].data)
print('all data loaded into RAM\n')
u_rho[u_rho>10] = .0
v_rho[v_rho>10] = .0
content2.close()
content.close()
# 读取网格大小信息
lonmin, lonmax, latmin, latmax = np.min(lon_rho), np.max(lon_rho), np.min(lat_rho), np.max(lat_rho)
lonarray,latarray = lon_rho[0,:],lat_rho[:,0]
len_xi,len_eta = len(lonarray),len(latarray) # grid size
length = ini_particle_matrix.shape 
N = length[0] # 读取质点总数
cpu_count = cpu_count() # 读取当前计算机核数
time_rho0 = time_rho[0] # 读取到u_rho和v_rho索引的初始时刻, 大部分时候给定的流场文件都是从中间截取的.
#####################################################################

# 返回任意位置，任意时间的速度;  单位：度，度，秒； 时间从0开始;
# @functools.lru_cache()
def V(lon,lat,time_): 
    # t1 = time.perf_counter()
    lon = 113.5
    lat = 22.5
    time_ = 1200
    flag = 0  # 返回状态flag, 1表示陆地, 2表示开边界

    # 开边界判定
    if not ((lonmin < lon < lonmax) and (latmin< lat < latmax)):
        flag = 2
        return [0.0,0.0,flag]
    
    # 将经纬度坐标 转化为 格点坐标
    lonstart = lonarray[0]
    latstart = latarray[0]
    xinterval = (len_xi - 1)/(lonarray[-1] - lonstart)
    yinterval = (len_eta - 1)/(latarray[-1] - latstart)
    x = (lon - lonstart)*xinterval
    y = (lat - latstart)*yinterval

    time_ = time_ - time_rho0 #!!!!!!!!!! 
    time_ = time_/3600  #刚好对应0秒-时间索引0；3600秒-时间索引1...
    time1 = int(time_//1)
    time2 = time1 + 1
    x1 = int(x//1)  #整数部分  
    x2 = x1 + 1
    y1 = int(y//1)
    y2 = y1 + 1
    
    def is_edge_or_boundary(x1,y1):
        return [(mask_rho[y1,x1] == 1 ), (mask_rho[y1,x1+1] == 1), (mask_rho[y1+1,x1+1] == 1), (mask_rho[y1+1,x1] == 1)]
    
    # 海陆判定
    _inocean = is_edge_or_boundary(x1,y1) # 周边四个点是否在海里的布尔值数组
    if (not inocean([lon,lat])): #!!!
        flag = 1
        return [0.0,0.0,flag]

    if inocean([lon,lat]) and (not _inocean[0]+_inocean[1]+_inocean[2]+_inocean[3]):
        print('出现了inocean判断海点而无uv值的情况!')
        raise Exception
    # t2 = time.perf_counter()
    # print("完成坐标转换的时间:",t2-t1)
    
    u = np.zeros([8])
    v = np.zeros([8])
    if _inocean[0]:
        u[0] = u_rho[time1,y1,x1]
        v[0] = v_rho[time1,y1,x1]
    if _inocean[1]:
        u[1] = u_rho[time1,y1,x2]
        v[1] = v_rho[time1,y1,x2]
    if _inocean[2]:
        u[2] = u_rho[time1,y2,x2]
        v[2] = v_rho[time1,y2,x2]
    if _inocean[3]:
        u[3] = u_rho[time1,y2,x1]
        v[3] = v_rho[time1,y2,x1]

    if _inocean[0]:
        u[4] = u_rho[time2,y1,x1]
        v[4] = v_rho[time2,y1,x1]
    if _inocean[1]:
        u[5] = u_rho[time2,y1,x2]
        v[5] = v_rho[time2,y1,x2]
    if _inocean[2]:
        u[6] = u_rho[time2,y2,x2]
        v[6] = v_rho[time2,y2,x2]
    if _inocean[3]:
        u[7] = u_rho[time2,y2,x1]
        v[7] = v_rho[time2,y2,x1]
    # t3 = time.perf_counter()
    # print("完成周边八点值读取耗时:",t3-t2)

    def distance(m1,n1,k1,m2,n2,k2):  # 欧几里得距离
        return math.sqrt((m1-m2)*(m1-m2)+(n1-n2)*(n1-n2)+(k1-k2)*(k1-k2))

    d = np.zeros([8])
    d[0] = distance(x,y,time_,x1,y1,time1)
    d[1] = distance(x,y,time_,x2,y1,time1)
    d[2] = distance(x,y,time_,x2,y2,time1)
    d[3] = distance(x,y,time_,x1,y2,time1)
    d[4] = distance(x,y,time_,x1,y1,time2)
    d[5] = distance(x,y,time_,x2,y1,time2)
    d[6] = distance(x,y,time_,x2,y2,time2)
    d[7] = distance(x,y,time_,x1,y2,time2)

    for point in range(len(d)):  # 正巧卡在了整数编号的网格点上，手动处理分母为零的情况
        if abs(d[point]) < 1e-15:
            return [1*u[point]*(360/(6371000*2*3.14159*math.cos(lat*(math.pi/180)))),1*v[point]*(360/(6371000*2*3.14159)),flag]
        
    def func_d(distance_list):  # 每点的权重因子, 平方反比关系
        fd = np.zeros([8])
        for point in range(len(distance_list)):
            fd[point] = (1/(distance_list[point]*distance_list[point]))
        return fd
    fd_list = func_d(d)

    def weight(fd_list):  # 每点的权
        weight = np.zeros([8])
        for point in range(len(fd_list)):
            weight[point] = ((fd_list[point])/(sum(fd_list)))
        return weight
    weight = weight(fd_list)
    
    # 加权平均, 即 权向量u/v 和 速率向量weight 点乘
    u = sum([i*j for i,j in zip(u,weight)])
    v = sum([i*j for i,j in zip(v,weight)])

    if diffusion == 'on':  #扩散项
        u = u + np.random.normal(0,random_speed_sigma) #正态随机扩散，0为速度均值， 0.05为标准差
        v = v + np.random.normal(0,random_speed_sigma)
        #u = u + random.uniform(-random_speed_sigma*2,random_speed_sigma*2) 
        #v = v + random.uniform(-random_speed_sigma*2,random_speed_sigma*2)

    # 改单位, 改成degrees/second
    # u = u*(360/(6371000*2*3.14159*math.cos(lat*(math.pi/180))))
    # v = v*(360/(6371000*2*3.14159))
    # t4 = time.perf_counter()
    # print("空间插值耗时:",t4-t3)
    return [u,v,flag]

# 某位置，某时刻，返回一个时间步长后的状态 
# 四阶龙格库塔算法  
def jump(position):
    lon_0 = position[0]
    lat_0 = position[1]
    time_0 = position[2]
    [u,v,flag] = V(lon_0,lat_0,time_0) 
    if flag == 0:
        [u1,v1] = [u,v]
        [u2,v2] = V( lon_0+dt*0.5*u1, lat_0+dt*0.5*v1, time_0+dt*0.5 )[0:2]
        [u3,v3] = V( lon_0+dt*0.5*u2, lat_0+dt*0.5*v2, time_0+dt*0.5 )[0:2]
        [u4,v4] = V( lon_0+dt*u3, lat_0+dt*v3, time_0+dt)[0:2]

        u,v = (1/6)*(u1+2*u2+2*u3+u4), (1/6)*(v1+2*v2+2*v3+v4)
        
        lon,lat = lon_0 + dt*u, lat_0 + dt*v
      
        newflag = V(lon,lat,time_0+dt)[2] # 检查新位置是否上陆或出界
        if newflag == 1: # 某一步跳上陆地返回上一步, 回到海里
            lon,lat = lon_0,lat_0
        if newflag == 2: # 某一步跳出边界,则抛弃粒子,设置成NaN
            lon,lat = 19981231,19981231
    if flag == 1: # 初始时刻就释放在陆地上的粒子, 从头到尾都别动
        lon,lat = lon_0,lat_0
    if flag == 2: # 
        lon,lat = 19981231,19981231
    return [lon,lat,time_0+dt]

# 定义每个质点的编号
def storage_filename(x):
    return 'particle_'+str(x+1)+'.dat'

# 释放指定编号范围的粒子
def release(basenum,NN):
    result = np.zeros([NN,3,timestep_numbers])
    for particle_num in range(basenum,basenum+NN):  #对每一个粒子
        print('now tracking: ','particle_'+str(particle_num+1))
        condition = list(ini_particle_matrix[particle_num][:])  #读取对应编号粒子的初始状态
        count = 0
        while count < timestep_numbers:
            result[particle_num,:,count] = condition
            condition = jump(condition)
            count = count + 1
    np.save(str(basenum)+'_'+str(basenum+NN)+'.npy',result)



if __name__ == '__main__':
    T1 = time.perf_counter()
    print('tracking {} particles\n'.format(N))
    if cpu_count >= 30: # 工作站
        interval = N//16
        print(interval)
        p1 = Process(target=release,args=(0,interval))
        p2 = Process(target=release,args=(interval,interval))
        p3 = Process(target=release,args=(interval*2,interval))
        p4 = Process(target=release,args=(interval*3,interval))
        p5 = Process(target=release,args=(interval*4,interval))
        p6 = Process(target=release,args=(interval*5,interval))
        p7 = Process(target=release,args=(interval*6,interval))
        p8 = Process(target=release,args=(interval*7,interval))
        p9 = Process(target=release,args=(interval*8,interval))
        p10 = Process(target=release,args=(interval*9,interval))
        p11 = Process(target=release,args=(interval*10,interval))
        p12 = Process(target=release,args=(interval*11,interval))
        p13 = Process(target=release,args=(interval*12,interval))
        p14 = Process(target=release,args=(interval*13,interval))
        p15 = Process(target=release,args=(interval*14,interval))
        p16 = Process(target=release,args=(interval*15,N-15*interval))
        print('Child process start')
        for process in {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16}:
            process.start()
            time.sleep(5)
        for process in {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16}:
            process.join()
        print('Child process end')
    else: # 笔记本
        interval = N//4
        print(interval)
        p1 = Process(target=release,args=(0,interval))
        p2 = Process(target=release,args=(interval,interval))
        p3 = Process(target=release,args=(interval*2,interval))
        p4 = Process(target=release,args=(interval*3,N-3*interval))
        print('Child process start')
        for process in {p1,p2,p3,p4}:
            process.start()
            time.sleep(5)
        for process in {p1,p2,p3,p4}:
            process.join()
        print('Child process end')

    print('正在删除冗余文件')
    for i in range(N):
        os.remove('particle_'+str(i+1)+'.dat')
    print('文件删除完成')
    
    T2 = time.perf_counter()
    print()
    print('共计耗时: ',(T2-T1)/1440,'小时')