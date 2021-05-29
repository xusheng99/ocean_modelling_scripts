# 8.0更新说明
# 使用scipy的插值函数提速

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
import multiprocessing
from scipy.interpolate import interpn 


############################ 控制参数 ##############################
## 工作和存储目录
os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\')

## 流场文件参数
# flowfield_prefix = 'D:\\application\\SCS\\temp\\'  # 南海
flowfield_prefix = 'C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\'
his_timestep = 3600 # history文件中, 两个时刻之间的时间差 单位秒
his_writestep = 24 # 非初始的history文件中的时刻数目, 单位个; 注意, 初始的his_0001.nc比后续的正常的文件多一个时刻. 

## 粒子行为控制
diffusion = 'off' # 是否开启扩散? 'on'/'off'
random_speed_sigma = 0.06 # 随机速度的标准差（正态分布模式） / 波动范围的四分之一（几何概型）
stick_stop = 'off'  # 碰到陆地或者是开边界后是否停止计算改粒子

## 积分步长和时间控制
dt = 1200  # 粒子运动步长, 单位秒
total_time = 86400*30 # 总的模拟时间, 86400秒为一天
release_time = 0 + 86400*180 # 释放的初始时刻, 单位秒
timestep_numbers = int(total_time/dt) # 总共运动多少个时间步长数目, 后面要用这个变量

## 初始释放粒子范围和数量控制 
# minX = 105.5
# maxX = 124.5
# minY = 5.5
# maxY = 24.5
# N_edge = 200 # 在上述矩形中均匀释放N_edge*N_edge个粒子

minX = 115.0
maxX = 122
minY = 15
maxY = 22
N_edge = 12 # 在上述矩形中均匀释放N_edge*N_edge个粒子
####################################################################


















####################################################################
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
# 使用cdo取出表面的子数据集为surface_uv.nc, 一次性全部读取到内存中, 低配电脑容易崩.
content = Dataset('./surface_uv.nc4') 
content2 = Dataset('./Grd.nc')
mask_rho = content2.variables['mask_rho'][:].data
lon_rho = content2.variables['lon_rho'][:].data
lat_rho = content2.variables['lat_rho'][:].data
time_rho = content.variables['ocean_time'][:].data
#2016年6月18日~7月17日, 前后闭区间, 每日24个时刻, 第一个时刻为凌晨1时, 最后一个时刻为当日24时
u_rho = content.variables['u_eastward'][:].data
v_rho = content.variables['v_northward'][:].data
u_rho[u_rho>10] = .0
v_rho[v_rho>10] = .0
u_rho,v_rho = np.squeeze(u_rho),np.squeeze(v_rho)
content.close()
# 读取网格大小信息
lonmin, lonmax, latmin, latmax = np.min(lon_rho), np.max(lon_rho), np.min(lat_rho), np.max(lat_rho)
lonarray,latarray = lon_rho[0,:],lat_rho[:,0]
len_xi,len_eta = len(lonarray),len(latarray) # grid size
length = ini_particle_matrix.shape 
N = length[0] # 读取质点总数
cpu_count = multiprocessing.cpu_count() # 读取当前计算机核数
#####################################################################

def in_ocean(x1,y1):
    return [(mask_rho[y1,x1] == 1 ), (mask_rho[y1,x1+1] == 1), (mask_rho[y1+1,x1+1] == 1), (mask_rho[y1+1,x1] == 1)]

def V(lon,lat,time_):
    flag = 0  # 返回状态flag, 1表示陆地, 2表示开边界
    
    # 开边界判定
    if not ((lonmin < lon < lonmax) and (latmin< lat < latmax)):
        flag = 2
        return [0.0,0.0,flag]

    # 海陆判定        
    lonstart = lonarray[0]
    latstart = latarray[0]
    xinterval = (len_xi - 1)/(lonarray[-1] - lonstart)
    yinterval = (len_eta - 1)/(latarray[-1] - latstart)
    x = (lon - lonstart)*xinterval
    y = (lat - latstart)*yinterval
    inocean = in_ocean(int(x),int(y))
    if  not inocean[0]+inocean[1]+inocean[2]+inocean[3]:
    # if not inocean[0]*inocean[1]*inocean[2]*inocean[3]:
        flag = 1
        return [0.0,0.0,flag]

    # 插值
    points = (time_rho,latarray,lonarray)
    xi = np.array([time_,lat,lon])
    ui = float(interpn(points,u_rho,xi,fill_value = False,bounds_error = False,method = 'linear'))
    vi = float(interpn(points,v_rho,xi,fill_value = False,bounds_error = False,method = 'linear'))
    
    if diffusion == 'on':  #扩散项
        u = u + np.random.normal(0,random_speed_sigma) #正态随机扩散，0为速度均值， 0.05为标准差
        v = v + np.random.normal(0,random_speed_sigma)
        # u = u + random.uniform(-random_speed_sigma*2,random_speed_sigma*2) 
        # v = v + random.uniform(-random_speed_sigma*2,random_speed_sigma*2)

    u = ui*(360/(6371000*2*3.14159*math.cos(lat*(math.pi/180))))
    v = vi*(360/(6371000*2*3.14159))
    # t4 = time.perf_counter()
    # print("空间插值耗时:",t4-t3)
    return [u,v,flag]

# V(115.6,15.4,86400*180)
# jump(jump([115.6,15.4,86400*180]))

# 某位置，某时刻，返回一个时间步长后的状态 
# 四阶龙格库塔算法  
def jump(position):
    # position = [114,12,4000]
    lon_0 = position[0]
    lat_0 = position[1]
    time_0 = position[2]
    [u,v,flag] = V(lon_0,lat_0,time_0) 
    if flag == 0:
        [u1,v1] = [u,v]
        [u2,v2] = V( lon_0+dt*0.5*u1, lat_0+dt*0.5*v1, time_0+dt*0.5 )[0:2]
        [u3,v3] = V( lon_0+dt*0.5*u2, lat_0+dt*0.5*v2, time_0+dt*0.5 )[0:2]
        [u4,v4] = V( lon_0+dt*u3, lat_0+dt*v3, time_0+dt)[0:2]

        u,v = (1/6)*(u1+2*u2+2*u3+u4), (1/6)*(v1+2*v2+2*v3+v4) ###
        
        lon,lat = lon_0 + dt*u, lat_0 + dt*v
      
        newflag = V(lon,lat,time_0+dt)[2] # 检查新位置是否上陆或出界
        if newflag == 1: # 某一步跳上陆地返回上一步, 回到海里
            lon,lat = lon_0,lat_0
        if newflag == 2: # 某一步跳出边界,则抛弃粒子,设置成NaN
            lon,lat = 19981231,19981231
    if flag == 1: # 初始时刻就释放在陆地上的粒子, 从头到尾都别动
        lon,lat = lon_0,lat_0
    if flag == 2: # 初始时刻就在外面的粒子, 抛弃, 话说这属于初始情况设置错误
        lon,lat = 19981231,19981231
    return [lon,lat,time_0+dt]

# 定义每个质点的编号
def storage_filename(x):
    return 'particle_'+str(x+1)+'.dat'

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

# release(0,20)

# 保存初始时刻和末尾时刻的状态
def save_ini_final():
    final = np.zeros([N,2]) #维度顺序: 时间,粒子编号,经或纬
    ini = np.zeros([N,2])
    for i in range(N):
        each_particle_data = np.loadtxt('particle_'+str(i+1)+'.dat')
        final[i,:] = each_particle_data[-1,:]
        ini[i,:] = each_particle_data[0,:]
    print('初始时刻: initial.dat \n最终时刻: final.dat')
    np.savetxt('final.dat',final)
    np.savetxt('initial.dat',ini)

# 保存指定时刻的状态, 单位秒
def snapshot(time):
    time = time//dt
    shot = np.zeros([N,2]) #维度顺序: 时间,粒子编号,经或纬
    for i in range(N):
        each_particle_data = np.loadtxt('particle_'+str(i+1)+'.dat')
        shot[i,:] = each_particle_data[time,:]
    print('第',str(time+1)+"时刻的状态已经保存到: ",str(time+1)+'.dat')
    np.savetxt(str(time+1)+'.dat',shot)

if __name__ == '__main__':
    T1 = time.perf_counter()
    
    if cpu_count >=16: # 工作站
        interval = N//16
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
        for process in {p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16}:
            process.join()
        print('Child process end')
    else: # 笔记本
        interval = N//4
        p1 = Process(target=release,args=(0,interval))
        p2 = Process(target=release,args=(interval,interval))
        p3 = Process(target=release,args=(interval*2,interval))
        p4 = Process(target=release,args=(interval*3,N-3*interval))
        print('Child process start')
        for process in {p1,p2,p3,p4}:
            process.start()
        for process in {p1,p2,p3,p4}:
            process.join()
        print('Child process end')

    save_ini_final()

    T2 = time.perf_counter()
    print()
    print('共计耗时: ',(T2-T1)/1440,'小时')