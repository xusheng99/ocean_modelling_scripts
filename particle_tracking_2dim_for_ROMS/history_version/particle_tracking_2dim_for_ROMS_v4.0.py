import os
import shutil
import sys
import numpy as np

import math
import random
import scipy
import time
import re
from multiprocessing import Process

# 4.0 更新说明
# 区域释放 更改存储方式, 使之支持ftle的计算
# 运行方式: 命令行运行 python ./particle_tracking_2dim_for_ROMS_v4.0.py 
# 此版本如有粒子触岸, 则停止运行, 导致各粒子文件长度不一致, 会在后面的reconstruct中报错
os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\')



################# 初始释放粒子范围和数量控制 ##############
minX = 112
maxX = 118
minY = 15
maxY = 21

minX = 112
maxX = 112
minY = 15
maxY = 15

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
##########################################################

################ 流场文件指定 ###########################
flowfield_file = 'D:\\scs_model\\scs_model_only_uv_in_rho_grid\\scs_his1.nc'  #南海
#flowfield_file = 'D:\\bhd_sea_model\\bhd3_refine_undone\\out_his.nc'  #东海

content = Dataset(flowfield_file)

mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
mask_rho = content.variables['mask_rho'][:].data
lon_rho = content.variables['lon_rho'][:].data
lat_rho = content.variables['lat_rho'][:].data
#timeseries = content.variables['ocean_time'][:].data  #0/3600/7200... 秒
lonmin = np.min(lon_rho)
lonmax = np.max(lon_rho)
latmin = np.min(lat_rho)
latmax = np.max(lat_rho)
lonarray = lon_rho[0,:]
latarray = lat_rho[:,0]
len_xi = len(lonarray)
len_eta = len(latarray)  #grid size
########################################################


# V(132.00051623080796,31.73222669809303,0)

# 返回任意位置，任意时间的速度;  单位：度，度，秒； 时间从0开始;
def V(lon,lat,time): #
    # 将经纬度转化为0~229,0~239之间的一个小数坐标
    #lon,lat,time = [110.53675635432762,15.668942270366207,234234]
    #lon,lat,time = [123,25.5,4000]
    #lon,lat,time = [120.3,34.11,4000]
    #lon,lat,time = [123.91852031920288,27.193463908687335,22646*3600+2592000*3]
    xinterval = lonarray[1] - lonarray[0]
    yinterval = latarray[1] - latarray[0]
    x = 99999
    y = 99999
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

    if x == 99999 or y == 99999: #碰到边界就被黏住
        print(lon,lat,'now arrived at edge')
        return [0,0]

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

    for point in range(8):  
        if u[point].mask:
            return [0,0]

    def distance(m1,n1,k1,m2,n2,k2):
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

    for point in range(len(d)):  # 卡在了整数网格上，分母为零的情况
        if abs(d[point]) < 1e-15:
            return [1*u[point],1*v[point]]
        
    def func_d(distance_list):
        fd = []
        for point in range(len(distance_list)):
            fd.append(1/distance_list[point])
        return fd
    fd_list = func_d(d)

    def weight(fd_list):
        weight = []
        for point in range(len(fd_list)):
            weight.append((fd_list[point])/(sum(fd_list)))
        return weight
    weight = weight(fd_list)

    for point in range(8):
        u[point] = 1*u[point].data #乘1改变类型
        v[point] = 1*v[point].data

    u_cross_weight = [i*j for i,j in zip(u,weight)]
    v_cross_weight = [i*j for i,j in zip(v,weight)]
    u = sum(u_cross_weight)
    v = sum(v_cross_weight)

    return [u,v]

# 某位置，某时刻，返回一个时间步长后的状态 
# 简单欧拉法   
def jump(position):
    #position = [112.98144484693191,10.994775553461952,234234]
    lon_0 = position[0]
    lat_0 = position[1]
    time_0 = position[2]
    [u,v] = V(lon_0,lat_0,time_0)
    if u != 0 and v != 0:
        #u = u + np.random.normal(0,random_speed_sigma) #正态随机扩散，0为速度均值， 0.05为标准差
        #v = v + np.random.normal(0,random_speed_sigma)
        #u = u + random.uniform(-random_speed_sigma*2,random_speed_sigma*2) #随机扩散，0为速度均值， 0.05为标准差
        #v = v + random.uniform(-random_speed_sigma*2,random_speed_sigma*2)
        lon = lon_0 + (360/(6371000*2*3.14159*math.cos(lat_0*(math.pi/180))))*dt*u  #默认弧度制
        lat = lat_0 + (360/(6371000*2*3.14159))*dt*v
    else:  #触岸黏附
        lon = lon_0
        lat = lat_0
    return [lon,lat,time_0+dt]



# 加载粒子初始位置时刻信息
ini_particle_matrix = np.loadtxt('particles_initial.txt')

dt = 360  #粒子运动步长, 单位秒
random_speed_sigma = 0.06 #随机速度的标准差（正态分布模式） / 波动范围的四分之一（几何概型）
length = ini_particle_matrix.shape  #质点个数，由于随机扩散的存在，故可以将所有点放置在同一个初始位置
length = length[0]
N = length  # 读取质点总数

total_time = 86400*90 #总的模拟时间, 86400秒为一天
timestep_numbers = int(total_time/dt) # 时间步长数目, 后面要用这个变量
#ini_particle = [109,20.45,0]  #北部湾内部
#ini_particle = [121.5,33,0] #苏北
#ini_particle = [123,25.5,0] #台湾附近，靠近黑潮
#ini_particle = [125,30,0] #东海中心
#ini_particle = [112,14,3600] #南海中心偏左
#ini_particle = [111,18,3600] #海南岛东南角陆架
#ini_particle = [115,15,3600] #挑了一个流场比较复杂的地方


def storage_filename(x):
    return 'particle_'+str(x+1)+'.dat'

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
        #print('该质点运动了',count,'步/',int((count*dt)/(3600*24)),'天, 耗时',(t2-t1),'秒','\n')

def reconstruct():
    print('正在重新写结果')
    interval = 900  #10个时间步长写入一次全场数据, 即3600s(1小时)

    # 读取数据, 可能会爆内存
    a = np.zeros([(timestep_numbers//interval),N,2]) #维度顺序: 时间,粒子编号,经或纬
    for i in range(N):
        each_particle_data = np.loadtxt('particle_'+str(i+1)+'.dat')
        a[:,i,:] = each_particle_data[0:-1:interval,:]
    print('all data loaded in,now reconstructing...')

    for each_timestep in range(timestep_numbers//interval):
        print(each_timestep)
        with open(str(each_timestep+1)+'.dat','w',encoding='utf-8') as f:
            for each_particle in range(N):
                f.write(str(a[each_timestep,each_particle,0])+' '+str(a[each_timestep,each_particle,1])+'\n')
        f.close()
    print('重写完成')


if __name__ == '__main__':
    #按粒子, 把任务分配到八个线程计算
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



