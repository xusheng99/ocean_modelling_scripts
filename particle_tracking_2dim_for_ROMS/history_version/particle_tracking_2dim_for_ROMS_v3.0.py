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

#v3.0 改动说明
#尝试使用多进程来加快速度

os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\')

'''
#读取配置文件
with open('config.txt','r',encoding='utf-8') as f:
    config = f.read()
f.close()
'''


################################################
#在命令行中键入 python ./particle_tarcking_2dim_for_ROMS.py N
################################################


flowfield_file = 'D:\\scs_model\\scs_model_only_uv_in_rho_grid\\scs_his1.nc'  #南海
#flowfield_file = 'D:\\bhd_sea_model\\bhd3_refine_undone\\out_his.nc'  #东海

content = Dataset(flowfield_file)

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
len_eta = len(latarray)  #网格大小

'''
a = round(((18-latmin)/(latmax-latmin))*238)
b =round(((117-lonmin)/(lonmax - lonmin))*238)
c = int(234234/3600)
uu = content.variables['u_eastward'][0:100,-1,:,:].data
vv = content.variables['v_northward'][0:100,-1,:,:].data
'''

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

    '''
    #老算法，简单的线性内插
    # 此处可考虑用其他算法增加精度
    u1 = content.variables['u_eastward'][time1,-1,y1,x1].data
    v1 = content.variables['v_northward'][time1,-1,y1,x1].data
    u_round_x = (content.variables['u_eastward'][time1,-1,y1,x2] - content.variables['u_eastward'][time1,-1,y1,x1])
    u_round_y = (content.variables['u_eastward'][time1,-1,y2,x1] - content.variables['u_eastward'][time1,-1,y1,x1])
    u_round_t = (content.variables['u_eastward'][time2,-1,y1,x1] - content.variables['u_eastward'][time1,-1,y1,x1])
    u = u1 + u_round_x*delta_x + u_round_y*delta_y + u_round_t*delta_time
    v_round_x = (content.variables['v_northward'][time1,-1,y1,x2] - content.variables['v_northward'][time1,-1,y1,x1])
    v_round_y = (content.variables['v_northward'][time1,-1,y2,x1] - content.variables['v_northward'][time1,-1,y1,x1])
    v_round_t = (content.variables['v_northward'][time2,-1,y1,x1] - content.variables['v_northward'][time1,-1,y1,x1])
    v = v1 + v_round_x*delta_x + v_round_y*delta_y + v_round_t*delta_time
    '''
    return [u,v]

# 某位置，某时刻，返回一个时间步长后的状态    
def jump(position):
    #position = [112.98144484693191,10.994775553461952,234234]
    lon_0 = position[0]
    lat_0 = position[1]
    time_0 = position[2]
    [u,v] = V(lon_0,lat_0,time_0)
    if u != 0 and v != 0:
        u = u + np.random.normal(0,random_speed_sigma) #正态随机扩散，0为速度均值， 0.05为标准差
        v = v + np.random.normal(0,random_speed_sigma)
        #u = u + random.uniform(-random_speed_sigma*2,random_speed_sigma*2) #随机扩散，0为速度均值， 0.05为标准差
        #v = v + random.uniform(-random_speed_sigma*2,random_speed_sigma*2)
        lon = lon_0 + (360/(6371000*2*3.14159*math.cos(lat_0*(math.pi/180))))*dt*u  #默认弧度制
        lat = lat_0 + (360/(6371000*2*3.14159))*dt*v
    else:  #触岸黏附
        lon = lon_0
        lat = lat_0
    return [lon,lat,time_0+dt]

'''
def initial()
with open('0config.txt','r',encoding='utf-8') as f:
    config = f.read()
f.close()
N = re.search(r'N\s=\s([0-9]+)',config).group[1]
'''

dt = 360  #单位秒
random_speed_sigma = 0.06 #随机速度的标准差（正态） / 波动范围的四分之一（线性）
N = int(sys.argv[1])  #质点个数，由于随机扩散的存在，故可以将所有点放置在同一个初始位置
#N = 30
total_time = 86400*30 #总的模拟时间
#ini_particle = [109,20.45,0]  #北部湾内部
#ini_particle = [121.5,33,0] #苏北
#ini_particle = [123,25.5,0] #台湾附近，靠近黑潮
#ini_particle = [125,30,0] #东海中心
#ini_particle = [112,14,3600] #南海中心偏左
#ini_particle = [111,18,3600] #海南岛东南角陆架
ini_particle = [115,15,3600] #挑了一个流场比较复杂的地方

def storage_filename(x):
    return str(x+1)+'.dat'

def release(basenum,NN):
    for particle_num in range(basenum,basenum+NN):
        t1 = time.perf_counter()
        with open(storage_filename(particle_num),'w',encoding='utf-8') as f:
            print('now operating: ',storage_filename(particle_num))
            condition = ini_particle
            count = 0
            while condition[2] < total_time+ini_particle[2]: #超过截至时间，停止运行
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



if __name__ == '__main__':
    
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
    print('Child process start.')
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
    T2 = time.perf_counter()
    print()
    print('共计耗时: ',T2-T1,'秒')


'''
d,e,f = [130,27,345]
V(d,e,f)
V2(d,e,f)
jump([121.70314514155001,33.1160913448275,12373200]) #这个时间超了，会奥错

mat = np.load('rlt_gene_features.npy-layer-3-train.npy')
scipy.io.savemat('gene_features.mat', {'gene_features': mat})

#这个地理坐标报错了，撞到了陆地？但似乎这里也没有陆地啊？
jump([110.53675635432762,15.668942270366207,234234])
V(110.53675635432762,15.668942270366207,234234)  
#原因找到了，是因为cos使用弧度制，输入了角度，由于cos在分母，某些值会产生巨大的波动，直接跳出南海
jump([112.98144484693191,10.994775553461952,234234])
V(112.98144484693191,10.994775553461952,234234)
s = np.random.normal(0,0.05)
s
'''
