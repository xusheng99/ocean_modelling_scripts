#! encoding='utf-8'
####    脚本说明    ####
#指定时间starttime和endtime，时间单位为"hours since201901010000"
#指定经纬度lon和lat,可以使用linspace创建一个坐标矩阵
#指定z分层深度序列，通过zarray
#主函数的参数顺序为def time_depth(lon,lat,,starttime=3624,endtime=4344,z_sequence=zarray):
#运行该脚本即在workdir下得到每个地点时间序列*深度序列的二维数组，以txt格式存储


from netCDF4 import Dataset
import os
import numpy as np 
import math
import shutil
#----------------函数中用到的参数----------------
    #经纬度范围
lon1 = 105
lon2 = 125
lat1 = 5
lat2 = 25

    #sigma分层参数
theta_s,theta_b,hc,N = [4.5,1.5,5,40]
sc = [theta_s,theta_b,hc,N]

    #指定的z分层每层的深度
zarray = [2,3,5,10,15,20,30,40,50,60,70,80,90,100,120,150,200,300,400,500,1000,2000]

    #工作目录和绘图目录
workdir = 'C:\\Users\\XUSHENG\\Desktop\\read_roms_result\\'    #将某些较小的nc文件直接放到这里
plotdir = 'C:\\Users\\XUSHENG\\Desktop\\作图用文件\\'    #读取的数据移动到这里供matlab使用
os.chdir(workdir)
    
    #netcdf文件的位置
hisname = 'C:\\Users\\XUSHENG\\Desktop\\SCS\\1.nc'    ## 不要打错字，好好检查一下路径
gridname = 'C:\\Users\\XUSHENG\\Desktop\\read_roms_result\\roms_grid.nc'

    #掩码和深度
grid = Dataset(gridname)
maskr = grid.variables['maskr'][:]
lon_length = maskr.shape[1]-1   #从mask读取网格的大小！
lat_length = maskr.shape[0]-1

    #sigma层数对应深度的三维数组
def zmaker(gridfilename,sc,zeta=0):
    #顺序给出网格文件,theta_s,theta_b,hc,N；返回指定位置指定sigma层数的深度
    nc=Dataset(gridfilename)
    h=nc.variables['h'][:]
    hmin=np.nanmin(h)
    nc.close()
    theta_s,theta_b,Tcline,N=sc
    hc=Tcline
    if hc>hmin and np.abs(hc-hmin)>1e-10:
        raise ValueError ("Critical depth exceeds minimum")
    M=h.shape
    sc=(np.arange(1.,N+0.1)-N-0.5)/N    #使用1.是为了指定浮点型的1.0数字，便于精确的除法运算
    #sc创建了一个从-9.8到-0.2的数量为25个的等差数列
    Cs=(1-theta_b)*np.sinh(theta_s*sc)/np.sinh(theta_s)+ \
        theta_b*(-0.5+0.5*np.tanh(theta_s*(sc+0.5))/np.tanh(0.5*theta_s))
    z=np.zeros([N,M[0],M[1]])    #创建等大小的空数组
    for k in range(N):    #对每一层来说
        z0=hc*sc[k]+(h-hc)*Cs[k]#
        z[k]=z0+zeta*(1.+z0/h)
    #搞出我想要的**单增恒正**的z函数
    z = z.transpose(2,1,0)
    z = z*(-1)    #把z翻成正值
    ztemp = np.empty(z.shape)
    for k in range(N):
        ztemp[:,:,N-k-1] = z[:,:,k]
    z = ztemp    #把z上下翻转一遍，使之成为海表为0层，海底为24层
    return z
z = zmaker(gridname,sc)

def sigmalevelfinder(x,y,sd):#sd为深度，正值，给定深度，返回相应的sigma层数
    for k in range(N):
        #先解决极端情况
        if sd < z[x,y,0]:
            re = 0
            break 
        if sd > z[x,y,N-1]:
            re = N-1
            break
        #再解决一般情况
        if z[x,y,k] < sd < z[x,y,k+1]:
            re = [k,k+1]
            break
    return re

def sigmalevelfinder2(x,y,sd):
    for k in range(N):
        if sd < z[x,y,0]:
            re = 0
            break 
        if sd > z[x,y,N-1]:
            re = N-1
            break
        if z[x,y,k] < sd < z[x,y,k+1]:
            if (sd - z[x,y,k]) < (z[x,y,k+1]-sd):#离k比较近
                re = k
                break
            else:
                re = k+1
                break
    return re

def coordinatetrans(x,y):#指定经纬度，自动选取对应最合适的网格点
    def cord_trans(x,y):
        #第一步，将给定形式的经纬度转化成标准的小数表示的经纬度
        a = np.arange(0.1,1.1,0.15)
        npfloattype = type(a[3])
        if type(x) == type(1) or type(x) == npfloattype or type(x) == type(1.1):    #整数或numpy浮点类型
            point_lon = x
        elif len(x) == 2:
            point_lon = x[0]+(x[1]/60)
        elif len(x) == 3:
            point_lon = x[0]+(x[1]/60)+(x[2]/3600)
        else:
            raise ValueError ('invalid input!')
        if type(y) == type(1) or type(y) == npfloattype or type(x) == type(1.1):
            point_lat = y
        elif len(y) == 2:
            point_lat = y[0]+(y[1]/60)
        elif len(y) == 3:
            point_lat = y[0]+(y[1]/60)+(y[2]/3600)
        else:
            raise ValueError ('invalid input!')
        x = (((point_lon-lon1)/(lon2-lon1))*lon_length)
        y = (((point_lat-lat1)/(lat2-lat1))*lat_length)
        return x,y
    def distance(p1,p2):
        #定义计算两点距离的函数，以两个二元列表的形式输入
        x = p1[0]
        y = p1[1]
        m = p2[0]
        n = p2[1]
        return math.sqrt((x-m)*(x-m)+(y-n)*(y-n))  
    def namegen(i):
        return 'P'+str(i)
    def trans_int(x,y):
        xtemp = round(x,0)
        ytemp = round(y,0)
        xtemp = int(xtemp)    #这里的int函数必不可少
        ytemp = int(ytemp)
        if maskr[ytemp,xtemp] == 1 : ###!!!!!! 必须先将ytemp和xtemp转化成整数才能作为数组的索引，1.0==1，但在数组索引里就不同。
        #切记nc文件里读出来的变量都是先lat后lon的！！！！！
            point_coordinate_int = [xtemp,ytemp]
        else:
            P1=[xtemp-1,ytemp-1]
            P2=[xtemp,ytemp-1]
            P3=[xtemp+1,ytemp-1]
            P4=[xtemp-1,ytemp]
            P5=[x,y]
            P6=[xtemp+1,ytemp]
            P7=[xtemp-1,ytemp+1]
            P8=[xtemp,ytemp+1]
            P9=[xtemp+1,ytemp+1]
            distance_list = []
            for i in [1,2,3,4,6,7,8,9]:
                if maskr[locals()[namegen(i)][1],locals()[namegen(i)][0]] == 1:
                    if i%2 != 0:
                        distance_list.append((distance(locals()[namegen(i)],P5)/1.414))
                    else:
                        distance_list.append(distance(locals()[namegen(i)],P5))
                else:
                    distance_list.append(1000)
            minvalue = np.min(distance_list)
            num = distance_list.index(minvalue)
            if num == 0:
                point_coordinate_int = P1
                print("取了左下角点作为海水点近似")
            elif num == 1:
                point_coordinate_int = P2
                print("取了下方点作为海水点近似")
            elif num ==2:
                point_coordinate_int = P3
                print("取了右下角点作为海水点近似")
            elif num == 3:
                point_coordinate_int = P4
                print("取了左边点作为海水点近似")
            elif num ==4:
                point_coordinate_int = P6
                print("取了右边点作为海水点近似")
            elif num == 5:
                point_coordinate_int = P7
                print("取了左上角点作为海水点近似")
            elif num == 6:
                point_coordinate_int = P8
                print("取了上方点作为海水点近似")
            elif num == 7:
                point_coordinate_int = P9
                print("取了右上角点作为海水点近似")
            else:
                print('务必选一个在海里的位置！')
        return point_coordinate_int[0],point_coordinate_int[1]
    x,y = cord_trans(x,y)
    x0=x
    y0=y
    print('格点精确值：',x0,y0)
    x,y = trans_int(x,y)
    print('最终近似格点：',x,y)
    return x,y






def time_depth_z(lon,lat,starttime=3624,endtime=4344,z_sequence=zarray):
    #读取
    x,y = coordinatetrans(lon,lat)
    content = Dataset(hisname)

    salt = content.variables['salt'][3624:4344,:,y,x].data    #直接从nc中读出的均为先lat后lon
    temp = content.variables['temp'][3624:4344,:,y,x].data
    u = content.variables['u'][3624:4344,:,y,x].data
    v = content.variables['v'][3624:4344,:,y,x].data

    #清洁
    def clean(x):
        x[abs(x>50)] = np.nan
    clean(u)
    clean(v)
    clean(salt)
    clean(temp)

    #表底反序
    salt1 = np.zeros((salt.shape))
    temp1 = np.zeros((temp.shape))
    u1 = np.zeros((u.shape))
    v1 = np.zeros((v.shape))
    for k in range(salt.shape[1]):
        salt1[:,k] = salt[:,salt.shape[1]-k]
        temp1[:,k] = temp[:,temp.shape[1]-k]
        u1[:,k] = u[:,u.shape[1]-k]
        v1[:,k] = v[:,v.shape[1]-k]
    salt = salt1
    temp = temp1
    u = u1
    v = v1


    #创建空数组用于存储
    saltz = np.zeros([salt.shape[0],len(z_sequence)])
    saltz[saltz == 0] = np.nan
    tempz = np.ones([temp.shape[0],len(z_sequence)])
    saltz[saltz == 1] = np.nan
    uz = np.zeros([u.shape[0],len(z_sequence)])
    vz = np.zeros([v.shape[0],len(z_sequence)])
        #这里不要图方便，写后面的几个tempz,uz,vz等于saltz，如果这样写了，这几个变量将指向同一片内存区域，后面的会覆盖前面的

    #线性插值计算并写入
    for t in range(saltz.shape[0]):
        for k in range(len(z_sequence)):
            
            sd = z_sequence[k]
            if sd > depth[y,x]:
                saltz[t,k] = np.nan
                tempz[t,k] = np.nan
                uz[t,k] = np.nan
                vz[t,k] = np.nan
            else:    #使用内置的插值函数要更简洁且方便
                xp = z[x,y,:]

                saltp = salt[t,:]
                saltz[t,k] = np.interp(sd,xp,saltp)
                tempp = temp[t,:]
                tempz[t,k] = np.interp(sd,xp,tempp)
                up = u[t,:]
                uz[t,k] = np.interp(sd,xp,up)
                vp = v[t,:]
                vz[t,k] = np.interp(sd,xp,vp)
            

            '''
            sd = z_sequence[k]
            if sd > depth[y,x]:
                saltz[t,k] = np.nan
                tempz[t,k] = np.nan
                uz[t,k] = np.nan
                vz[t,k] = np.nan
            else:
                up_down_level = sigmalevelfinder(x,y,sd)
                re = sigmalevelfinder2(x,y,sd)
                if type(up_down_level) == type([1,2]):
                    m = up_down_level[0]
                    n = up_down_level[1]
                    saltz[t,k] = salt[t,m]+((salt[t,n]-salt[t,m])/\
                        (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                    tempz[t,k] = temp[t,m]+((temp[t,n]-temp[t,m])/\
                        (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                    uz[t,k] = u[t,m]+((u[t,n]-u[t,m])/\
                        (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                    vz[t,k] = v[t,m]+((v[t,n]-v[t,m])/\
                        (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                else: #如果是表底层
                    if up_down_level == 0:  #表层  
                        m = up_down_level
                        n = m + 1
                        saltz[t,k] = salt[t,m]+((salt[t,n]-salt[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                        tempz[t,k] = temp[t,m]+((temp[t,n]-temp[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                        uz[t,k] = u[t,m]+((u[t,n]-u[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                        vz[t,k] = v[t,m]+((v[t,n]-v[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                    else:    #底层
                        n = up_down_level
                        m = n - 1
                        saltz[t,k] = salt[t,m]+((salt[t,n]-salt[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                        tempz[t,k] = temp[t,m]+((temp[t,n]-temp[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                        uz[t,k] = u[t,m]+((u[t,n]-u[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
                        vz[t,k] = v[t,m]+((v[t,n]-v[t,m])/\
                            (z[x,y,n]-z[x,y,m]))*(sd-z[x,y,m])
           '''
    return saltz,tempz,uz,vz


def rrr(x,y):
    
x,y = coordinatetrans(115.5,20.5)
content = Dataset(hisname)
u = content.variables['u_eastward'][:,:,y,x].data
v = content.variables['v_northward'][:,:,y,x].data
temp = content.variables['temp'][:,:,y,x].data
salt = content.variables['salt'][:,:,y,x].data

#inverse
u1 = np.zeros(u.shape)
v1 = np.zeros(v.shape)
temp1 = np.zeros(temp.shape)
salt1 = np.zeros(salt.shape)
for k in range(u1.shape[1]):
    u1[:,k] = u[:,u.shape[1]-k]
    v1[:,k] = v[:,v.shape[1]-k]
    temp1[:,k] = temp[:,temp.shape[1]-k]
    salt1[:,k] = salt[:,salt.shape[1]-k]
u = u1
v = v1
salt = salt1
temp = temp1






lonlist = np.linspace(115.5,116.0,6)
latlist = np.linspace(20.5,21.0,6)

def r(lon,lat):
    salt,temp,u,v = time_depth_z(lon,lat)
    saltname = str(lat)+'N_'+str(lon)+'E'+'_盐度.txt'
    tempname = str(lat)+'N_'+str(lon)+'E'+'_温度.txt'
    uname = str(lat)+'N_'+str(lon)+'E'+'_u.txt'
    vname = str(lat)+'N_'+str(lon)+'E'+'_v.txt'
    np.savetxt(saltname,salt)
    np.savetxt(tempname,temp)
    np.savetxt(uname,u)
    np.savetxt(vname,v)


lonlist = np.linspace(115.5,116.0,6)
latlist = np.linspace(20.5,21.0,6)

for i in lonlist:
    for j in latlist:
        r(i,j)
        print(i,j,'done')

'''
lon = 120.5
lat = 20.5
r(lon,lat)
'''
'''
x,y = coordinatetrans(115.5,20.5)
lonlist = np.linspace(115.5,116.0,6)
latlist = np.linspace(20.5,21.0,6)
for lon in lonlist:
    for lat in latlist:
   '''     