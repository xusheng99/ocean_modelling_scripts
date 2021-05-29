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
zarray = [2,3,5,10,15,20,30,40,50,60,70,80,90,100,120,150,200,300,400,500,600,700,800,900,1000]

    #工作目录和绘图目录
workdir = 'C:\\Users\\XUSHENG\\Desktop\\read_roms_result\\'    #将某些较小的nc文件直接放到这里
plotdir = 'C:\\Users\\XUSHENG\\Desktop\\作图用文件\\'    #读取的数据移动到这里供matlab使用
os.chdir(workdir)
    
    #netcdf文件的位置
hisname = 'C:\\Users\\XUSHENG\\Desktop\\SCS\\1.nc'    ## 不要打错字，好好检查一下路径

    #掩码和深度
lon_length = 360   
lat_length = 360

nc = Dataset(hisname)
depth = nc.variables['h'][:]

#---------------以上为基本参数----------------------------



def move(sourcefile,targetdir):
    pwd = os.getcwd()
    os.chdir(targetdir)
    filelist = os.listdir(targetdir)
    if sourcefile in filelist:
        os.remove(sourcefile)
    os.chdir(pwd)
    shutil.move(sourcefile,targetdir)
    print(str(sourcefile),'has been moved to',str(targetdir))

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
        if 1 == 1 : ###!!!!!! 必须先将ytemp和xtemp转化成整数才能作为数组的索引，1.0==1，但在数组索引里就不同。
        #切记nc文件里读出来的变量都是先lat后lon的！！！！！
            point_coordinate_int = [xtemp,ytemp]
            '''
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
            '''
        return point_coordinate_int[0],point_coordinate_int[1]
    x,y = cord_trans(x,y)
    x0=x
    y0=y
    print('格点精确值：',x0,y0)
    x,y = trans_int(x,y)
    print('最终近似格点：',x,y)
    return x,y



def rrr(lon,lat):    #sigma版
    x,y = coordinatetrans(lon,lat)
    content = Dataset(hisname)
    u = content.variables['u_eastward'][:,:,y,x].data
    v = content.variables['v_northward'][:,:,y,x].data
    temp = content.variables['temp'][:,:,y,x].data
    salt = content.variables['salt'][:,:,y,x].data

    z = np.loadtxt(zname(lon,lat))

    #inverse
    u1 = np.zeros(u.shape)
    v1 = np.zeros(v.shape)
    temp1 = np.zeros(temp.shape)
    salt1 = np.zeros(salt.shape)
    for k in range(u1.shape[1]):
        u1[:,k] = u[:,u.shape[1]-1-k]
        v1[:,k] = v[:,v.shape[1]-1-k]
        temp1[:,k] = temp[:,temp.shape[1]-1-k]
        salt1[:,k] = salt[:,salt.shape[1]-1-k]
    u = u1
    v = v1
    salt = salt1
    temp = temp1

    uname = 'u_'+str(lat)+'N_'+str(lon)+'E.txt'
    vname = 'v_'+str(lat)+'N_'+str(lon)+'E.txt'
    tempname = 'temp_'+str(lat)+'N_'+str(lon)+'E.txt'
    saltname = 'salt_'+str(lat)+'N_'+str(lon)+'E.txt'
    depthinfoname = str(lat)+'N_'+str(lon)+'E'+'sigma分层信息.txt'
    np.savetxt(depthinfoname,z[:])
    np.savetxt(uname,u)
    np.savetxt(vname,v)
    np.savetxt(tempname,temp)
    np.savetxt(saltname,salt)

    #单独归类文件夹
    dirname = str(lat)+'E_'+str(lon)+'N'
    os.mkdir(dirname)
    move(uname,workdir+dirname)
    move(vname,workdir+dirname)
    move(saltname,workdir+dirname)
    move(tempname,workdir+dirname)
    move(depthinfoname,workdir+dirname)

def rrrz(lon,lat,z_sequence=zarray):    #z分层数据
    x,y = coordinatetrans(lon,lat)
    content = Dataset(hisname)
    u = content.variables['u_eastward'][:,:,y,x].data
    v = content.variables['v_northward'][:,:,y,x].data
    temp = content.variables['temp'][:,:,y,x].data
    salt = content.variables['salt'][:,:,y,x].data

    z = np.loadtxt(zname(lon,lat))

    #inverse
    u1 = np.zeros(u.shape)
    v1 = np.zeros(v.shape)
    temp1 = np.zeros(temp.shape)
    salt1 = np.zeros(salt.shape)
    for k in range(u1.shape[1]):
        u1[:,k] = u[:,u.shape[1]-1-k]
        v1[:,k] = v[:,v.shape[1]-1-k]
        temp1[:,k] = temp[:,temp.shape[1]-1-k]
        salt1[:,k] = salt[:,salt.shape[1]-1-k]
    u = u1
    v = v1
    salt = salt1
    temp = temp1
    
    #插值
    saltz = np.zeros([salt.shape[0],len(z_sequence)])
    saltz[saltz == 0] = np.nan
    tempz = np.ones([temp.shape[0],len(z_sequence)])
    saltz[saltz == 1] = np.nan
    uz = np.zeros([u.shape[0],len(z_sequence)])
    vz = np.zeros([v.shape[0],len(z_sequence)])

    for t in range(saltz.shape[0]):
        for k in range(len(z_sequence)):
            
            sd = z_sequence[k]
            if sd > depth[y,x]:    #判断水深
                saltz[t,k] = np.nan
                tempz[t,k] = np.nan
                uz[t,k] = np.nan
                vz[t,k] = np.nan
            else:    #使用内置的插值函数要更简洁且方便
                xp = z[:]

                saltp = salt[t,:]
                saltz[t,k] = np.interp(sd,xp,saltp)
                tempp = temp[t,:]
                tempz[t,k] = np.interp(sd,xp,tempp)
                up = u[t,:]
                uz[t,k] = np.interp(sd,xp,up)
                vp = v[t,:]
                vz[t,k] = np.interp(sd,xp,vp)

    #按点位构造文件名
    uname = 'uz_'+str(lat)+'N_'+str(lon)+'E.txt'
    vname = 'vz_'+str(lat)+'N_'+str(lon)+'E.txt'
    tempname = 'tempz_'+str(lat)+'N_'+str(lon)+'E.txt'
    saltname = 'saltz_'+str(lat)+'N_'+str(lon)+'E.txt'
    
    np.savetxt(uname,uz)
    np.savetxt(vname,vz)
    np.savetxt(tempname,tempz)
    np.savetxt(saltname,saltz)

    #按点位单独归类文件夹
    dirname = str(lat)+'E_'+str(lon)+'N'
    os.mkdir(dirname)
    move(uname,workdir+dirname)
    move(vname,workdir+dirname)
    move(saltname,workdir+dirname)
    move(tempname,workdir+dirname)

    print(lon,lat,'done','\n')

def zname(lon,lat):    #因为python中sigma分层的计算公式有点问题，所以从matlab中读出相应的列表在这里使用，靠名字区分点位
    zname = str(lon)+'E_'+str(lat)+'N.txt'
    return zname


lonlist = np.linspace(115.5,116.0,6)
latlist = np.linspace(20.5,21.0,6)

for i in lonlist:
    for j in latlist:
        rrrz(i,j)


with open('z分层信息.txt','w',encoding='utf-8') as f:
    for i in range(len(zarray)):
        f.write(str(i+1)+'层'+' ')
        f.write(str(zarray[i])+'米'+'\n')
f.close()

move('z分层信息.txt',workdir+'z分层数据')


