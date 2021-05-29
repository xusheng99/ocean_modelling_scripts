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

def move(sourcefile,targetdir):
    pwd = os.getcwd()
    os.chdir(targetdir)
    filelist = os.listdir(targetdir)
    if sourcefile in filelist:
        os.remove(sourcefile)
    os.chdir(pwd)
    shutil.move(sourcefile,targetdir)
    print(str(sourcefile),'has been moved to',str(targetdir))

def point_uv2txt(ncfile,x,y,unamesuffix='_u_his_point.txt',vnamesuffix='_v_his_point.txt',dirname='C:\\Users\\XUSHENG\\Desktop\\作图用文件\\'):
    x,y = coordinatetrans(x,y)
    content = Dataset(ncfile)
    u = content.variables['u_eastward'][:,-1,y,x]
    v = content.variables['v_northward'][:,-1,y,x]
    uname = '['+str(x)+','+str(y)+']'+unamesuffix
    vname = '['+str(x)+','+str(y)+']'+vnamesuffix
    np.savetxt(uname,u)
    np.savetxt(vname,v)
    move(uname,dirname)
    move(vname,dirname)

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


point_uv2txt(hisname,115.5,20.5)
