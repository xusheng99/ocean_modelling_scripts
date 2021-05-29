from netCDF4 import Dataset
import os
import numpy as np 
import math
import shutil
import matlab
import matlab.engine

#----------------函数中用到的参数----------------   
    #经纬度范围
lon1 = 105
lon2 = 125
lat1 = 5
lat2 = 25

    #sigma分层参数
theta_s,theta_b,hc,N = [7,2,15,40]
sc = [theta_s,theta_b,hc,N]

    #指定的z分层每层的深度
zarray = [3,5,10,15,20,40,50,100,150,200,300,500,1000,1500,2000,3000]

    #工作目录和绘图目录
workdir = 'C:\\Users\\XUSHENG\\Desktop\\read_roms_result\\'    #将某些较小的nc文件直接放到这里
plotdir = 'C:\\Users\\XUSHENG\\Desktop\\作图用文件\\'    #读取的数据移动到这里供matlab使用
os.chdir(workdir)
    
    #netcdf文件的位置
ininame = 'roms_ini.nc'
#因为这个文件一般比较大，不建议搬运到本文件夹后再处理
hisname = 'd:\\application\\SCS\\his_0003.nc'    ## 不要打错字，好好检查一下路径
gridname = 'd:\\application\\SCS\\DATA\\Grd.nc'
bryname = 'roms_bry.nc'
grid = Dataset(gridname)
his = Dataset(hisname)

    #掩码和深度
maskr = grid.variables['mask_rho'][:]
masku = grid.variables['mask_u'][:]
maskv = grid.variables['mask_v'][:]
depth = grid.variables['h'][:]
depth = depth*maskr 
lon_length = maskr.shape[0]-1    #从mask读取网格的大小！
lat_length = maskr.shape[1]-1
depth.fill_value = 0
'''
eng = matlab.engine.start_matlab()
def zmaker(gridfilename,sc,zeta=0):
    nc = Dataset(gridfilename)
    h = nc.variables['h'][:]
    h.fill_value = 0 
    h = h.data
    h = h.tolist()
    theta_s,theta_b,Tcline,N=sc
    z = eng.zlevs(h,zeta,theta_s,theta_b,Tcline,N,'r',2)
    return z
    #sigma层数对应深度的三维数组
'''
def zmaker(gridfilename,sc,zeta=0):
    #顺序给出网格文件,theta_s,theta_b,hc,N；返回指定位置指定sigma层数的深度
    nc=Dataset(gridfilename)
    h=nc.variables['h'][:]
    hmin=np.nanmin(h)
    nc.close()
    theta_s,theta_b,Tcline,N=sc
    hc=Tcline
    '''
    if hc>hmin and np.abs(hc-hmin)>1e-10:
        raise ValueError ("Critical depth exceeds minimum")
    '''
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
z = zmaker(gridname,sc)    #重要的sigma层数-实际深度的三维数组

#-------------------以上参数不可缺失--------------------------------------------------




#--------------------一些用到的函数-----------------------------------------------
def date2time(date0):
    pass

def coordinatetrans(x,y):#指定经纬度，自动选取对应最合适的网格点
    def cord_trans(x,y):
        #第一步，将给定形式的经纬度转化成标准的小数表示的经纬度
        a = np.arange(0.1,1.1,0.15)
        npfloattype = type(a[3])
        if type(x) == type(1) or type(x) == npfloattype or type(x) == type(1.11):    #整数或numpy浮点类型
            point_lon = x
        elif len(x) == 2:
            point_lon = x[0]+(x[1]/60)
        elif len(x) == 3:
            point_lon = x[0]+(x[1]/60)+(x[2]/3600)
        else:
            raise ValueError ('invalid input!')
        if type(y) == type(1) or type(y) == npfloattype or type(x) == type(1.11):
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
            for i in [1,2,3,4,6,7,8,9]:    #生成变量名的方法!
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
            elif num == 1:
                point_coordinate_int = P2
            elif num ==2:
                point_coordinate_int = P3
            elif num == 3:
                point_coordinate_int = P4
            elif num ==4:
                point_coordinate_int = P6
            elif num == 5:
                point_coordinate_int = P7
            elif num == 6:
                point_coordinate_int = P8
            elif num == 7:
                point_coordinate_int = P9
            else:
                pass
        return point_coordinate_int[0],point_coordinate_int[1]
    x,y = cord_trans(x,y)
    x0=x
    y0=y
    print('格点精确值：',x0,y0)
    x,y = trans_int(x,y)
    print('最终近似格点：',x,y)
    return x,y

def move(sourcefile,targetdir):    #移动sourcefile至targetdir,如果已存在同名文件则覆盖
    pwd = os.getcwd()
    os.chdir(targetdir)
    filelist = os.listdir(targetdir)
    if sourcefile in filelist:
        os.remove(sourcefile)
    os.chdir(pwd)
    shutil.move(sourcefile,targetdir)
    print(str(sourcefile),'has been moved to',str(targetdir))

#按顺序指定ncfile，经度，纬度，u存储文件，v存储文件，文件保存位置；来生成全时间序列，全深度的二维数组
#默认的存储文件夹为桌面的“作图用文件”，默认的txt名为u_his_point.txt以及v_his_point.txt
def point_uv2txt(ncfile,x,y,unamesuffix='_u_his_point.txt',vnamesuffix='_v_his_point.txt',dirname='C:\\Users\\XUSHENG\\Desktop\\作图用文件\\'):
    x,y = coordinatetrans(x,y)
    content = Dataset(ncfile)
    u = content.variables['u_eastward'][:,24,y,x]
    v = content.variables['v_northward'][:,24,y,x]
    uname = '['+str(x)+','+str(y)+']'+unamesuffix
    vname = '['+str(x)+','+str(y)+']'+vnamesuffix
    np.savetxt(uname,u)
    np.savetxt(vname,v)
    move(uname,dirname)
    move(vname,dirname)

point_uv2txt(hisname,108,19)
def zetaslice(ncfile,time):#读取指定时间的全场水位数据并存储到作图文件夹
    nc = Dataset(ncfile)
    zeta = nc.variables['zeta'][time,:,:]
    zeta = zeta.data
    clean(zeta)
    zeta = zeta.transpose(1,0)
    a,b = zeta.shape
    with open('zeta.txt','w',encoding='utf-8') as f:
        for i in range(a):
            for j in range(b):
                f.write(str(i)+' ')
                f.write(str(j)+' ')
                f.write(str(zeta[i,j]))
                f.write('\n')
            print(i)
        f.close()
        print('all done!')
    move('zeta.txt',plotdir)
    print('zeta.txt 已经存储到指定文件夹')
    return zeta

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
    return re  #返回两个数, 用于Linear method

def sigmalevelfinder2(x,y,sd):  #返回最靠近的层数, 用于nearest method
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

def clean(x):  #清洗数据, 删除异常点
    x[abs(x)>50] = np.nan

#指定深度的全部物理量的切片
def depth_slice(ncfile,time,sd,interp_method='linear'):
    #method可选 linear,nearest    线性插值或最近邻插值
    #读取u_his,v_his中的value，不含掩码
    pass
    content = Dataset(ncfile)
    u = content.variables['u_eastward'][time,:,:,:].data    #！！！务必要写.data,否则读出的是带mask的数据，会极大拖慢插值运算
    v = content.variables['v_northward'][time,:,:,:].data    #masku是二维的，u_his是三维的，无法直接相乘！
    salt = content.variables['salt'][time,:,:,:].data
    temp = content.variables['temp'][time,:,:,:].data
    N = salt.shape[0]    #层数
    for k in range(u.shape[0]):    #这里不要写成u_his.shape(0),会报元组不可遍历的错
        u[k,:,:] = u[k,:,:]*maskr
    for k in range(v.shape[0]):
        v[k,:,:] = v[k,:,:]*maskr    
    clean(u)    #删除异常点
    clean(v)    #比如u[16,7,190]就是一个异常大的数值，导致作图时正常数据被淹没
    clean(salt)
    clean(temp)
 
    salt = salt.transpose(2,1,0)
    temp = temp.transpose(2,1,0)
    u = u.transpose(2,1,0)
    v = v.transpose(2,1,0)

    #咸鱼翻身~
    salt1 = np.zeros((salt.shape))
    temp1 = np.zeros((temp.shape))
    u1 = np.zeros((u.shape))
    v1 = np.zeros((v.shape))
    for k in range(salt.shape[2]):
        salt1[:,:,k] = salt[:,:,N-1-k]
        temp1[:,:,k] = temp[:,:,N-1-k]
        u1[:,:,k] = u[:,:,N-1-k]
        v1[:,:,k] = v[:,:,N-1-k]
    salt = salt1
    temp = temp1
    u = u1
    v = v1

    #新建二维空数组 用于存储
    uz = np.zeros((depth.shape[1],depth.shape[0]))
    vz = np.zeros((depth.shape[1],depth.shape[0]))
    Vz = np.zeros((depth.shape[1],depth.shape[0]))
    saltz = np.zeros((depth.shape[1],depth.shape[0]))
    tempz = np.zeros((depth.shape[1],depth.shape[0]))

    if interp_method == 'linear':
        pass
        for j in range(depth.shape[0]):
            for i in range(depth.shape[1]):
                #判断我想要的深度值是不是真的在水里
                if sd > depth[j,i]: #如果在海底，则填nan值
                    uz[i,j] = np.nan
                    vz[i,j] = np.nan
                    saltz[i,j] = np.nan
                    tempz[i,j] = np.nan
                else:  #如果不在海底，则根据sigmalevelfinder找到的层数插值
                    up_down_level = sigmalevelfinder(i,j,sd)
                    if type(up_down_level) == type([1,2]):#如果非表层且非底层
                        m = up_down_level[0]
                        n = up_down_level[1]

                        uz[i,j] = u[i,j,m]+((u[i,j,n]-u[i,j,m])/\
                            (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                        vz[i,j] = v[i,j,m]+((v[i,j,n]-v[i,j,m])/\
                            (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                        saltz[i,j] = salt[i,j,m]+((salt[i,j,n]-salt[i,j,m])/\
                            (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                        tempz[i,j] = temp[i,j,m]+((temp[i,j,n]-temp[i,j,m])/\
                            (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                        Vz[i,j] = math.sqrt((uz[i,j])*(uz[i,j])+(vz[i,j])*(vz[i,j]))
                        Vz[Vz==0] = np.nan
                    else: #如果是表底层
                        if up_down_level == 0:  #表层
                            m = up_down_level
                            n = m + 1
                            uz[i,j] = u[i,j,m]+((u[i,j,n]-u[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                            vz[i,j] = v[i,j,m]+((v[i,j,n]-v[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                            saltz[i,j] = salt[i,j,m]+((salt[i,j,n]-salt[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                            tempz[i,j] = temp[i,j,m]+((temp[i,j,n]-temp[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                            Vz[i,j] = math.sqrt((uz[i,j])*(uz[i,j])+(vz[i,j])*(vz[i,j]))
                            Vz[Vz==0] = np.nan
                        else:  #底层
                            n = up_down_level
                            m = n - 1
                            uz[i,j] = u[i,j,m]+((u[i,j,n]-u[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                            vz[i,j] = v[i,j,m]+((v[i,j,n]-v[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m]) 
                            saltz[i,j] = salt[i,j,m]+((salt[i,j,n]-salt[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                            tempz[i,j] = temp[i,j,m]+((temp[i,j,n]-temp[i,j,m])/\
                                (z[i,j,n]-z[i,j,m]))*(sd-z[i,j,m])
                            Vz[i,j] = math.sqrt((uz[i,j])*(uz[i,j])+(vz[i,j])*(vz[i,j]))   
                            Vz[Vz==0] = np.nan

    if interp_method == 'nearest':
        pass
        for j in range(depth.shape[0]):
            for i in range(depth.shape[1]):
                #判断我想要的深度值是不是真的在水里
                if sd > depth[j,i]: #如果在海底，则填nan值
                    uz[i,j] = np.nan
                    vz[i,j] = np.nan
                    saltz[i,j] = np.nan
                    tempz[i,j] = np.nan
                    Vz[i,j] = np.nan
                else:  #如果不在海底，则根据sigmalevelfinder找到的层数插值
                    level = sigmalevelfinder2(i,j,sd)
                    uz[i,j] = u[i,j,level]
                    vz[i,j] = u[i,j,level]
                    saltz[i,j] = salt[i,j,level]
                    tempz[i,j] = temp[i,j,level]
                    Vz[i,j] = math.sqrt((uz[i,j])*(uz[i,j])+(vz[i,j])*(vz[i,j]))
                    Vz[Vz==0] = np.nan

    #按规则命名文件并存储
    def zu_field_namegen(x,t):
        return 'u_'+str(x)+'meters'+str(t)+'hours_flowfield.txt'
    def zv_field_namegen(x,t):
        return 'v_'+str(x)+'meters'+str(t)+'hours_flowfield.txt'
    def zVelocity_field_namegen(x,t):
        return 'Velocity_'+str(x)+'meters'+str(t)+'hours_flowfield.txt'
    def zsalt_field_namegen(x,t):
        return 'salt_'+str(x)+'meters'+str(t)+'hours_field.txt'
    def ztemp_field_namegen(x,t):
        return 'temp_'+str(x)+'meters'+str(t)+'hours_field.txt'
    np.savetxt(zu_field_namegen(sd,time),uz)    #勿写成np.save，那样就存成了numpy独占的npy格式数据了
    np.savetxt(zv_field_namegen(sd,time),vz)
    np.savetxt(zVelocity_field_namegen(sd,time),Vz)
    np.savetxt(zsalt_field_namegen(sd,time),saltz)
    np.savetxt(ztemp_field_namegen(sd,time),tempz)

    move(zu_field_namegen(sd,time),plotdir)
    move(zv_field_namegen(sd,time),plotdir)
    move(zVelocity_field_namegen(sd,time),plotdir)
    move(zsalt_field_namegen(sd,time),plotdir)
    move(ztemp_field_namegen(sd,time),plotdir)

    if interp_method == 'linear':
        print('\n','方法：线性插值'+'\n'+'相关文件已存储至'+plotdir)
    if interp_method == 'nearest':
        print('\n','方法：最近邻插值'+'\n'+'相关文件已存储至'+plotdir)
    return uz,vz,Vz,saltz,tempz
#--------------------分割线-------------------------------------------------------


if __name__ == "__main__":
    try:
        pass
        #for sd in [10,50,100,200,300,400,500]:
        #for sd in [20]:
            #depth_slice_linear(hisname,720,sd)    #文件，时间（小时），深度
            #depth_slice_linear(hisname,4344,sd)
            #print(sd,'done')
        for i in np.linspace(1900,2620,31):
            i = int(i)
            depth_slice(hisname,i,10)
        
        '''
        for lon in np.arange(112.5,115,0.36):
            for lat in np.arange(18.5,20.5,0.36):
                point_uv2txt(hisname,lon,lat)
        '''
    finally:
        print('\n','done!')
        pass



