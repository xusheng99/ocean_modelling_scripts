from netCDF4 import Dataset
import os
import numpy as np 
import math
import shutil
#----------------函数中用到的参数----------------
    #经纬度范围
lon1 = 105
lon2 = 124
lat1 = 5
lat2 = 24

    #sigma分层参数
theta_s,theta_b,hc,N = [1.0,0.4,1.0,25]
sc = [theta_s,theta_b,hc,N]

    #指定的z分层每层的深度
zarray = [3,5,10,15,20,40,50,100,150,200,300,500,1000,1500,2000,3000]

    #工作目录和绘图目录
workdir = 'C:\\Users\\XUSHENG\\Desktop\\read_roms_result\\'    #将某些较小的nc文件直接放到这里
plotdir = 'C:\\Users\\XUSHENG\\Desktop\\作图用文件\\'    #读取的数据移动到这里供matlab使用
os.chdir(workdir)
    
    #netcdf文件的位置
ininame = 'roms_ini.nc'
avgname = 'd:\\scs_model\\scs_avg.nc'    #因为这个文件一般比较大，不建议搬运到本文件夹后再处理
hisname = 'd:\\scs_model\\scs_his.nc'    ## 不要打错字，好好检查一下路径
gridname = 'roms_grid.nc'
bryname = 'roms_bry.nc'
#fcaininame = 'H:\\NHS\\nh_ini.nc'
ini = Dataset(ininame)
grid = Dataset(gridname)
bry = Dataset(bryname)
avg = Dataset(avgname)
his = Dataset(hisname)
#fcaini = Dataset(fcaininame)

    #掩码和深度
maskr = grid.variables['mask_rho'][:]
masku = grid.variables['mask_u'][:]
maskv = grid.variables['mask_v'][:]
depth = grid.variables['h'][:]
depth = depth*maskr 
lon_length = maskr.shape[0]    #从mask读取网格的大小！
lat_length = maskr.shape[1]

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

u = his.variables['u'][720,:,:,:].data
v = his.variables['v'][720,:,:,:].data
temp = his.variables['temp'][720,:,:,:].data
salt = his.variables['salt'][720,:,:,:].data

#将u网格和v网格上的值插到rho网格上
uap = np.empty((depth.shape[1],depth.shape[0],N))
vap = np.zeros((depth.shape[1],depth.shape[0],N))
Va = np.ones((depth.shape[1],depth.shape[0],N))
for i in range(depth.shape[1]):
    for j in range(depth.shape[0]):
        for k in range(N):
            if i == 0:#自左向右第一个rho点
                uap[i,j,k] = 1.5*u[k,j,i]-0.5*u[k,j,i+1]
            if i == (depth.shape[1]-1):#自左向右最后一个rho点
                uap[i,j,k] = 1.5*u[k,j,i-1]-0.5*u[k,j,i-2]
            else:
                uap[i,j,k] = 0.5*(u[k,j,i-1]+u[k,j,i])
            if j == 0:
                vap[i,j,k] = 1.5*v[k,j,i]+0.5*v[k,j+1,i]
            if j == (depth.shape[0]-1):
                vap[i,j,k] = 1.5*v[k,j-1,i]-0.5*v[k,j-2,i]
            else:
                vap[i,j,k] = 0.5*(v[k,j-1,i]+v[k,j,i])
            uap[i,j,k] = uap[i,j,k]    #这里就先不写*maskr，看看怎么样
            vap[i,j,k] = vap[i,j,k]
            Va[i,j,k] = math.sqrt(uap[i,j,k]*uap[i,j,k]+vap[i,j,k]*vap[i,j,k]) 
            
temp = temp.transpose(2,1,0)
salt = salt.transpose(2,1,0)

def cleanuv(x):
    x[abs(x>10)] = np.nan

cleanuv(uap)
cleanuv(vap)
temp[abs(temp>50)] = np.nan
salt[abs(temp)>50] = np.nan



with open('salt.txt','w') as f:
    for i in range(uap.shape[0]):
        for j in range(uap.shape[1]):
            for k in range(uap.shape[2]):
                f.write(str(i)+' '+str(j)+' '+str(z[i,j,k])+' '+str(salt[i,j,k])+'\n')
print('done')
f.close()


         
