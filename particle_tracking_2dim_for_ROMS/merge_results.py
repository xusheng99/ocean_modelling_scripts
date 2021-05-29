import sys
import numpy as np
import os

os.chdir(sys.argv[2])

N = int(sys.argv[1])
interval = N//8

## numpy在对npy格式的数组操作时, 压根不需要将数据内容载入到内存里, 可以直接操作metainfo, 极大节省了内存
a1 = np.load(str(0)+'_'+str(interval)+'.npy')
a2 = np.load(str(interval)+'_'+str(2*interval)+'.npy')
a3 = np.load(str(2*interval)+'_'+str(3*interval)+'.npy')
a4 = np.load(str(3*interval)+'_'+str(4*interval)+'.npy')
a5 = np.load(str(4*interval)+'_'+str(5*interval)+'.npy')
a6 = np.load(str(5*interval)+'_'+str(6*interval)+'.npy')
a7 = np.load(str(6*interval)+'_'+str(7*interval)+'.npy')
a8 = np.load(str(7*interval)+'_'+str(N)+'.npy')    

a = np.concatenate((a1,a2,a3,a4,a5,a6,a7,a8),axis=0)

np.save('result.npy',a)

for i in range(a.shape[2]):
    np.savetxt(str(int(a[0,2,i]))+'.dat',a[:,:,i])