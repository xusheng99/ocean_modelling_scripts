import os
# 1. 指定ncdir为存有若干hycom数据的文件夹
# 2. 运行该程序, 在上述文件夹中会多出一个'a.txt'
# 3. 打开a.txt, 复制里面的内容

ncdir = 'D:\\harmonic_ana'
os.chdir(ncdir)
aaa = os.listdir(ncdir)

with open('a.txt','w',encoding='utf-8') as f:
    for eachfile in aaa:
        f.write(str(eachfile)+' ')
f.close()
