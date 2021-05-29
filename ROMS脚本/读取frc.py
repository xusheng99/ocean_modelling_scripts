import os

ncdir = 'D:\\application\\PLANT\\DATA\\frc'
os.chdir(ncdir)
aaa = os.listdir(ncdir)

with open('a.txt','w',encoding='utf-8') as f:
    for eachfile in aaa:
        f.write('./DATA/frc/'+str(eachfile)+' '+' \ '+'\n')
f.close()
