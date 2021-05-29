import os
import shutil

def move(sourcefile,targetdir):    #移动sourcefile至targetdir,如果已存在同名文件则覆盖
    pwd = os.getcwd()
    os.chdir(targetdir)
    filelist = os.listdir(targetdir)
    if sourcefile in filelist:
        os.remove(sourcefile)
    os.chdir(pwd)
    shutil.move(sourcefile,targetdir)
    print(str(sourcefile),'has been moved to',str(targetdir))

