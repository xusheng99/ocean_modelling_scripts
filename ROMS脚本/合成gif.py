import imageio
import os
import time

os.chdir('C:\\Users\\XUSHENG\\Desktop\\particle_tracking_results\\plant2\\gif')
filelist = os.listdir()
a = []
for eachfile in filelist:
    if '.gif' in eachfile:
        pass
    else:
        a.append(eachfile)
filelist = a
print(filelist)

#name = time.asctime( time.localtime(time.time()) )
#name = name + '.gif'
name = 'ddd.gif'

img_paths = filelist
gif_images = []
for path in img_paths:
    gif_images.append(imageio.imread(path))
imageio.mimsave(name,gif_images,fps=10)
