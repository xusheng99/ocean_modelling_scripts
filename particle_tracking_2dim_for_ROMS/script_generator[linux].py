import os 

startday = 180
endday = 190

os.chdir('/mnt/c/Users/XUSHENG/Documents/scripts/paticle_tracking_2dim_for_ROMS')

def w(num):
    return 'python3 ./particle_tracking_2dim_for_ROMS_v8.0.py  '+'/mnt/c/Users/XUSHENG/Desktop/particle_tracking_results/tdn/'+str(num)+' '+str(num)

with open('run.sh','w',encoding='utf-8') as f:
    f.write('# /usr/bin/bash\n')
    f.write('    \n')
    for i in range(startday,endday,3):
        f.write('if !( -e {} ); then\n'.format(str(i)))
        f.write('mkdir {}\n'.format(str(i)))
        f.write('fi\n')
        f.write('    \n')
        f.write(w(i)+'\n')
        f.write('    \n')
    f.write('echo \"done\" ')
f.close()

