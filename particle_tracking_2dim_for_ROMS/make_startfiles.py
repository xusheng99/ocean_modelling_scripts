import os
import sys

N = int(sys.argv[1])
initime = sys.argv[2]
workdir = sys.argv[3]

with open('1core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 1core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py 0 '+str(N//8)+' '+initime+' '+workdir)
f.close()

with open('2core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 2core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py '+str(N//8)+' '+str(N//8)+' '+initime+' '+workdir)
f.close()

with open('3core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 3core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py '+str(2*(N//8))+' '+str(N//8)+' '+initime+' '+workdir)
f.close()

with open('4core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 4core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py '+str(3*(N//8))+' '+str(N//8)+' '+initime+' '+workdir)
f.close()

with open('5core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 5core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py '+str(4*(N//8))+' '+str(N//8)+' '+initime+' '+workdir)
f.close()

with open('6core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 6core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py '+str(5*(N//8))+' '+str(N//8)+' '+initime+' '+workdir)
f.close()

with open('7core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 7core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py '+str(6*(N//8))+' '+str(N//8)+' '+initime+' '+workdir)
f.close()

with open('8core.bat','w',encoding='utf-8') as f:
    f.write('TITLE 8core\n')
    f.write('python ./particle_tracking_2dim_for_ROMS_b2.py '+str(7*(N//8))+' '+str(N-7*(N//8))+' '+initime+' '+workdir)
f.close()