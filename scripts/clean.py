import os

np = 10
fnameList = [['gm',''], ['rsa','.mat'], ['frt.txt',''], 
        ['surface_coor.txt','']]

for f in fnameList:
    for i in range(np):
        fname = f[0]+str(i)+f[1]
        if os.path.exists(fname):
            print('Removing ', fname, ' ...')
            os.remove(fname) 
