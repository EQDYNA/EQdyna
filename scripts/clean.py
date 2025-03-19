import os

np = 10000
fnameList = [['gm',''], ['rsa','.mat'], ['frt.txt',''], 
        ['surface_coor.txt',''], ['pstr.txt',''], ['src_evol',''],['src','']]

for f in fnameList:
    for i in range(np):
        fname = f[0]+str(i)+f[1]
        if os.path.exists(fname):
            print('Removing ', fname, ' ...')
            os.remove(fname) 

os.system("rm -rf body*st*dp*.txt")
os.system("rm -rf faultst*dp*.txt")
os.system("rm -rf src*")
