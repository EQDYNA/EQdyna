import numpy as np
from numpy.fft import fft2, ifft2
import matplotlib.pyplot as plt
from user_defined_params import par

def generateFractalSurface(lx, hurst_exponent, seed=None):
    np.random.seed(seed)
    
    lh = lx//2
    pony = hurst_exponent + 1.
    hrttwo = 0.5*np.sqrt(2.0)
    kcut = lx/4
    fcsq = kcut*kcut
    npole = 4
    
    w = np.zeros((lx,lx),dtype=np.float64)
    transf = np.zeros((lx,lx),dtype=np.complex128)
    
    transf = hrttwo*(np.random.normal(size=(lx,lx)) + 
                     1j*np.random.normal(size=(lx,lx)))
    for i in range(1,lx+1):
        for j in range(1,lx+1):
            # calculate square of wavenumber
            ik = i-1 
            jk = j-1
            if ik>lh:
                ik = ik - lx
            if jk>lh:
                jk = jk - lx
            ksq = ik*ik+jk*jk
            fsq = ksq
            if ksq==0:
                transf[i-1,j-1] = 0.0
                pspec = 0.0
            #elif ksq == 1:
            #    if jk == 0:
            #        transf[i-1,j-1] = -2.0
            #        pspec = 1.0/fsq**pony
            #    else:
            #        transf[i-1,j-1] = 0.0
            #        pspec = 1.0/fsq**pony
            #elif ksq == 2:
            #    transf[i-1,j-1] = 0.0
            #    pspec = 1.0/fsq**pony
            else:
                pspec = 1.0/fsq**pony
                       
            filt  = 1.0/(1.0+(fsq/fcsq)**npole)
            transf[i-1,j-1] = np.sqrt(pspec*filt)*transf[i-1,j-1]
    
    transf = ifft2(transf)
    coef = 1.0/4.0
    w = coef * np.real(transf)

    return w

hurst_exponent = 0.7
seedId = 48
alpha0 = 0.005


w = generateFractalSurface(par.nfx, hurst_exponent, seedId)
mean_w = np.mean(w)
std_w  = np.std(w)
print(mean_w, std_w)
fractalSurface     = w*alpha0/(std_w/par.nfx/par.dx)
dxFractalSurface, dyFractalSurface = np.gradient(fractalSurface/par.dx, axis=(0,1))

subSurface = fractalSurface[:par.nfx, :par.nfz]
subDx      = dxFractalSurface[:par.nfx, :par.nfz]
subDy      = dyFractalSurface[:par.nfx, :par.nfz]

result = np.zeros((par.nfx*par.nfz+2,3))

result[0,0] = par.nfx
result[0,1] = par.nfz
result[1,0] = par.dx
result[1,1] = par.fxmin
result[1,2] = par.fzmin
ntmp = 2 # skipping the first two lines
for i in range(par.nfx):
    for j in range(par.nfz):
        result[ntmp,0] = subSurface[i,j]
        result[ntmp,1] = subDx[i,j]
        result[ntmp,2] = subDy[i,j]
        ntmp = ntmp + 1

np.savetxt('bFault_Rough_Geometry.txt', result, fmt='%f', delimiter='\t')

print('Peak roughness of the fractal surface is ', np.max(fractalSurface),' m')

fig = plt.figure(figsize=(24,6), dpi=300, facecolor='w', edgecolor='w')
fig.add_subplot(131)
plt.contourf(fractalSurface)
ax = plt.gca()
ax.set_aspect('equal')
plt.colorbar()
plt.title(f'2D Fractal Surface (H={hurst_exponent})')

fig.add_subplot(132)
plt.contourf(dxFractalSurface)
ax = plt.gca()
ax.set_aspect('equal')
plt.colorbar()
plt.title('Derivative along x')

fig.add_subplot(133)
plt.contourf(dyFractalSurface)
ax = plt.gca()
ax.set_aspect('equal')
plt.colorbar()
plt.title('Derivative along y')

plt.savefig('cFractalSurface.png', dpi = 600)
