from numpy import *
from pylab import *

n = 100
nlat = 300
src = rand(2,n)
charges = rand(n)-0.5 + 1j*(rand(n)-0.5)

# distance between source and target box
rdis = 2

# generate grid of targets 
x = linspace(rdis,rdis+1,nlat)
y = linspace(0,1,nlat)

xx,yy = meshgrid(x,y)

# complexify source and target locations
z = xx[:] + 1j*yy[:]
zsrc = src[0,:] + 1j*src[1,:]

pot_lap = zeros(shape(z),dtype='complex')
pot_gauss = zeros(shape(z),dtype='complex')

for i in range(n):
    pot_lap = pot_lap + charges[i]/(z-zsrc[i])
    pot_gauss = pot_gauss + charges[i]*exp(-abs(z-zsrc[i])**2)


vm1_lap= max(-10,amin(real(pot_lap[:])))
vm2_lap = min(10,amax(real(pot_lap[:])))

vm1_gauss= max(-10,amin(real(pot_gauss[:])))
vm2_gauss = min(10,amax(real(pot_gauss[:])))


vm1 = vm1_lap
vm2 = vm2_lap
figure(1)
clf()
plot(src[0,:],src[1,:],'b.',alpha=0.5)
plot([0,1,1,0,0],[0,0,1,1,0],'k-',linewidth=1.2)
imshow(real(pot_lap),extent=[rdis,rdis+1,1,0],vmin=vm1,vmax=vm2)
colorbar(shrink=0.9)
xlim(-0.2,rdis+1.2) 
ylim(-0.2,1.2)
show()
     

vm1 = vm1_gauss
vm2 = vm2_gauss
figure(2)
clf()
plot(src[0,:],src[1,:],'b.',alpha=0.5)
plot([0,1,1,0,0],[0,0,1,1,0],'k-',linewidth=1.2)
imshow(real(pot_gauss),extent=[rdis,rdis+1,1,0],vmin=vm1,vmax=vm2)
colorbar(shrink=0.9)
xlim(-0.2,rdis+1.2) 
ylim(-0.2,1.2)
show()
     
