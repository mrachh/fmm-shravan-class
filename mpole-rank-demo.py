from numpy import *
from pylab import *

def get_pts(idist=0,n=400):
#   sources and targets distributed in squares 
#   well separated from each other
    if(idist == 0):
        src = rand(2,n)
        targ = rand(2,n) 
        targ[0,:] = targ[0,:]+2
# sources and targets contained in annular discs well separated
# from each other
    if(idist == 1):
        thet = rand(n)*2*pi
        r = rand(n)
        src = zeros((2,n))
        targ = zeros((2,n))
        src[0,:] = r*cos(thet)
        src[1,:] = r*sin(thet)

        thet = rand(n)*2*pi
        r = rand(n) + 2
        targ[0,:] = r*cos(thet)
        targ[1,:] = r*sin(thet)
    return src,targ


# Get matrix of interaction between sources and targets
def get_lap_mat(src,targ):
    [k,n] = shape(src)
    [k,m] = shape(targ)
    A = zeros((m,n),dtype='complex')
    z = src[0,:] + 1j*src[1,:]
    y = targ[0,:] + 1j*targ[1,:]
    for i in range(m):
        A[i,:] = 1.0/(y[i]-z)
    return A

# Compute low rank apporximation of matrix using multipole 
# expansion
def get_mpole_approx_lap(src,targ,p):
    [k,n] = shape(src)
    [k,m] = shape(targ)
    umat = zeros((m,p+1),dtype='complex')
    vmat = zeros((p+1,n),dtype='complex')
    z = src[0,:] + 1j*src[1,:]
    y = targ[0,:] + 1j*targ[1,:]
    for i in range(p+1):
        vmat[i,:] = z**(i)
        umat[:,i] = (1.0/y**(i+1))
    amat = matmul(umat,vmat) 

    return umat,vmat,amat

# main program begins here
[src,targ] = get_pts(idist=1)
figure(1)
plot(src[0,:],src[1,:],'r.')
plot(targ[0,:],targ[1,:],'b.')
axis('equal')
show()

A = get_lap_mat(src,targ)

pmax = 100
errs = zeros(pmax+1)

for p in range(pmax+1):
    [up,vp,Ap] = get_mpole_approx_lap(src,targ,p)
    errs[p] = norm(Ap[:]-A[:])
    
figure(2)
semilogy(errs,'k.',basey=10)
show()

input("press any key to compare to svd:\n") 

[u,s,vh] = svd(A)

figure(3)
semilogy(errs,'k.',basey=10)
semilogy(s[0:pmax],'kx',basey=10,ms=5)

show()
