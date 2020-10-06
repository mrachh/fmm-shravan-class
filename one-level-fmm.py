from numpy import *
from pylab import *
from collections import Counter
import time

def one_level_fmm(n=1000,p=40,nbox1d=5):

# Generate sources in the unit square
    src = rand(2,n)

# Tree sort sources into a uniform grid of boxes
    srcsort = zeros((2,n))
    nboxes = nbox1d*nbox1d
    centers = zeros(nboxes,dtype='complex')
    h = 1.0/nbox1d

    boxind = zeros(n,dtype=int)
    for i in range(n):
        iind = int(src[0,i]/h)
        jind = int(src[1,i]/h)
        boxind[i] = iind*nbox1d + jind
    
    boxindsort = sort(boxind)
    isrcsort = argsort(boxind)
    isrcsortinv = argsort(isrcsort)
    srcsort = src[:,isrcsort]


    [vc,nc] = unique(boxindsort,return_counts=True) 
    nsrc = zeros(nboxes,dtype='int')
    for i in range(size(vc)):
        nsrc[vc[i]] = nc[i]
    isrcstart = zeros(nboxes+1,dtype='int')
    isrcstart[1:n] = cumsum(nsrc)


    mpole = zeros((p+1,nboxes),dtype='complex')

    zsrc = src[0,:] + 1j*src[1,:]
    zsrcsort = srcsort[0,:] + 1j*srcsort[1,:]

    t0 = time.time()
    # Form mulitpole expansion
    for ibox in range(nboxes):
        i1 = int(ibox/nbox1d)
        i2 = ibox%nbox1d
        centers[ibox] = i1*h + h/2 + 1j*(i2*h+h/2)
    
        for isrc in range(isrcstart[ibox],isrcstart[ibox+1]):
            mpole[:,ibox] = mpole[:,ibox] + (zsrcsort[isrc]-centers[ibox])**(range(p+1))
    t1 = time.time()
    tformmp = t1-t0

    # evaluate multipole expansions of boxes that are well separated

    potsort = zeros(n,dtype='complex')
    t0 = time.time()
    for ibox in range(nboxes):
        i1 = int(ibox/nbox1d)
        i2 = ibox%nbox1d
        for jbox in range(nboxes):
            j1 = int(jbox/nbox1d)
            j2 = jbox%nbox1d


            if(abs(i1-j1)>1 or abs(i2-j2)>1):
                for isrc in range(isrcstart[ibox],isrcstart[ibox+1]):
                    ytmp = 1.0/(zsrcsort[isrc]-centers[jbox])**range(1,p+2)
                    potsort[isrc] = potsort[isrc] + dot(ytmp,mpole[:,jbox])

    t1 = time.time()
    tmpeval = t1-t0
    # compute directly the contribution from nearby boxes
    t0 = time.time()
    for ibox in range(nboxes):
        i1 = int(ibox/nbox1d)
        i2 = ibox%nbox1d
        for jbox in range(nboxes):
            j1 = int(jbox/nbox1d)
            j2 = jbox%nbox1d
            if(abs(i1-j1)<=1 and abs(i2-j2)<=1):
                for isrc in range(isrcstart[ibox],isrcstart[ibox+1]):
                    for jsrc in range(isrcstart[jbox],isrcstart[jbox+1]):
                        if(isrc != jsrc):
                            potsort[isrc] = potsort[isrc] + 1.0/(zsrcsort[isrc]-zsrcsort[jsrc])

    t1 = time.time()
    tnear = t1-t0

    ttot = tnear + tformmp + tmpeval
    tfrac = zeros(3)
    tfrac[0] = tformmp/ttot
    tfrac[1] = tmpeval/ttot
    tfrac[2] = tnear/ttot


    pot = potsort[isrcsortinv]       
    # test accuracy of the computed potential at the first 10
    # targets
    ntest = min(10,n)
    potex = zeros(ntest,dtype='complex')
    for i in range(ntest):
        for j in range(n):
            if(i != j):
                potex[i] = potex[i] + 1.0/(zsrc[i]-zsrc[j])
    erra = norm(potex-pot[0:ntest])/norm(potex)
    print('relative error=',erra)

    return erra,ttot,tfrac

# Main program begins here
ntest = 4
ttot = zeros(ntest)
erra = zeros(ntest)
tfrac = zeros((3,ntest))
print("=====")
p = 20
for i in range(ntest):
    n = int(1000*2**(i))

    nbox1d = int((n)**(0.25))
    print("n=",n,"  nbox1d=",nbox1d)t 
    
    erra[i],ttot[i],tfrac[:,i] = one_level_fmm(n=n,p=p,nbox1d=nbox1d)
    print("erra=",erra[i])
    print("ttot=",ttot[i])
    print("tfrac=",tfrac[:,i])
    print(" ")
    print(" ")

