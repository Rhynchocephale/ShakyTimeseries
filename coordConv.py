import numpy as np
import os
import ellipsoid
import shutil

pi = np.pi
a = float(6378137)
f = 1/298.257222101
phi0 = float(0)
lamb0 = float(173)*pi/180
N0 = float(10000000)
E0 = float(1600000)
k0 = float(0.9996)
b = a*(1-f)
e2 = 1-(b/a)**2
e4 = e2**2
e6 = e2**3

A0 = 1 - (e2/4) - (3*e4/64) - (5*e6/256)
A2 = (3.0/8) * (e2+e4/4+15*e6/128)
A4 = (15.0/256) * (e4 + 3*e6/4)
A6 = 35*e6/3072

n = f/(2-f)
n2 = n*n
n3 = n2*n
n4 = n2*n2
G = a*(1-n)*(1-n2)*(1+9*n2/4+225*n4/64)

B2 = 3*n/2 - 27*n3/32
B4 = 21*n2/16 - 55*n4/32
B6 = 151*n3/96
B8 = 1097*n4/512

def meridianDist(phi):
    return  a*(A0*phi-A2*np.sin(2*phi)+A4*np.sin(4*phi)-A6*np.sin(6*phi))

m0 = meridianDist(phi0)

def footpoint(mPrime):

    sigma = mPrime/G
    phiPrime = sigma + B2*np.sin(2*sigma) + B4*np.sin(4*sigma) + B6*np.sin(6*sigma) + B8*np.sin(8*sigma)
    
    return phiPrime

def trans2geo(E,N):

    cn1  =  (N - N0)/k0 + m0 
    phiPrime = footpoint(cn1)
    sinPhiPrime = np.sin(phiPrime)
    cosPhiPrime = np.cos(phiPrime)

    eslt = (1-e2*sinPhiPrime*sinPhiPrime)
    nu = a/np.sqrt(eslt)
    rho = nu * (1-e2) / eslt
    psi = nu/rho

    EPrime = E-E0
    x = EPrime/(k0*nu)
    x2 = x*x

    t = sinPhiPrime/cosPhiPrime
    t2 = t*t
    t4 = t2*t2

    trm1 = 0.5
    trm2 = ((-4*psi+9*(1-t2))*psi+12*t2)/24
    trm3 = ((((8*(11-24*t2)*psi - 12*(21-71*t2))*psi + 15*((15*t2-98)*t2+15))*psi + 180*((-3*t2+5)*t2))*psi + 360*t4)/720
    trm4 = (((1575*t2+4095)*t2+3633)*t2+1385)/40320

    phi = phiPrime+(t*x*EPrime/(k0*rho))*(((trm4*x2-trm3)*x2+trm2)*x2-trm1)

    trm1 = 1
    trm2 = (psi+2*t2)/6
    trm3 = (((-4*(1-6*t2)*psi + (9-68*t2))*psi + 72*t2)*psi + 24*t4)/120
    trm4 = (((720*t2+1320)*t2+662)*t2+61)/5040

    lamb = lamb0 - (x/cosPhiPrime)*(((trm4*x2-trm3)*x2+trm2)*x2-trm1)

    return lamb,phi
    

def geo2trans(lamb,phi):

    dlon  =  lamb - lamb0
    while dlon > pi:
        dlon -= 2*pi 
    while dlon < -pi:
        dlon += 2*pi

    m = meridianDist(phi)

    sinPhi = np.sin(phi)

    eslt = (1-e2*sinPhi*sinPhi)
    nu = a/np.sqrt(eslt)
    rho = nu * (1-e2) / eslt
    psi = nu/rho

    cosPhi = np.cos(phi)
    w = dlon

    wc = cosPhi*w
    wc2 = wc*wc

    t = sinPhi/cosPhi
    t2 = t*t
    t4 = t2*t2
    t6 = t2*t4

    trm1 = (psi-t2)/6
    trm2 = (((4*(1-6*t2)*psi + (1+8*t2))*psi - 2*t2)*psi+t4)/120
    trm3 = (61 - 479*t2 + 179*t4 - t6)/5040

    E = (k0*nu*dlon*cosPhi)*(((trm3*wc2+trm2)*wc2+trm1)*wc2+1) 
    E = E+E0

    trm1 = 0.5
    trm2 = ((4*psi+1)*psi-t2)/24
    trm3 = ((((8*(11-24*t2)*psi - 28*(1-6*t2))*psi + (1-32*t2))*psi - 2*t2)*psi + t4)/720
    trm4 = (1385-3111*t2+543*t4-t6)/40320

    N = (nu*t)*((((trm4*wc2+trm3)*wc2+trm2)*wc2+trm1)*wc2)
    N = (N+m-m0)*k0+N0

    return E,N
   
def xyz2trans(xyz):
    lamb,phi,h = ellipsoid.grs80.geodetic(xyz) #in xyz, out degrees
    E,N = geo2trans(np.radians(lamb),np.radians(phi)) #in rad out enu
    
    return [E,N,h]

def trans2xyz(enu):
    lamb,phi = trans2geo(enu[0],enu[1])  #in m out rad
    xyz = ellipsoid.grs80.xyz(np.degrees(lamb),np.degrees(phi),enu[2])  #in degrees out m

    return xyz
   
def all_xyz2trans(folder):
    
    if not os.path.exists(folder+'/mercator/'): #creates the /transverseMercator directory if it does not already exists
            os.makedirs(folder+'/mercator/')

    #BROWSE ALL FILES IN FOLDER TO GET THEIR NAME
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    staCount = 1
    for sta in fileList:
        print('xyz->enu: '+sta.split('_')[0]+' ('+"{:.0f}".format(staCount)+'/'+"{:.0f}".format(len(fileList))+')')

        f = open(folder +'/'+ sta,'r')
        g = open(folder+'/mercator/'+ sta.replace('xyz','enu'),'w')
        
        f.readline()
        g.write('name epoch x y z\n') #header
        for line in f:
            line = line.split()
            xyz = map(float,line[2:5]) #out xyz
            ENU = xyz2trans(xyz)
            
            g.write(line[0]+' '+line[1]+' '+"{:.6f}".format(ENU[0])+' '+"{:.6f}".format(ENU[1])+' '+"{:.6f}".format(ENU[2])+'\n')
        
        f.close()
        g.close()
    
        staCount += 1
    
    return 0
    
def all_trans2xyz(folder):
    
    if not os.path.exists(folder+'/xyz/'): #creates the /transverseMercator directory if it does not already exists
            os.makedirs(folder+'/xyz/')

    #BROWSE ALL FILES IN FOLDER TO GET THEIR NAME
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    staCount = 1
    for sta in fileList:
        print('enu->xyz: '+sta.split('_')[0]+' ('+"{:.0f}".format(staCount)+'/'+"{:.0f}".format(len(fileList))+')')

        f = open(folder +'/'+ sta,'r')
        g = open(folder+'/xyz/'+ sta.replace('xyz','enu'),'w')

        f.readline()
        g.write('name epoch               x             y           z            flag stations_for_smoothing\n') #header
        for line in f:
            line = line.split()
            enu = map(float,line[2:5]) #get transverse mercator out enu
            xyz = trans2xyz(enu)
            g.write(line[0]+' '+line[1]+' '+"{:.6f}".format(xyz[0])+' '+"{:.6f}".format(xyz[1])+' '+"{:.6f}".format(xyz[2])+'\n')
        
        f.close()
        g.close()
    
        staCount += 1
    
    
#simplest way to go from enu to xyz
def all_xyz2enu(folder):
    #BROWSE ALL FILES IN FOLDER TO GET THEIR NAME
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    if os.path.exists(folder+'/enu/'): #clears the directory if it already exists
        shutil.rmtree(folder+'/enu/')
    os.makedirs(folder+'/enu/')
    os.makedirs(folder+'/enu/param/')
    
    i=1
    for fileName in fileList:
        print('xyz->enu: '+fileName.split('/')[-1].split('_')[0]+' ('+str(i)+'/'+str(len(fileList))+')')
        f=open(folder+'/'+fileName,'r')
        g=open(folder+'/enu/'+fileName.split('_xyz')[0]+'_enu.dat','w')
        f.readline() #header
        
        firstLine = f.readline().split()
        xyz0 = np.array(map(float,firstLine[2:5])) #gets the first xyz position
        lamb,phi,h = ellipsoid.grs80.geodetic(xyz0)
        enu_axes_matrix = ellipsoid.grs80.enu_axes(lamb,phi)
        h = open(folder+'/enu/param/'+fileName.split('_')[0]+'.param','w')
        h.write('xyz0 '+"{:.6f} {:.6f} {:.6f}".format(xyz0[0],xyz0[1],xyz0[2])+'\n')
        h.write('enu_axes '+'{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}'.format(enu_axes_matrix[0][0],enu_axes_matrix[0][1],enu_axes_matrix[0][2],enu_axes_matrix[1][0],enu_axes_matrix[1][1],enu_axes_matrix[1][2],enu_axes_matrix[2][0],enu_axes_matrix[2][1],enu_axes_matrix[2][2]))
        h.close()
        g.write('name	epoch	e   n  u\n')
        g.write(firstLine[0]+'\t'+firstLine[1]+'\t0.000000 0.000000 0.000000\n')
        for line in f:
            splittedLine = line.split()
            xyz=np.expand_dims(np.array(map(float,splittedLine[2:5]))-xyz0,1)
            #print(xyz.shape)
            enu=1000*enu_axes_matrix.dot(xyz) #-xyz0)
            #print(enu[0][0])
            g.write(splittedLine[0]+'\t'+splittedLine[1]+'\t'+"{:.6f} {:.6f} {:.6f}".format(enu[0][0],enu[1][0],enu[2][0])+'\n')
        
        g.close()
        f.close()
        i+=1

def all_enu2xyz(folder):
    #BROWSE ALL FILES IN FOLDER TO GET THEIR NAME
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f.split('_')[0])
    fileList.sort()
    
    if os.path.exists(folder+'/smoothed/xyz'): #clears the directory if it already exists
        shutil.rmtree(folder+'/smoothed/xyz')
    os.makedirs(folder+'/smoothed/xyz')

    u=1
    for sta in fileList:
        print('enu->xyz: '+sta+' ('+str(u)+'/'+str(len(fileList))+')')
        
        f2 = open(folder+'/enu/param/'+sta+'.param')
        xyz0 = np.array(map(float,f2.readline().split()[1:]))
        enu_axes_matrix=np.reshape(map(float,f2.readline().split()[1:]),(3,3)).transpose() #creates the matrix with the axes from line read
        
        f2.close()
        
        f = open(folder+'/smoothed/'+sta+'_igs08_enu_smoothed.dat','r')
        g = open(folder+'/smoothed/xyz/'+sta+'_igs08_xyz_smoothed.dat','w')
        g.write('name	epoch	x	y	z\n')
        f.readline() #header
        for line in f:
            splittedLine = line.split()
            enu = np.array(map(float,splittedLine[2:5]))
            if any(enu):
                xyz = (xyz0 + enu_axes_matrix.dot(enu/1000))
                g.write(splittedLine[0]+'\t'+splittedLine[1]+'\t'+"{:.6f} {:.6f} {:.6f}".format(xyz[0],xyz[1],xyz[2])+'\n')
        f.close()
        g.close()
        u+=1
        
'''
#print('geo2trans')
E,N = geo2trans(np.radians(172),np.radians(-36))
lamb,phi = trans2geo(E,N)
print(np.degrees(lamb),np.degrees(phi))


lamb,phi = trans2geo(1750000,5500000)
#ellipsoid.grs80.xyz(lamb,phi,300)
x,y = geo2trans(lamb,phi)
print(x,y)

folder = '/home/cdrouadaine/timeseries/timeseries/test'
all_xyz2trans(folder)
all_trans2xyz(folder+'/mercator/')
'''

#all_xyz2enu('/home/cdrouadaine/timeseries/timeseries')
#print(ellipsoid.grs80.geodetic([-4052051.9525,4212836.0919,-2545105.6746]))
