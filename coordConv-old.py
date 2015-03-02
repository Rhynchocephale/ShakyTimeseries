import numpy as np
import os
import ellipsoid

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

#meridian distance
A0 = 1-(e2/4)-(3*e2**2)/64-(5*e2**3)/256
A2 = (3/8)*(e2+(e2**2)/4+(15*e2**3)/128)
A4 = (15/256)*(e2**2+(3*e2**3)/4)
A6 = (35*e2**3)/3072

#Radius of curvature. In radians
def radius(phi):
    under = 1-e2*(np.sin(phi))**2
    rho = a*(1-e2)/(under**(3/2))
    nu = a/np.sqrt(under)
    psi = nu/rho
    #r2 = rho*nu*k0**2

    #print(rho,psi,nu)
    return rho,psi,nu

#in radians
def meridianDist(phi):
    return a*(A0*phi-A2*np.sin(2*phi)+A4*np.sin(4*phi)-A6*np.sin(6*phi))

#foot-point latitude
def footpoint(N):
    n = (a-b)/(a+b)
    G = a*(1-n)*(1-n**2)*(1+(9*n**2)/4+(225*n**4)/64)*(pi/180)
    m0 = meridianDist(phi0)
    NPrime = N-N0
    mPrime = m0 + NPrime/k0
    sigma = mPrime*np.pi/(180*G)  ###RAD OR DEG ???
    sigmaRad = sigma*pi/180
    B2 = (3*n)/2-(27*n**3)/32
    B4 = (21*n**2)/16-(55*n**4)/32
    B6 = (151*n**3)/96
    B8 = (1097*n**4)/512
    phiPrime = sigma + B2*np.sin(2*sigmaRad)+B4*np.sin(4*sigmaRad)+B6*np.sin(6*sigmaRad)+B8*np.sin(8*sigmaRad)   ###RAD OR DEG ???
    
    return phiPrime

#in degrees, out metres
def geo2trans(lamb,phi):
    lambRad = lamb*pi/180
    phiRad = phi*pi/180
    m = meridianDist(phi)
    m0 = meridianDist(phi0)
    rho, psi, nu = radius(phi)
    t = np.tan(phiRad)
    omega = lamb-lamb0
    cosphi = np.cos(phiRad)
    nusinphi = nu*np.sin(phiRad)
    
    #longitude
    T1 = ((omega**2)/6)*(cosphi**2)*(psi-t**2)
    T2 = ((omega**4)/120)*(cosphi**4)*(4*(1-6*t**2)*psi**3+(1+8*t**2)*psi**2-2*psi*t**2+t**4)
    T3 = ((omega**6)/5040)*(cosphi**6)*(61-479*t**2+179*t**4-t**6)
    EPrime = k0*nu*omega*cosphi*(1+T1+T2+T3)
    E = EPrime + E0

    #latitude
    T1 = ((omega**2)/2)*nusinphi*cosphi
    T2 = ((omega**4)/24)*nusinphi*(cosphi**3)*(4*psi**2+psi-t**2)
    T3 = ((omega**6)/720)*nusinphi*(cosphi**5)*(8*(11-24*t**2)*psi**4-28*(1-6*t**2)*psi**3+(1-32*t**2)*psi**2-2*psi*t**2+t**4)
    T4 = ((omega**8)/40320)*nusinphi*(cosphi**7)*(1385-3111*t**2+543*t**4-t**6)
    NPrime = k0*(m-m0+T1+T2+T3+T4)
    N = NPrime + N0

    return E,N

#in metres, out degrees    
def trans2geo(E,N):
    phiPrime = footpoint(N)
    rho, psi, nu = radius(phiPrime)
    t = np.tan(phiPrime)
    EPrime = E-E0
    x = EPrime/(k0*rho*nu)
    y = x*(EPrime/(k0*rho))
    
    #northing
    B = t*EPrime/(k0*rho)
    T1 = B*x/2
    T2 = B*((x**3)/24)*(-4*psi**2+9*(1-t**2)*psi+12*t**2)
    T3 = B*((x**5)/720)*(8*(11-24*t**2)*psi**4-12*(21-71*t**2)*psi**3+15*(15-98*t**2+15*t**4)*psi**2+180*(5*t**2-3*t**4)*psi+360*t**4)
    T4 = B*((x**7)/40320)*(1385+3633*t**2+4095*t**4+1575*t**6)
    phi = phiPrime-T1+T2-T3+T4
        
    #easting 
    T1 = x/np.cos(phiPrime)
    T2 = T1*((x**2)/6)*(psi+2*t**2)
    T3 = T1*((x**4)/120)*(-4*(1-6*t**2)*psi**3+(9-68*t**2)*psi**2+72*psi*t**2+24*t**4)
    T4 = T1*((x**6)/5040)*(61+662*t**2+1320*t**4+720*t**6)
    lamb = lamb0+T1-T2+T3-T4

    return lamb,phi

def all_xyz2trans(folder):
    
    if not os.path.exists(folder+'transverseMercator/'): #creates the /transverseMercator directory if it does not already exists
            os.makedirs(folder+'transverseMercator/')

    #BROWSE ALL FILES IN FOLDER TO GET THEIR NAME
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f)

    staCount = 1
    for sta in fileList:
        print('converting station: '+sta.split('_')[0]+' ('+"{:.0f}".format(staCount)+'/'+"{:.0f}".format(len(fileList))+')')

        f = open(folder + sta,'r')
        g = open(folder+'transverseMercator/'+ sta.split('.')[0]+'_merc.dat','w')
        
        f.readline()
        g.write('name epoch               e             n           u            flag stations_for_smoothing\n') #header
        for line in f:
            line = line.split()
            xyz = map(float,line[2:5]) #out xyz
            lamb,phi,h = ellipsoid.grs80.geodetic(xyz) #in xyz, out degrees
            E,N = geo2trans(lamb,phi) #in degrees out enu
            
            g.write(line[0]+' '+line[1]+' '+"{:.4f}".format(E)+' '+"{:.4f}".format(N)+' '+"{:.4f}".format(h)+' '+line[5]+'\n')
        
        f.close()
        g.close()
    
        staCount += 1
    
    return 0
    
def all_trans2xyz(folder):
    
    if not os.path.exists(folder+'xyz/'): #creates the /transverseMercator directory if it does not already exists
            os.makedirs(folder+'xyz/')

    #BROWSE ALL FILES IN FOLDER TO GET THEIR NAME
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f)

    staCount = 1
    for sta in fileList:
        print('converting station: '+sta.split('_')[0]+' ('+"{:.0f}".format(staCount)+'/'+"{:.0f}".format(len(fileList))+')')

        f = open(folder + sta,'r')
        g = open(folder+'xyz/'+ sta.split('.')[0]+'_xyz.dat','w')

        f.readline()
        g.write('name epoch               x             y           z            flag stations_for_smoothing\n') #header
        for line in f:
            line = line.split()
            xyz = map(float,line[2:5]) #get transverse mercator out enu
            lamb,phi = trans2geo(xyz[0],xyz[1]) #trans2geo  in enu out degrees
            xyz2 = ellipsoid.grs80.xyz(lamb,phi,xyz[2]) #geo2xyz   in degrees out xyz

            #g.write(line[0]+' '+line[1]+' '+"{:.4f}".format(xyz2[0])+' '+"{:.4f}".format(xyz2[1])+' '+"{:.4f}".format(xyz2[2])+' '+line[5]+' '+line[6]+'\n')
            g.write(line[0]+' '+line[1]+' '+"{:.4f}".format(xyz2[0])+' '+"{:.4f}".format(xyz2[1])+' '+"{:.4f}".format(xyz2[2])+' '+line[5]+'\n')
        
        f.close()
        g.close()
    
        staCount += 1
    
    return 0

#all_xyz2trans('/home/cdrouadaine/timeseries/TEST/')
#all_trans2xyz('/home/cdrouadaine/timeseries/TEST/transverseMercator/')
E,N = geo2trans(-41,174)
print(E,N)
lamb,phi = trans2geo(E,N)
print(lamb,phi)
