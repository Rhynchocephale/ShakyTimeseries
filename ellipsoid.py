import numpy as np
import math

class ellipsoid( object ):

    convergence=1.0e-10

    @staticmethod
    def _cossin( angle ):
        angle=np.radians(angle)
        return np.cos(angle),np.sin(angle)

    @staticmethod
    def enu_axes( lon, lat ):
        '''
        Returns an array defining the east, north, and up unit vectors
        at a specified latitude and longitude
        '''
        cln,sln = ellipsoid._cossin(lon)
        clt,slt = ellipsoid._cossin(lat)
        ve=np.array([-sln,cln,0])
        vn=np.array([-cln*slt,-sln*slt,clt])
        vu=np.array([clt*cln,clt*sln,slt])
        return np.vstack((ve,vn,vu))

    def __init__( self, a, rf ):
        '''
        Initiallize an ellipsoid based on semi major axis and inverse flattening
        '''
        self.a=float(a)
        self.rf=float(rf)
        self.b=a-a/rf if rf else a
        self.a2=a*a
        self.b2=self.b*self.b
        self.a2b2=self.a2-self.b2

    def xyz( self, lon, lat, hgt=0 ):
        '''
        Calculate the geocentric X,Y,Z coordinates at a longitude 
        and latitude
        '''
        cln,sln = ellipsoid._cossin(lon)
        clt,slt = ellipsoid._cossin(lat)
        bsac=np.hypot(self.b*slt,self.a*clt)
        p = self.a2*clt/bsac + hgt*clt
        xyz=[p*cln,p*sln,self.b2*slt/bsac+hgt*slt]
        return xyz

    def geodetic( self, xyz ):
        '''
        Calculate the longitude, latitude, and height corresponding 
        to a geocentric XYZ coordinate
        '''
        x,y,z = xyz[0:3]
        ln=np.arctan2(y,x)
        p=np.hypot(x,y)
        lt=np.arctan2(self.a2*z,self.b2*p)
        for i in range(10):
            lt0=lt
            slt=np.sin(lt)
            clt=np.cos(lt)
            bsac=np.hypot(self.b*slt,self.a*clt)
            lt=np.arctan2(z+slt*self.a2b2/bsac,p)
            if abs(lt-lt0) < self.convergence:
                break
        h=p*clt+z*slt-bsac
        return np.degrees(ln),np.degrees(lt),h

grs80 = ellipsoid(6378137.0,298.257222101)

'''
xyz=grs80.xyz(-39,73,350)
print(xyz)
print(grs80.geodetic(xyz))


lon,lat,h=grs80.geodetic([-1000000,500000,-6000000])
print(lon,lat,h)
print(grs80.xyz(lon,lat,h))
'''

