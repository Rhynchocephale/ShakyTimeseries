import pylab as pl
#from transformInput import *
import datetime as dt
import time
import os
import shutil

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt.datetime(year=year, month=1, day=1)
    startOfNextYear = dt.datetime(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


def fromYearFraction(number):
    year = int(np.floor(number))
    secondsInYear = (365+calendar.isleap(year)) * 24*3600
    secondsHere = float(number-year)*secondsInYear
    date = dt.datetime(year,1,1)+dt.timedelta(seconds=secondsHere)

    return date


def drawStuff():
    
    folder = "/media/dpts/3775-F21C/timeseries/timeseries_allPositioNZfitted"
    
    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    staNames = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            staNames.append(f.split('_')[0])
    staNames.sort()
    
    if os.path.exists(folder+'/plots'): #clears the /QRFolder/output directory if it already exists
        shutil.rmtree(folder+'/plots')
    os.makedirs(folder+'/plots')
    
    for staName in staNames:
        print('Printing '+staName)
        #inp = open('folder/smoothed/xyz/'+staName+'_igs08_xyz_smoothed.dat','r')
        inp = open(folder+'/smoothed/xyz/'+staName+'_igs08_xyz_smoothed.dat','r')
        outp = open(folder+'/withOffsets/'+staName+'_igs08_xyz.dat','r')
        #outp = open('folder/withSlowSlips/'+staName+'_igs08_xyz.dat','r')
        base = open(folder+'/'+staName+'_igs08_xyz.dat','r')
        
        #skips headers
        inp.readline()
        outp.readline()
        base.readline()
        
        iX = []
        iY = []
        iZ = []
        oX = []
        oY = []
        oZ = []
        bX = []
        bY = []
        bZ = []
        dates = []
        #reads line by line, both files at once
        while True:
            iLine = inp.readline().split()
            oLine = outp.readline().split()
            bLine = base.readline().split()
            if not oLine:
                break
            
            iDate = dt.datetime.strptime(iLine[1].split('T')[0],'%Y-%m-%d').date()
            oDate = dt.datetime.strptime(oLine[1].split('T')[0],'%Y-%m-%d').date()
            bDate = dt.datetime.strptime(oLine[1].split('T')[0],'%Y-%m-%d').date()
            
            while iDate != oDate or bDate != oDate or iDate != bDate: #extremely ugly, but works
                while oDate < iDate:
                    oLine = outp.readline().split()
                    oDate = dt.datetime.strptime(oLine[1].split('T')[0],'%Y-%m-%d').date()
                    
                while iDate < oDate:
                    iLine = inp.readline().split()
                    iDate = dt.datetime.strptime(iLine[1].split('T')[0],'%Y-%m-%d').date()
                    
                while iDate < bDate:
                    iLine = inp.readline().split()
                    iDate = dt.datetime.strptime(iLine[1].split('T')[0],'%Y-%m-%d').date()
                
                while bDate < iDate:
                    bLine = base.readline().split()
                    bDate = dt.datetime.strptime(bLine[1].split('T')[0],'%Y-%m-%d').date()
                
                while bDate < oDate:
                    bLine = base.readline().split()
                    bDate = dt.datetime.strptime(bLine[1].split('T')[0],'%Y-%m-%d').date()
                    
                while oDate < bDate:
                    oLine = outp.readline().split()
                    oDate = dt.datetime.strptime(oLine[1].split('T')[0],'%Y-%m-%d').date()
            
            iXYZ = map(float,iLine[2:5])
            oXYZ = map(float,oLine[2:5])
            bXYZ = map(float,bLine[2:5])
            
            iX.append(iXYZ[0])
            iY.append(iXYZ[1])
            iZ.append(iXYZ[2])
            oX.append(oXYZ[0])
            oY.append(oXYZ[1])
            oZ.append(oXYZ[2])
            bX.append(bXYZ[0])
            bY.append(bXYZ[1])
            bZ.append(bXYZ[2])
            dates.append(toYearFraction(dt.datetime.strptime(iLine[1].split('T')[0],'%Y-%m-%d').date())-2000)
        
        f = open(folder+'/JPL_offsets/'+staName+'_jumps.dat')
        #f = open('folder/smallOffsetsRemoved/'+staName+'_jumps.dat')
        jumps = []
        for line in f:
            line2 = line.replace("\n","")
            jumps.append(toYearFraction(dt.datetime.strptime(line2,'%Y-%m-%dT%H:%M:%S').date())-2000)
        f.close()
        
        print(iX[200])
        
        a = (iX[-1]-iX[0])/(dates[-1]-dates[0])
        b = iX[0] - a*dates[0]
        for i in range(0,len(iX)):
            iX[i] = (1-a)*iX[i] - b
            oX[i] = (1-a)*oX[i] - b
            bX[i] = (1-a)*bX[i] - b
        
        a = (iY[-1]-iY[0])/(dates[-1]-dates[0])
        b = iY[0] - a*dates[0]
        for i in range(0,len(iY)):
            iY[i] = (1-a)*iY[i] - b
            oY[i] = (1-a)*oY[i] - b
            bY[i] = (1-a)*bY[i] - b
            
        a = (iZ[-1]-iZ[0])/(dates[-1]-dates[0])
        b = iZ[0] - a*dates[0]
        for i in range(0,len(iX)):
            iZ[i] = (1-a)*iZ[i] - b
            oZ[i] = (1-a)*oZ[i] - b
            bZ[i] = (1-a)*bZ[i] - b
        
        print(iX[200])
        
        pl.subplot(3,1, 1)
        pl.plot(dates,bX,'+y')
        pl.plot(dates,iX,'+b')
        pl.plot(dates,oX,',-r')
        pl.ylim(min(min(iX),min(oX),min(bX)),max(max(iX),max(oX),max(bX)))
        pl.xlim(min(dates),max(dates))
        for jump in jumps:
            pl.axvline(x=jump,color='b')
        pl.ylabel('x')
    
        pl.subplot(3, 1, 2)
        pl.plot(dates,bY,'+y')
        pl.plot(dates,iY,'+b')
        pl.plot(dates,oY,',-r')
        pl.ylim(min(min(iY),min(oY),min(bY)),max(max(iY),max(oY),max(bY)))
        pl.xlim(min(dates),max(dates))
        for jump in jumps:
            pl.axvline(x=jump,color='b')
        pl.ylabel('y')
    
        pl.subplot(3, 1, 3)
        pl.plot(dates,bZ,'+y')
        pl.plot(dates,iZ,'+b')
        pl.plot(dates,oZ,',-r')
        pl.ylim(min(min(iZ),min(oZ),min(bZ)),max(max(iZ),max(oZ),max(bZ)))
        pl.xlim(min(dates),max(dates))
        for jump in jumps:
            pl.axvline(x=jump,color='b')
        pl.ylabel('z')

        #pl.show()
        pl.savefig(folder+'/plots/'+staName+'.png',bbox_inches='tight')
        pl.close()
        
drawStuff()
