import datetime as dt
import drawStuff
import pylab as pl

def plot1file(filePath):
    f=open(filePath)
    f.readline()
    myX = []
    myY = []
    myZ = []
    myDates = []
    for line in f:
        myLine = line.split()
        myDate = dt.datetime.strptime(myLine[1].split('T')[0],'%Y-%m-%d').date()
        myXYZ = map(float,myLine[2:5])
        myX.append(myXYZ[0])
        myY.append(myXYZ[1])
        myZ.append(myXYZ[2])
        myDates.append(drawStuff.toYearFraction(myDate)-2000)
    f.close()
    
    pl.subplot(3,1, 1)
    pl.plot(myDates,myX,'-b')
    pl.ylim(min(myX),max(myX))
    pl.xlim(min(myDates),max(myDates))
    pl.ylabel('x')
    
    pl.subplot(3,1, 2)
    pl.plot(myDates,myY,'-b')
    pl.ylim(min(myY),max(myY))
    pl.xlim(min(myDates),max(myDates))
    pl.ylabel('y')
    
    pl.subplot(3,1, 3)
    pl.plot(myDates,myZ,'-b')
    pl.ylim(min(myZ),max(myZ))
    pl.xlim(min(myDates),max(myDates))
    pl.ylabel('z')
    
    pl.savefig('/home/cdrouadaine/timeseries/timeseries/GISBplot.png')
    pl.show()
    
plot1file('/home/cdrouadaine/timeseries/timeseries/GISB_igs08_xyz.dat')
#plot1file('/home/cdrouadaine/timeseries/timeseries/smoothed/xyz/BLUF_igs08_xyz_smoothed.dat')
    
