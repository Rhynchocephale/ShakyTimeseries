import sys
sys.path.insert(0,'interpolation')
import interp_spl as spl
import init_data as data
import ellipsoid
import os
import os.path
import numpy as np
import datetime
import time

def getDates(currentFolder,trusted):
    #GETS THE BEGINNING OF OBS FOR EARLIER STATION
    i = 0
    minStart = datetime.date.max
    maxEnd = datetime.date.min
    for sta in trusted:
        staFile = currentFolder + '/enu/' + sta + '_igs08_enu.dat'    #gets to the file
        f = open(staFile,'r')
        f.readline()    #skips header
        
        firstObsDay = datetime.datetime.strptime(f.readline().split()[1].split('T')[0],'%Y-%m-%d').date()   #gets string containing date & creates a date from it
        
        #reads all dates until the last
        for line in f:
            #print(line)
            currentDay = datetime.datetime.strptime(line.split()[1].split('T')[0],'%Y-%m-%d').date()
        
        f.close()
        
        if firstObsDay < minStart:
            minStart = firstObsDay
        if currentDay > maxEnd:
            maxEnd = currentDay
        i+=1
        
    print("beginning of obs: "+minStart.isoformat())
    print("end of obs:       "+maxEnd.isoformat())

    #CREATION OF A LIST CONTAING ALL THE DATES
    currentDate = minStart
    dateList = []
    while currentDate <= maxEnd:
        dateList.append(currentDate)
        currentDate += datetime.timedelta(days=1)

    return [minStart,maxEnd,dateList]


def remAvgRes(currentFolder):
    nbDaysCalc = 4 #nb of days taken into account for running avg (*2+1)
    
    #currentFolder = os.path.dirname(os.path.abspath(__file__))
    trusted = []
    for f in os.listdir(currentFolder):
        if os.path.isfile(os.path.join(currentFolder,f)) and f.endswith('.dat'):
            trusted.append(f.split('_')[0])
    trusted.sort()
    lenTrusted = len(trusted)
    fileLen = np.empty(lenTrusted)
   
    minStart,maxEnd,dateList = getDates(currentFolder,trusted)
   
    if not os.path.exists(currentFolder+'/corrections/'): #creates the /enu directory if it does not already exists
        os.makedirs(currentFolder+'/corrections/')
   
    nbSta = 0                                      #stations counter
    deltaT = (maxEnd-minStart).days+1              #number of days of obs
    enu = np.zeros(shape=(lenTrusted,deltaT,4))    #stores residuals. 3D tensor: x (e n u)+if faraway point (1 yes, 0 no), y stations, z days.
    pos = {} #records approximate position of every station
    for sta in trusted:  #for each trusted station                                   
        print('processing ref station: '+sta+' ('+"{:.0f}".format(nbSta+1)+'/'+"{:.0f}".format(lenTrusted)+')')
        staFile = currentFolder + '/enu/' + sta + '_igs08_enu.dat'    #gets to the file
        posPerDay = np.zeros(shape=(deltaT,3))              #stores pos in xyz for every day
        
        f = open(staFile,'r')
        f.readline()       #skips header
        firstObsDay = datetime.datetime.strptime(f.readline().split()[1].split('T')[0],'%Y-%m-%d').date() #gets string containing date & creates a date from it
        f.close()
        
        f2 = open(currentFolder+'/enu/param/'+sta+'.param','r')
        pos[sta] = map(float,f2.readline().split()[1:4])      #get pos
        f2.close()
        
        f = open(staFile,'r')
        f.readline()       #skips header
        for line in f: 
            line = line.split()
            currentDay = dateList.index(datetime.datetime.strptime(line[1].split('T')[0],'%Y-%m-%d').date())    #gets string containing date, creates a date from it, and gets the index of this date
            xyz = map(float,line[2:5])                                     #array with x,y,z in it
            posPerDay[currentDay,:] = xyz
        f.close()

        pos[sta] = xyz

        o = open(currentFolder+'/corrections/'+sta+'.enu','w')
        o.write("Rough_station_position_xyz" + " {:.0f} {:.0f} {:.0f}\n".format(xyz[0],xyz[1],xyz[2]))
        
        # DIFFERENCE BETWEEN ACTUAL POSITION AND AVERAGE OF THE STATION'S POSITION
        beginDate = dateList.index(firstObsDay)-1
        r=0
        currentDate = minStart
        while r < posPerDay.shape[0]:       #reads the new matrix line by line to compute avg
            row = posPerDay[r,:]

            if any(row) and r>beginDate:    #if not full of zeros and after beginning date
                k = 0
                tmpsum = np.zeros(3)
                valuesForMedian = []
                #with an absolute rejection level
                for j in range(max(r-nbDaysCalc,0),min(r+nbDaysCalc,posPerDay.shape[0])):  #makes sure not to get out of the matrix, takes nbDaysCalc lines after and nbDaysCalc before, when possible
                    if j!=r:
                        posPerDayLine = posPerDay[j,:]
                        if any(posPerDayLine):    #if not full of zeros
                            tmpsum += posPerDayLine
                            k += 1
                            valuesForMedian.append(posPerDayLine)
                    
                if k>3:
                    tmpsum /= float(k)
                    #tmpsum = (tmpsum*k+row)/(k+1) #mean
                    tmpsum = np.median(valuesForMedian,axis=0) #median
                    enu[nbSta,r,0:3] = row-tmpsum
 
                o.write(currentDate.isoformat()+" {:.11f} {:.11f} {:.11f} ".format(enu[nbSta,r,0],enu[nbSta,r,1],enu[nbSta,r,2 ])+" {:.11f} {:.11f} {:.11f} ".format(row[0],row[1],row[2])+'\n')
            else:
                if enu[nbSta,r,3] == 0 and r>0:
                    enu[nbSta,r,3] = 2 #if !=0, then not taken into account when splining
            
            currentDate += datetime.timedelta(days=1)
            r+=1
        print('gaps in obs: '+str(np.count_nonzero(enu[nbSta,:,3]==2)))
        print('outliers: '+str(np.count_nonzero(enu[nbSta,:,3]==1))+'\n')
        nbSta+=1
        o.close()
  
    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(currentFolder + '/enu/'):
        if os.path.isfile(os.path.join(currentFolder + '/enu/',f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    if not os.path.exists(currentFolder+'/smoothed/'): #creates the /smoothed directory if it does not already exists
        os.makedirs(currentFolder+'/smoothed/')
    
    u=0
 
    staCount = 0
    for sta in fileList:
        staName = sta.split('_')[0]
        print('processing obs station: '+staName+' ('+"{:.0f}".format(staCount+1)+'/'+"{:.0f}".format(len(fileList))+')')
                
        absSta = currentFolder+'/enu/'+sta
                        
        f = open(absSta,'r')
        f.readline()       #skips header
        g = open(currentFolder+'/smoothed/'+ sta.split('.')[0]+'_smoothed.dat','w')
        g.write('name\tepoch\tx\ty\tz\n')
        for line in f:            
            line = line.split()
            xyz = map(float,line[2:5])
            today = datetime.datetime.strptime(line[1].split('T')[0],'%Y-%m-%d').date()
            if today >= minStart:
                i = dateList.index(today)
            else:
                i = -1
            isWritten = False
            if any(xyz):                   
                if i>=0 and i<deltaT:
                    corr = enu[staCount,i,:]
                    if (corr[3]== 0):
                        xyz2 = xyz-corr[0:3]
                        '''if any(corr > 20): #if correction > 1cm
                            print("too big correction:")
                            print(today)
                            print(corr[0:3])
                            xyz2 == xyz'''
                        g.write(line[0]+'\t'+line[1]+'\t'+"{:.4f}".format(xyz2[0])+'\t'+"{:.4f}".format(xyz2[1])+'\t'+"{:.4f}".format(xyz2[2])+'\n')
                        isWritten = True
                else:
                    u+=1
            if not isWritten:
                g.write(line[0]+'\t'+line[1]+'\t0.0000\t0.0000\t0.0000\n')
        staCount+=1
        
    return enu
