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
import coordConv


def makeAllTrusted( fileName ):
    currentFolder = os.path.dirname(os.path.abspath(__file__))
    
    #BROWSE ALL FILES IN FOLDER TO GET THEIR NAME
    #only takes the files, not the folders, and checks that the file is a .dat
    staList = []
    for f in os.listdir(currentFolder + '/timeseries'):
        if os.path.isfile(os.path.join(currentFolder + '/timeseries',f)) and f.endswith('.dat'):
            staList.append(f.split('_')[0])  #only keeps the 4-letter name of station
    
    g = open(fileName,'w')
    for sta in staList:
        g.write(sta+'\t')
    g.close()
    
    
'''
def kriging(enu,pos,step,dateList):
    X,Y,Z,X0,Y0 = data.init_data(enu,pos,step,0,0) #one first data init, just to get the size of X0 
    T = np.zeros((len(dateList),3))  #one grid per day per dimension
    for i in range(0,len(dateList)): #for each day
        print('processing day: '+' ('+"{:.0f}".format(i+1)+'/'+"{:.0f}".format(len(dateList))+')')
        for j in [0,1,2]: #for each dimension
            X,Y,Z,X0,Y0 = data.init_data(enu,pos,step,j,i)
            T = np.asarray(krg.interp_krg(X,Y,Z,X0,Y0,6000))

    return T
'''

def splines(enu,pos,dateList,currentFolder):

    posShape0 = len(pos)
    lenDateList = len(dateList)
    spl_param = np.zeros((posShape0+3,lenDateList,3))  #one grid per day per dimension
    nbOfUsedStations = np.zeros(lenDateList)
    modDateList = np.mod(lenDateList-1,13)
    for i in range(0,lenDateList): #for each day
        if np.mod(i,13)==modDateList:
            print('processing day: '+"{:.0f}".format(i+1)+'/'+"{:.0f}".format(lenDateList))
            #print('processing days: '+"{:.3f}".format(float(i+1)*100/lenDateList)+'%')
        
        for j in [0,1,2]: #for each dimension

            X,Y,Zalti,Zinterp,indices = data.init_data_nogrid(enu,pos,j,i)
            if X.size > 6: #if enough stations
                spli = np.asarray(spl.interp_spl(X,Y,Zalti,Zinterp,10^3))
                for indice in indices:
                    spli = np.insert(spli,indice+3,0)
            else:
                spli = np.zeros(posShape0+3)
            spl_param[:,i,j] = spli
        nbOfUsedStations[i] = len(indices) #counts the zeros
    
    nbOfUsedStations = posShape0 - nbOfUsedStations
    
    #np.save('splines.npy',spl_param)
    
    return [spl_param,nbOfUsedStations]

def corrSpline(enu,staCount,position,dist,spl,day,errors):
    enu_corr = np.zeros(3)
    K = dist*np.log(np.sqrt(dist))
    problem = 0   
    for dim in [0,1,2]:
        spldaydim = spl[:,day,dim]
        if any(spldaydim):
            enu_corr[dim] = spldaydim[0] + spldaydim[1]*position[0] + spldaydim[2]*position[1] + sum(spldaydim[3:]*K)
        else: #if no correction, running avg. COULD HAPPEN, BUT SHOULD BE RARE.
            enu_corr[dim] = enu[staCount,day,dim]
            problem = 1
    if problem:
        errors += 1
    return enu_corr,errors

def inverseDist(dist,enu,day):
    
    #first, compute the inverse of distances with other stations. If dist=0, then 1/dist=1. Then, ponderate the errors with it
    rowNb = 0
    corr_enu = np.zeros(3)
    sumInvDist = 0
    for row in dist:
        if row != 0 and any(enu[rowNb,day,0:3]):
            corr_enu += enu[rowNb,day,0:3]*row
            sumInvDist += row
            rowNb += 1
    
    if sumInvDist>0:
        corr_enu = corr_enu/sumInvDist
    
    return corr_enu,rowNb
    
def getDates(currentFolder,trusted):
    #GETS THE BEGINNING OF OBS FOR EARLIER STATION
    i = 0
    minStart = datetime.date.max
    maxEnd = datetime.date.min
    for sta in trusted:
        staFile = currentFolder + '/timeseries/' + sta + '_igs08_xyz.dat'    #gets to the file
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



def remAvgRes_nb( trustedList ):
    
    currentFolder = os.path.dirname(os.path.abspath(__file__))
    f = open(trustedList,'r')
    trusted = f.readline().split()                          #returns list of names of trusted stations
    f.close()
    lenTrusted = len(trusted)
    fileLen = np.empty(lenTrusted)
   
    minStart,maxEnd,dateList = getDates(currentFolder,trusted)
   
    if not os.path.exists(currentFolder+'/timeseries/enu/'): #creates the /enu directory if it does not already exists
        os.makedirs(currentFolder+'/timeseries/enu/')
   
    nbSta = 0                                      #stations counter
    deltaT = (maxEnd-minStart).days+1              #number of days of obs
    enu = np.zeros(shape=(lenTrusted,deltaT,4))    #stores residuals. 3D tensor: x (e n u)+if faraway point (1 yes, 0 no), y stations, z days.
    pos = {} #records approximate position of every station
    for sta in trusted:  #for each trusted station
        asd=0                                    
        print('processing ref station: '+sta+' ('+"{:.0f}".format(nbSta+1)+'/'+"{:.0f}".format(lenTrusted)+')')
        staFile = currentFolder + '/timeseries/' + sta + '_igs08_xyz.dat'    #gets to the file
        #nbLines = lineCount(staFile)-1                     #counts the lines in the file, minus header
        #llh = np.zeros(shape=(deltaT,3))                   #stores lon, lat, h
        posPerDay = np.zeros(shape=(deltaT,3))              #stores pos in xyz for every day
        
        f = open(staFile,'r')
        f.readline()       #skips header
        firstObsDay = datetime.datetime.strptime(f.readline().split()[1].split('T')[0],'%Y-%m-%d').date() #gets string containing date & creates a date from it
        f.close()
        
        f = open(staFile,'r')
        f.readline()       #skips header
        for line in f: 
            line = line.split()
            currentDay = dateList.index(datetime.datetime.strptime(line[1].split('T')[0],'%Y-%m-%d').date())    #gets string containing date, creates a date from it, and gets the index of this date
            xyz = map(float,line[2:5])                                     #array with x,y,z in it
            posPerDay[currentDay,:] = xyz
            #lon,lat,h = ellipsoid.grs80.geodetic(xyz)
            #llh[currentDay,:] = [lon,lat,h,]
        f.close()
                
        #pos[nbSta,:] = [merc[0],merc[1],merc[2],lon,lat,h]
        #pos[nbSta,:] = [merc[0],merc[1],merc[2]]
        merc = coordConv.xyz2trans(xyz)
        pos[sta] = merc

        o = open(currentFolder+'/timeseries/enu/'+sta+'.enu','w')
        o.write("Rough_station_position_xyz" + " {:.0f} {:.0f} {:.0f}\n".format(xyz[0],xyz[1],xyz[2]))
        
        # DIFFERENCE BETWEEN ACTUAL POSITION AND AVERAGE OF THE STATION'S POSITION
        beginDate = dateList.index(firstObsDay)-1
        r=0
        currentDate = minStart
        nbDaysCalc = 8 #nb of days taken into account for running avg (*2+1)
        while r < posPerDay.shape[0]:       #reads the new matrix line by line to compute avg
            row = posPerDay[r,:]
            if any(row) and r>beginDate:    #if not full of zeros and after beginning date
                k = 0
                tmpsum = np.zeros(3)               
                neighbours = np.zeros(3)
                #based on number of neighbours
                for j in range(max(r-nbDaysCalc,0),min(r+nbDaysCalc,posPerDay.shape[0])): 
                        posPerDayLine = posPerDay[j,:]
                        if any(posPerDayLine) and enu[nbSta,r,3] == 0:    #if not full of zeros and valid point
                            tmpsum += posPerDayLine
                            k += 1
                            neighbours += abs(row-posPerDayLine) < 0.005
                        
                if k>3:
                    if all(neighbours >= k/4) and sum(neighbours)>=3*k/4:
                        if enu[nbSta,r,3] == 0:
                            tmpsum /= k
                            enu[nbSta,r,0:3] = row-tmpsum
                    else:
                        enu[nbSta,r,3] = 1
                        r-=nbDaysCalc+1
                    
                o.write(currentDate.isoformat()+" {:.11f} {:.11f} {:.11f} ".format(enu[nbSta,r,0]*1000,enu[nbSta,r,1]*1000,enu[nbSta,r,2]*1000)+" {:.11f} {:.11f} {:.11f} ".format(row[0],row[1],row[2])+'\n')
            else:
                if enu[nbSta,r,3] == 0 and r>0:
                    enu[nbSta,r,3] = 2 #if !=0, then not taken into account when splining
            
            currentDate += datetime.timedelta(days=1)
            r+=1
        print(np.count_nonzero(enu[nbSta,:,3]==1))
        nbSta+=1
        o.close()
        
    ''' 
#METHOD USED HERE-------------------------------------------------
    # AVERAGE OF THE ERRORS PER DAY
    enu_corr = np.zeros(shape=(deltaT,4))
    for i in range(enu.shape[1]):       #per day
        tmpsum = np.zeros(3)
        k=0
        for j in range(lenTrusted):   #per station
            line = enu[j,i,0:3]
            if any(line):               #if not full of zeros
                tmpsum+=line
                k+=1
        if k>0:
            tmpsum /= float(k)
            enu_corr[i,:] = [tmpsum[0],tmpsum[1],tmpsum[2],k]
#------------------------------------------------------------------
    '''
    
    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(currentFolder + '/timeseries'):
        if os.path.isfile(os.path.join(currentFolder + '/timeseries',f)) and f.endswith('.dat'):
            fileList.append(f)
    
    if not os.path.exists(currentFolder+'/timeseries/smoothed/'): #creates the /smoothed directory if it does not already exists
        os.makedirs(currentFolder+'/timeseries/smoothed/')
    
    #ZZ0 = kriging(enu,pos,10000,dateList)
    spl_param,nbOfUsedStations = splines(enu,pos,dateList,currentFolder)
 
    u=0
 
    staCount = 0
    for sta in fileList:
        staName = sta.split('_')[0]
        print('processing obs station: '+staName+' ('+"{:.0f}".format(staCount+1)+'/'+"{:.0f}".format(len(fileList))+')')
                
        absSta = currentFolder+'/timeseries/'+sta
        
        #IF USING DISTANCES---------------------------------------------------------------
        position = pos[staName]

        dist = np.zeros(len(pos))       
        distRowNb = 0
        x=position[0]
        y=position[1]
        z=position[2]
        for key in pos:
            distance = np.sqrt( (x-pos[key][0])**2 + (y-pos[key][1])**2 + (z-pos[key][2])**2 )
            dist[distRowNb] = (distance,3)[distance==0] #result = (on_false,on_true)[condition]
            distRowNb += 1
        #---------------------------------------------------------------------------------
        
                
        f = open(absSta,'r')
        f.readline()       #skips header
        g = open(currentFolder+'/timeseries/smoothed/'+ sta.split('.')[0]+'_smoothed.dat','w')
        g.write('name epoch               x             y           z            flag stations_for_smoothing\n')
        errors = 0
        for line in f:            
            line = line.split()
            xyz = map(float,line[2:5])
            today = datetime.datetime.strptime(line[1].split('T')[0],'%Y-%m-%d').date()
            if today >= minStart:
                i = dateList.index(today)
            else:
                i = -1
            if any(xyz):                   
                if i>=0 and i<deltaT:
                    #lon,lat,h = ellipsoid.grs80.geodetic(xyz)
                    #enu_corr,staNb = inverseDist(dist,enu,i)
                    if (enu[staCount,i,3] == 0):
                        #enu_corr,errors = corrSpline(enu,staCount,position,dist,spl_param,i,errors)
                        xyz_corr,errors = corrSpline(enu,staCount,position,dist,spl_param,i,errors)
                        nbOfStationsUsedToday = nbOfUsedStations[i]
                        #nbOfStationsUsedToday = enu_corr[i,3]
                        #enu_pos = coordConv.xyz2trans(xyz)
                        #enu_pos_corr = enu_pos - enu_corr
                        xyz_pos_corr = xyz - xyz_corr
                        #xyz2 = coordConv.trans2xyz([enu_pos_corr[0],enu_pos_corr[1],enu_pos_corr[2]])
                        xyz2 = xyz_pos_corr
                        xyz_diff = np.array(xyz) - np.array(xyz2)
                        if any(xyz_diff > 0.03): #if correction > 3cm
                            xyz2 = np.array(xyz) - np.array(enu[staCount,i,0:3])
                            print("too big correction:")
                            print(xyz_diff)
                            
                        g.write(line[0]+' '+line[1]+' '+"{:.4f}".format(xyz2[0])+' '+"{:.4f}".format(xyz2[1])+' '+"{:.4f}".format(xyz2[2])+' '+line[5]+' '+"{:.0f}".format(nbOfStationsUsedToday)+'\n')
                    #else:
                        #xyz2 = xyz - enu[staCount,i,0:3]  #if bad point, back to avg and that's all
                    
                    #if np.random.rand()<0.0005:
                    #    print((np.array(xyz2)-np.array(xyz))*1000)
                else:
                    u+=1
            
        print("days without enough stations: "+str(errors))
        staCount+=1
        
    return enu

'''
t = time.time()
remAvgRes('/home/cdrouadaine/timeseries/StationsWeTrust.txt')
elapsed = time.time() - t
print(elapsed)
print(elapsed/222)
'''
