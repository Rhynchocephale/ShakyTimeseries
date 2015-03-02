import os
import numpy as np
import stn_pred_model as stn
#from remAvgRes_withoutSplines import *
import datetime as dt
import shutil

def toText(dates,data,fileName,staName):
    f = open(fileName,'w')
    i=0
    f.write('name\tepoch\tx\ty\tz\n')
    for date in dates:
        f.write(staName+'\t'+date.strftime('%Y-%m-%d')+'T12:00:00\t'+"{:.6f}\t{:.6f}\t{:.6f}\n".format(data[i,0],data[i,1],data[i,2]))
        i+=1
    f.close()
    return 0
    
     
def getxyz(fileName):
    f = open(fileName,'r')
    f.readline() #header
    xyz = map(float,f.readline().split()[2:5])
    return xyz

#used for jump detection only. Everything after that will be done with the real xyz timeseries
def getResiduals(folder):
    
    if os.path.exists(folder+'/residuals/'): #clears the /residuals directory if it already exists
        shutil.rmtree(folder+'/residuals/')
    os.makedirs(folder+'/residuals/')
        
    if os.path.exists(folder+'/xml/'): #clears the /xml directory if it already exists
        shutil.rmtree(folder+'/xml/')
    os.makedirs(folder+'/xml/')    
    
    #WITH SMOOTHED
    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder+'/smoothed/xyz/'):
        if os.path.isfile(os.path.join(folder+'/smoothed/xyz/',f)) and f.endswith('.dat'):
            fileList.append(f.split('_')[0])
    fileList.sort()
    '''
    #WITH RAW
    fileList = [x+'_igs08_xyz.dat' for x in listOfStations]
    fileList.sort()
    '''
    i=1
    for staName in fileList:
        print('computing residuals: '+staName+' ('+str(i)+'/'+str(len(fileList))+')')
        xyz = getxyz(folder+'/'+staName+'_igs08_xyz.dat')
        m = stn.model(station=staName, xyz=xyz, filename=folder+'/xml/'+staName+'.xml', loadfile=False)
        
        m.loadTimeSeries(filename=folder+'/smoothed/xyz/'+staName+'_igs08_xyz_smoothed.dat')
        m.fitAllLinear()
        m.save()
        dates,enu,useobs = m.getObs()
        enu_calc = m.calc(m.dates)
        residuals = 1000*(enu-enu_calc)
        toText(dates,residuals,folder+'/residuals/'+staName+'_residuals.dat',staName)
        i+=1

def addResiduals2Model(folder):
    if os.path.exists(folder+'/output/'): #clears the /residuals directory if it already exists
        shutil.rmtree(folder+'/output/')
    os.makedirs(folder+'/output/')
    
    #BROWSE ALL FILES IN FOLDER
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder+'/withSlowSlips/'):
        if os.path.isfile(os.path.join(folder+'/withSlowSlips/',f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    p=1
    for fileName in fileList:
        staName = fileName.split('_')[0]
        print('Final output: '+staName+' ('+str(p)+'/'+str(len(fileList)))
        
        f = open(folder+'/withSlowSlips/'+fileName,'r')
        g = open(folder+'/output/'+staName+'_igs08_enu.dat','w')
        xyz = getxyz(folder+'/smoothed/xyz/'+fileName)
        m = stn.model(station=staName, xyz=xyz, filename=folder+'/xmlWithSlowSlips/'+staName+'.xml', loadfile=True)
        m.loadTimeSeries(filename=folder+'/withSlowSlips/'+fileName)
        dates,enu,useobs = m.getObs()
        enu_calc = m.calc(m.dates)
        
        g.write('name epoch               x             y           z            flag\n')
        f.readline() #header
        currentDate = dt.date.min
        i = 0
        while i<len(dates):
            
            #if desync
            while currentDate < dates[i]:
                line = f.readline().split()
                currentDate = dt.datetime.strptime(line[1].split('T')[0],'%Y-%m-%d').date()
            while dates[i] < currentDate:
                i+=1
                
            xyz = map(float,line[2:5]) + enu_calc[i]
                
            g.write(line[0]+' '+line[1]+' '+"{:.5f}".format(xyz[0])+' '+"{:.5f}".format(xyz[1])+' '+"{:.5f}".format(xyz[2])+'\n')
            i+=1
            
        f.close()
        g.close()
        p+=1

def putOffsets(folder,removeSmall=False):
    
    #BROWSE ALL FILES IN FOLDER
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder+'/smoothed/'):
        if os.path.isfile(os.path.join(folder+'/smoothed/',f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    #fileList = [x+'_igs08_xyz_smoothed.dat' for x in listOfStations] #smoothed
    #fileList = [x+'_igs08_xyz.dat' for x in listOfStations] #raw
    #fileList.sort()
    
    if os.path.exists(folder+'/withOffsets/'): #clears the /withOffsets directory if it already exists
        shutil.rmtree(folder+'/withOffsets/')
    os.makedirs(folder+'/withOffsets/')
    
    if os.path.exists(folder+'/xmlWithOffsets/'): #clears the /xmlWithOffsets directory if it already exists
        shutil.rmtree(folder+'/xmlWithOffsets/')
    os.makedirs(folder+'/xmlWithOffsets/')
    
    if removeSmall:
        if os.path.exists(folder+'/smallOffsetsRemoved/'): #clears the /withOffsets directory if it already exists
            shutil.rmtree(folder+'/smallOffsetsRemoved/')
        os.makedirs(folder+'/smallOffsetsRemoved/')


    subfolder = '/smoothed/xyz/' #smoothed
    #subfolder = '/' #raw
    #oneMonthAndAHalf = dt.timedelta(days=45)
    i=1
    for fileName in fileList:
        staName = fileName.split('_')[0]
        print("==================\nFitting model for: "+staName+' ('+str(i)+'/'+str(len(fileList))+')\n==================')
        
        xyz = getxyz(folder+subfolder+fileName.replace('enu','xyz'))
        m = stn.model(station=staName, xyz=xyz)
        m.loadTimeSeries(filename=folder+subfolder+fileName.replace('enu','xyz'))
        m.fitAllLinear()
        
        #adds offsets from file
        f = open(folder+"/JPL_offsets/"+staName+"_jumps.dat",'r')
        offsetDates = f.readlines()
        
        deleted = 0
        previousDate = dt.date.min
        for date in offsetDates:
            if date != "\n":
                properDate = date.split('\n')[0] #+'T12:00:00'
                dtDate = dt.datetime.strptime(properDate,'%Y-%m-%dT%H:%M:%S').date()
                #if dtDate >= previousDate+oneMonthAndAHalf:
                if True:
                    offset = stn.equipment_offset(m,properDate)
                    m.addComponent(offset)
                    #previousDate = dtDate
                #else:
                    #deleted += 1
                    
        f = open(folder+"/JPL_offsets/"+staName+"_slowslips.dat",'r')
        offsetDates = f.readlines()
        
        #adds slowslips from file
        deleted = 0
        previousDate = dt.date.min
        for date in offsetDates:
            if date != "\n":
                properDate = date.split('\n')[0].split() #+'T12:00:00'
                dtDate1 = dt.datetime.strptime(properDate[0],'%Y-%m-%dT%H:%M:%S').date()
                dtDate2 = dt.datetime.strptime(properDate[1],'%Y-%m-%dT%H:%M:%S').date()
                #if dtDate >= previousDate+oneMonthAndAHalf:
                if True:
                    slowSlip = stn.slow_slip(m,date=dtDate1,duration=(dtDate1-dtDate2).days)
                    m.addComponent(slowSlip)
                    #previousDate = dtDate
        
        m.fitAllLinear()
        m.save(folder+'/xmlWithOffsets/'+staName+'.xml')
        dates = m.dates
        #dates,enu,useobs = m.getObs()
        enu_calc = m.calc(dates,enu=False)
        toText(dates,enu_calc,folder+'/withOffsets/'+staName+'_igs08_xyz.dat',staName)
        
        if removeSmall:
            #iterationNb = 1
            #currentDeleted = 1
            initialComponentsNb = len(np.nonzero([x.componentType() in ["tectonic_offset", "equipment_offset"] for x in m.components])[0]) + deleted
            dateOfPrevious = dt.date.min
            print('Removing smallest offsets')
            
            #removes the offsets smaller than noise / 1mm
            #noise = stn.model.robustStandardError(enu_calc)

            for c in m.components:
                if c.componentType() in ("tectonic_offset", "equipment_offset"):
                    '''if c.eventDate().date() < dateOfPrevious+oneMonthAndAHalf: #if two offsets at the same time. Happens sometimes, no idea why.
                        m.removeComponent(c)
                        deleted+=1
                        continue
                    dateOfPrevious = c.eventDate().date()'''
                    values = c.paramValues()
                
                    isNull = 0
                    for u in [0,1,2]:
                        #if values[u+1] < noise[u]*2:
                        if values[u+1] < 0.02:
                            #c.setComponent(u,0,fixed=False)
                            isNull += 2*u+1
                            '''c.parameters[u].setFitValue(0)
                            c.parameters[u].setFixed(True)
                            print('1 offset parameter to 0')
                            if u>0:
                            print(c.parameters)
                            '''
                    if isNull == 9: #if e&n&u == 0
                        m.removeComponent(c)
                        deleted += 1
                    elif isNull in [6,8,5] : #if u == 0 or e&u == 0 or n&u == 0
                        c.setComponent(2,0,fixed=True)
                    elif isNull == 4: #if e&n == 0
                        c.setComponent(0,0,fixed=True)
                        c.setComponent(1,0,fixed=True)
            
            g = open(folder+'/smallOffsetsRemoved/'+staName+'_jumps.dat','w')
            for c in m.components:
                if c.componentType() in ("tectonic_offset", "equipment_offset"):
                    g.write(dt.datetime.strftime(c.eventDate().date(),'%Y-%m-%dT%H:%M:%S')+'\n')
            g.close()

            print('Deleted offsets: '+str(deleted)+'/'+str(initialComponentsNb))
            
            m.fitAllLinear()
            m.save(folder+'/xmlWithOffsets/'+staName+'.xml')
            #dates,enu,useobs = m.getObs()
            enu_calc = m.calc(dates,enu=False)
            toText(dates,enu_calc,folder+'/withSmallOffsetsRemoved/'+staName+'_igs08_xyz.dat',staName)
        i+=1

def sortFile(fileName):
        #ORDERS FILE ALPHABETICALY
        unsortedFile = open(fileName, "r")
        sortedData = sorted(unsortedFile, key = str.lower) #sorts the file, stores it
        os.remove(fileName) #removes file
        sortedFile = open(fileName,"w") #recreates file
        for line in sortedData:
            sortedFile.write(line)
        sortedFile.close()

def getSlowSlips(folder):
    
    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder+'/withOffsets/'):
        if os.path.isfile(os.path.join(folder+'/withOffsets/',f)) and f.endswith('.dat'):
            fileList.append(f)
    
    #fileList = [x+'_igs08_xyz.dat' for x in listOfStations]
    fileList.sort()
    
    if os.path.exists(folder+'/withSlowSlips/'): #clears the /output directory if it already exists. Creates it otherwise
        shutil.rmtree(folder+'/withSlowSlips/')
    os.makedirs(folder+'/withSlowSlips/')
    
    if os.path.exists(folder+'/xmlWithSlowSlips/'): #clears the /output directory if it already exists. Creates it otherwise
        shutil.rmtree(folder+'/xmlWithSlowSlips/')
    os.makedirs(folder+'/xmlWithSlowSlips/')
    
    fourMonths = dt.timedelta(days=122)
    
    staNb=1
    oneMonth=dt.timedelta(days=30)
    for fileName in fileList:
        print(fileName+': '+str(staNb)+'/'+str(len(fileList)))
        staName = fileName.split('_')[0]
        #sortFile(folder+'/offsets/'+staName+"_jumps.dat") #sorts offsets by name
        xyz = getxyz(folder+'/withOffsets/'+fileName)
        m = stn.model(station=staName, xyz=xyz, filename=folder+'/xmlWithOffsets/'+staName+'.xml', loadfile=True)
        m.loadTimeSeries(filename=folder+'/withOffsets/'+fileName)
        
        offsetList=[]
        for c in m.components:
            if c.componentType() in ("tectonic_offset", "equipment_offset"):
                offsetList.append([c.paramValues(),c.eventDate(),c])

        
        i=1
        while i<len(offsetList)-1:
            j=i
            print('Offset '+str(i)+'/'+str(len(offsetList)-1))
            while (offsetList[j+1][1]-offsetList[j][1]) < fourMonths:
                offset_j = offsetList[j]
                offset_jplus1 = offsetList[j+1]
                if offset_j[1].date() == offset_jplus1[1].date(): #if two offsets at the same date. Happens sometimes, no idea why.
                    m.removeComponent(offset_jplus1[2])
                    m.save(folder+'/xmlWithSlowSlips/'+staName+'.xml')
                print(offset_jplus1[1])
                #if (cmp(offset_j[0][0],0) in [cmp(offset_jplus1[0][0],0),0] and cmp(offset_j[0][1],0) in [cmp(offset_jplus1[0][1],0),0]): #if offsets in same direction (or no offset)
                if ((cmp(offset_j[0][0],0) in [cmp(offset_jplus1[0][0],0),0]) or (abs(offset_j[0][0]) > 5*abs(offset_jplus1[0][0])) or (abs(offset_jplus1[0][0]) > 5*abs(offset_j[0][0]))) and ((cmp(offset_j[0][1],0) in [cmp(offset_jplus1[0][1],0),0]) or (abs(offset_j[0][1]) > 5*abs(offset_jplus1[0][1])) or (abs(offset_jplus1[0][1]) > 5*abs(offset_j[0][1]))):
                    j+=1
                    if j>=len(offsetList)-1:
                        break
                else:
                    break
            
            if j-i > 1:
                m2 = m.copy()
                m2.loadTimeSeries(filename=folder+'/withOffsets/'+fileName)
                #slowSlip = stn.slow_slip_ramp(m2,date=offsetList[i][1],end_date=offsetList[j][1]) #slow slip ramp
                slowSlip = stn.slow_slip(m2,date=offsetList[i][1],duration=(offsetList[j][1]-offsetList[i][1]).days) #slow slip with error function
                m2.addComponent(slowSlip)
                for row in offsetList[i:j+1]:
                    m2.removeComponent(row[2])
                m2.fitAllLinear()
            
                if isbetter(m2,m,m.getObs()[1],offsetList[i][1]-oneMonth,offsetList[j][1]+oneMonth):
                #if True:
                    del offsetList[i:j+1]
                    m = m2.copy()
                    m.loadTimeSeries(filename=folder+'/withOffsets/'+fileName)
                    m.save(folder+'/xmlWithSlowSlips/'+staName+'.xml')
                    print('SLOW SLIP DETECTED')
                    i-=1
            i+=1
        
        dates = m.dates
        #dates,enu,useobs = m.getObs()
        enu_calc = m.calc(dates,enu=False)
        toText(dates,enu_calc,folder+'/withSlowSlips/'+staName+'_igs08_xyz.dat',staName)
        staNb+=1

def isbetter(m2,m,obs,beginDate,endDate):

    all_days = m.dates
    days = []
    for day in all_days:
        if day > endDate:
            break
        if day >= beginDate:
            days.append(day)
    
    enu_m = m.calc(days)
    enu_m2 = m2.calc(days)
    
    dist_m = 0
    dist_m2 = 0
    for i in range(0,len(days)):
        dist_m += np.power(abs(enu_m[i]-obs[i]),2)
        dist_m2 += np.power(abs(enu_m2[i]-obs[i]),2)
    
    return np.linalg.norm(dist_m[0:2]) > np.linalg.norm(dist_m2[0:2])
