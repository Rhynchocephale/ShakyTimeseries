import numpy as np
import datetime as dt
import time
import os
import stn_pred_model as pred
import usePredModel as use
import shutil
import scipy.stats
import sys


def JPL_STP1(currentFolder,hole=True):

    minimumDaysBetweenConsecutiveEarthquakes = int(20)
    minimumSignificance = 0.999
    windowSizes = [640,320,160,80,40,20]
    #windowSizes = [25]
    minimumJumpSizePlani = 2 #in millimetres
    minimumJumpSizeAlti = 4
    minimumJumpSizes = np.array([minimumJumpSizePlani,minimumJumpSizeAlti])
    minimumDaysForASlowSlip = 4

    inputFolder = currentFolder + '/residuals/'
    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(inputFolder):
        if os.path.isfile(os.path.join(inputFolder,f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    if os.path.exists(currentFolder+'/xmlWithJPLOffsets/'): #clears the /xml directory if it already exists
        shutil.rmtree(currentFolder+'/xmlWithJPLOffsets/')
    os.makedirs(currentFolder+'/xmlWithJPLOffsets/')
    '''
    if os.path.exists(currentFolder+'/significance/'): #clears the /xml directory if it already exists
        shutil.rmtree(currentFolder+'/significance/')
    os.makedirs(currentFolder+'/significance/')
    '''
    outputFolder = currentFolder + '/JPL_offsets/'
    if os.path.exists(outputFolder): #clears the directory if it already exists
        shutil.rmtree(outputFolder)
    os.makedirs(outputFolder)
    
    lenFileList = '/'+str(len(fileList))+')'
    
    fileNb = 1
    twelveHours = dt.timedelta(hours=12)
    jumpMagnitude = [0,0]
    for fileName in fileList:
        staName = fileName.split('_')[0]
        
        print("=======================\nBEGINNING OF COMPUTATION FOR "+staName+' ('+str(fileNb)+lenFileList)
        subfolder = '/smoothed/xyz/'
        xyz = use.getxyz(currentFolder+subfolder+staName+'_igs08_xyz_smoothed.dat')
        m = pred.model(station=staName, xyz=xyz)
        m.loadTimeSeries(filename=currentFolder+subfolder+staName+'_igs08_xyz_smoothed.dat')
        m.fitAllLinear()
        iterationNb = 1
        overallJumps = [] #list of all jumps found so far
        overallSlowSlips = []
        
        f = open(currentFolder+'/'+staName+'_igs08_xyz.dat')
        f.readline() #header
        positionsForNoise = []
        for line in f:
            positionsForNoise.append( map(float,line.split()[2:5]) )
        f.close()
        
        errors = m.robustStandardError(np.array(positionsForNoise))
        
        isFirstDay = True
        dateList=[]
        daysToCompute_original=[] #if =1, we'll compute. If =0, skip
        
        f = open(inputFolder+fileName)
        f.readline() #header
            
        #reads all lines in the file, gets the dates, for later use
        for line in f:
            splitted = line.split()
            currentDate = int(pred.asday(dt.datetime.strptime(splitted[1].split('T')[0],'%Y-%m-%d').date()))
            if isFirstDay:
                prevDate = currentDate-1
                isFirstDay = False
            while currentDate-prevDate > 1:
                prevDate += 1
                dateList.append(prevDate)
                daysToCompute_original.append(0)
            prevDate += 1
            dateList.append(currentDate)
            if any(map(float,splitted[2:5])):
                daysToCompute_original.append(1)
            else:
                daysToCompute_original.append(0)
        f.close()
        
        
        #for halfLengthOfRamp in [1.0/60,1.0/20,1.0/10,1.0/5]:
        for halfLengthOfRamp in [1.0/60,1.0/60,1.0/10]:
            isNotEnd = True
            daysToCompute = list(daysToCompute_original)
            while isNotEnd: #finds biggest jumps, add them to model, recompute residuals, and again, and again, until jumps too uncertain
            
                print("ITERATION "+str(iterationNb)+' (STATION '+str(fileNb)+lenFileList)
                print('Reading file')
    
                dailyPos=np.zeros((len(dateList),3))
                f = open(inputFolder+fileName)
                f.readline() #header
                #reads all lines in the file, stores it in dailyPos, for later use
                for line in f:
                    splitted = line.split()
                    currentDate = dateList.index(int(pred.asday(dt.datetime.strptime(splitted[1].split('T')[0],'%Y-%m-%d').date())))
                    dailyPos[currentDate,:] = map(float,splitted[2:5])
                print('File successfully read')
                f.close()
                
                dailyPos /= errors #weighting of the obs. The noisier the less weight.
                
                rowNb=0
                lenDailyPos = len(dailyPos)
                significanceStorage = np.zeros((len(dateList),len(windowSizes))) #dates, windows. Stores result of FTest for each day, each window, each dim
                ftestStorage = np.zeros((len(dateList),len(windowSizes)))
                jumpSizeStorage = np.zeros((len(dateList),len(windowSizes),2)) #2 extra dim for plani/alti
                #jumpMagnitude = np.zeros((len(dateList),10)) #stores magnitude of jump everyday, every window
                for row in dailyPos:
                    if not (rowNb+1)%200:
                    #if True:
                        print('Computing day '+str(rowNb+1)+'/'+str(lenDailyPos))
                    if daysToCompute[rowNb]:   #if the day is to be computed
                        windowNb = 0
                        for N in windowSizes:
                            if N < rowNb < len(daysToCompute)-N:
                                minIndex = max(rowNb-N,0)
                                maxIndex = min(rowNb+N,lenDailyPos)
                                window = dailyPos[minIndex:maxIndex,:] #size-changing window
                                window_dates = np.array(dateList[minIndex:maxIndex])
                                lenWindow = len(window)
                                
                                hole0=int(np.floor(lenWindow*(1.0/2-halfLengthOfRamp)))
                                hole1=int(np.ceil(lenWindow*(1.0/2+halfLengthOfRamp)))
                                date0=window_dates[hole0]
                                date1=window_dates[hole1]
                                zeros = np.all(window[:,1:]==0,axis=1) #reject days without obs
                                
                                #window = np.ma.array(window,mask=np.vstack((zeros,zeros,zeros)).transpose()) #masking array to get only the days with obs
                                #window_dates = np.ma.array(window_dates,mask=zeros)
                                
                                window = np.delete(window,np.where(zeros==1),axis=0)
                                window_dates = np.delete(window_dates,np.where(zeros==1))
                                
                                if date1 != date0:
                                    ramp=np.fmax(0.0,np.fmin(1.0,(window_dates-date0).astype(float)/(date1-date0)))
                                else:
                                    ramp=[0]*len(window_dates[:hole0])+[1]*len(window_dates[hole0:])
                                    
                                params1=np.vstack((window_dates,np.ones(window_dates.shape)))
                                params=np.vstack((params1,ramp))
                                params1=params1.T
                                params=params.T
                                
                                #if np.ma.MaskedArray.count(firstHalf_dates)>N/2 and np.ma.MaskedArray.count(secondHalf_dates)>N/2:
                                if window_dates[:hole0].size > N/2 and window_dates[hole1:].size > N/2 :
                                    pT, residualsT, rankT, singular_valuesT = np.linalg.lstsq(params1,window) #least-square estimation for the full set of obs
                                    p1, residuals1, rank1, singular_values1 = np.linalg.lstsq(params,window)
                                    
                                    ssrT=np.sum(residualsT)
                                    ssr1=np.sum(residuals1)
                                    
                                    dofT=window.size-6
                                    dof1=window.size-9
                                    
                                    p1 = abs(p1[2]*errors)
                                    
                                    if ssr1!=0:
                                        fstat=((ssrT-ssr1)*dof1)/(3*ssr1)
                                        #significanceStorage[rowNb,windowNb] = scipy.stats.f.cdf(fstat,dof1+dof2,window.count()-dof1-dof2-1)
                                        #significanceStorage[rowNb,windowNb]=scipy.stats.f(3,dof1).cdf(fstat)
                                        significanceStorage[rowNb,windowNb]=scipy.stats.f.cdf(fstat,3,dof1)
                                        ftestStorage[rowNb,windowNb] = fstat
                                        jumpSizeStorage[rowNb,windowNb,0] = np.power(abs(p1[0])**3+abs(p1[1])**3,1.0/3)
                                        jumpSizeStorage[rowNb,windowNb,1] = abs(p1[2])
                                    else:
                                        significanceStorage[rowNb,windowNb] = np.inf
                                        ftestStorage[rowNb,windowNb] = np.inf
                                        jumpSizeStorage[rowNb,windowNb,:] = np.inf
                            windowNb+=1
                    rowNb+=1
                    
                #np.savetxt('/home/cdrouadaine/timeseries/timeseries/significanceStorage.txt', 100*np.mean(significanceStorage,axis=0))                    
                    
                    #np.savetxt(currentFolder+'/significance/'+staName+'_'+str(iterationNb)+'.txt',significanceStorage[:,:,1],delimiter='\t')
                #significanceAverage = np.mean(significanceStorage,axis=1)
                    
                #now going through the calculated data to find the most significant jumps
                jumpFound = False
                print('==========================')
                #print(np.where(significanceStorage==np.amax(significanceStorage))[0][0],np.amax(significanceStorage))
                print(np.where(significanceStorage==np.amax(significanceStorage))[0],np.amax(significanceStorage))
                '''
                print(np.where(significanceStorage==np.max(significanceStorage)))
                time.sleep(15)
                '''
                if np.amax(significanceStorage) < minimumSignificance:
                    isNotEnd = False
  
                while not jumpFound and isNotEnd:
                    #dataOnMostSignificant = np.where(significanceStorage > np.amax(0.98*significanceStorage))
                    dataOnMostSignificant = np.where(ftestStorage > np.amax(0.98*ftestStorage))
                    #indexOfMostSignificant = dataOnMostSignificant[0][0] #gets the day number
                    #windowOfMostSignificant = dataOnMostSignificant[1][0] #gets the number of the window size where max was found
                    
                    #IF SEVERAL DAYS WITH MAX SIGNIFICANCE
                    #First, choose 
                    windowOfMostSignificant = dataOnMostSignificant[1]
                    tmpmax = 0
                    for x in np.unique(windowOfMostSignificant):
                        numberOfOccurences = np.count_nonzero(windowOfMostSignificant == x)
                        if numberOfOccurences >= tmpmax:
                            tmpmax = numberOfOccurences
                            possibleWindow = x
                    windowOfMostSignificant = possibleWindow
                    
                    #possibleIndexes = np.where(significanceStorage[:,windowOfMostSignificant]==np.amax(significanceStorage[:,windowOfMostSignificant]))
                    possibleIndexes = np.where(ftestStorage[:,windowOfMostSignificant]==np.amax(ftestStorage[:,windowOfMostSignificant]))
                    maxOfFtestStorage = np.amax(ftestStorage[possibleIndexes,windowOfMostSignificant])
                    indexOfMostSignificant = np.where(ftestStorage[:,windowOfMostSignificant] == maxOfFtestStorage)[0][0]
                
                    dateOfMostSignificant = dateList[indexOfMostSignificant]
                    mySignificance = significanceStorage[indexOfMostSignificant,windowOfMostSignificant]
                    #significanceStorage[indexOfMostSignificant,:] = 0
                    significanceStorage[indexOfMostSignificant,:] = 0
                    ftestStorage[indexOfMostSignificant,:] = 0
                    
                    if mySignificance < minimumSignificance:
                        break
                            
                    print("Possible jump: "+str(pred.fromday(dateOfMostSignificant).date()))
                    '''if any(jumpSizeStorage[indexOfMostSignificant,windowOfMostSignificant,:] > 0.5):
                        print([abs(x-dateOfMostSignificant)>minimumDaysBetweenConsecutiveEarthquakes for x in overallJumps])
                        print([(dateOfMostSignificant > x[0] and dateOfMostSignificant < x[1]) for x in overallSlowSlips])
                        print(jumpSizeStorage[indexOfMostSignificant,windowOfMostSignificant,:])'''
                        
                    if all([abs(x-dateOfMostSignificant)>minimumDaysBetweenConsecutiveEarthquakes for x in overallJumps]) and any([jumpSizeStorage[indexOfMostSignificant,windowOfMostSignificant,w] > minimumJumpSizes[w] for w in [0,1]]): #if jump far enough from the jumps already found, add it, otherwise, just delete it and continue the search
                    #and not any([(dateOfMostSignificant > x[0] and dateOfMostSignificant < x[1]) for x in overallSlowSlips])
                    #if True:
                        daysToCompute[indexOfMostSignificant] = 0
                        
                        #fits to find the date
                        m2 = m.copy()
                        
                        #np.savetxt('/home/cdrouadaine/timeseries/timeseries/m.txt',m.calc(m.dates))
                        #m.save('/home/cdrouadaine/m.xml')
                        m2.loadTimeSeries(filename=currentFolder+subfolder+staName+'_igs08_xyz_smoothed.dat')
                        #np.savetxt('/home/cdrouadaine/timeseries/timeseries/m2_obs.txt',m2.getObs()[1])
                        myOffset = pred.equipment_offset(m2,pred.fromday(dateOfMostSignificant)) #add jump to model
                        m2.addComponent(myOffset)
                        ok,mesg = m2.fitAllLinear(myOffset)
                        if not ok:
                            print mesg
                            
                        
                        '''m2.save('/home/cdrouadaine/m2.xml')
                        myOffset.parameters[0].setFixed(False)
                        #print myOffset
   
                        np.savetxt('/home/cdrouadaine/timeseries/timeseries/m2bis.txt',m2.calc(m2.dates))
                        m2.save('/home/cdrouadaine/m2bis.xml')
                        ok,mesg=m2.fitAllLinear(myOffset)
                        if not ok:
                            print "2",mesg
                        m2.save('/home/cdrouadaine/m2.xml')
                        print(currentFolder+subfolder+staName+'_igs08_xyz_smoothed.dat')
                        np.savetxt('/home/cdrouadaine/timeseries/timeseries/m2.txt',m2.calc(m2.dates))
                        sys.exit()
                        time.sleep(500)
                        myOffset.parameters[0].setFixed(True)'''
                        #offsetDate = pred.asday(myOffset.eventDate())
                
                        jumpValues = myOffset.paramValues()
                        print(jumpValues)
                        jumpMagnitude[0] = np.power(abs(jumpValues[1])**3+abs(jumpValues[2])**3,1.0/3)*1000
                        jumpMagnitude[1] = abs(jumpValues[3])*1000
                        print('Magnitude: +'+str(jumpMagnitude))
                        
                        #GET SIZE OF OFFSET HERE, AND CHECK AGAIN
                        if not any([jumpMagnitude[w] > minimumJumpSizes[w]/2.0 for w in [0,1]]):
                            break
                        
                        if not any([jumpMagnitude[w] > minimumJumpSizes[w] for w in [0,1]]):
                            continue
                
                        '''if offsetDate != dateOfMostSignificant:
                            mySignificance = np.amax(significanceStorage[offsetDate,:])
                            #dataOnMostSignificant = np.where(significanceStorage[offsetDate,:]==np.amax(significanceStorage[offsetDate,:]))
                            dataOnMostSignificant = np.where(ftestStorage[offsetDate,:]==np.amax(ftestStorage[offsetDate,:]))
                            indexOfMostSignificant = offsetDate #gets the day number
                            windowOfMostSignificant = dataOnMostSignificant[0][0] #gets the number of the window size where max was found
                            dimOfMostSignificant = dataOnMostSignificant[1][0] #gets dim
                            dateOfMostSignificant = dateList[indexOfMostSignificant]
                            significanceStorage[indexOfMostSignificant,:] = 0
                            ftestStorage[indexOfMostSignificant,:] = 0'''
                        
                        jumpFound = True
    
                        windowSizes[ windowOfMostSignificant ]      #KINDA WEIRD THAT WINDOW LENGTH IS FIXED
    
                        lenWindow = len(window)
                        jumpSize = np.ceil(lenWindow*(1.0/2+halfLengthOfRamp)) - np.floor(lenWindow*(1.0/2-halfLengthOfRamp))
                        
                        minIndex = max(indexOfMostSignificant-lenWindow/2,0)
                        maxIndex = min(indexOfMostSignificant+lenWindow/2,lenDailyPos)
                        window = dailyPos[minIndex:maxIndex,:] #size-changing window
                        window_dates_full = np.array(dateList[minIndex:maxIndex])
                        
                        zeros = np.all(window[:,1:]==0,axis=1) #reject days without obs
                        window = np.delete(window,np.where(zeros==1),axis=0)
                        window_dates = np.delete(window_dates_full,np.where(zeros==1))

                        #testOffest = pred.equipment_offset(m2,pred.fromday(dateOfMostSignificant))                        
                        m3 = m.copy()
                        m3.loadTimeSeries(filename=currentFolder+subfolder+staName+'_igs08_xyz_smoothed.dat')
                        
                        testSlowSlip = pred.slow_slip(m3,date=pred.fromday(dateOfMostSignificant),duration=jumpSize)
                        testSlowSlip.parameters[0].setFixed(True)
                        for paramNb in [1,2,3,4]:
                            testSlowSlip.parameters[paramNb].setFixed(False)
                        m3.addComponent(testSlowSlip)
                        m3.fitAllLinear(testSlowSlip)
                        testSlowSlip.parameters[1].setFixed(True)
                        
                        slowSlipLength = testSlowSlip.paramValues()[1]
                        if slowSlipLength < minimumDaysForASlowSlip: #offset
                            m = m2.copy()
                            m.loadTimeSeries(filename=currentFolder+subfolder+staName+'_igs08_xyz_smoothed.dat')
                            #ADD TO THE LISTS
                            print(u"\u2191 Taken as a jump")
                            print("Significance: "+str(np.amax(mySignificance)))
                            overallJumps.append(dateOfMostSignificant)
                            for x in range(max(indexOfMostSignificant-minimumDaysBetweenConsecutiveEarthquakes,0),min(indexOfMostSignificant+minimumDaysBetweenConsecutiveEarthquakes,len(daysToCompute))):
                                daysToCompute[x] = 0
                                
                        else: #slow slip
                            m = m3.copy()
                            m.loadTimeSeries(filename=currentFolder+subfolder+staName+'_igs08_xyz_smoothed.dat')
                            print(u"\u2191 Taken as a slow slip")
                            print("Significance: "+str(np.amax(mySignificance)))
                            overallSlowSlips.append([dateOfMostSignificant-slowSlipLength,dateOfMostSignificant+slowSlipLength])
                            for x in range(max(int(indexOfMostSignificant-slowSlipLength-minimumDaysBetweenConsecutiveEarthquakes),int(0)),min(int(indexOfMostSignificant+slowSlipLength+minimumDaysBetweenConsecutiveEarthquakes),int(len(daysToCompute)))):
                                daysToCompute[x] = 0

                    else:
                        significanceStorage[indexOfMostSignificant,:] = 0
                        daysToCompute[indexOfMostSignificant] = 0
    
                if not jumpFound:
                    isNotEnd = False
                
                if isNotEnd:
                    #now we have found all most significant jumps for this iteration
                    overallJumps.sort()
                    print(["%i" %(x-dateList[0]) for x in overallJumps])
                    print(overallSlowSlips)
                    print('Fitting')
                    #fit that, compute residuals, and do it again
                    m.save(currentFolder+'/xmlWithJPLOffsets/'+staName+'.xml')
                    print('Computing residuals')
                    dates,enu,useobs = m.getObs()
                    enu_calc = m.calc(m.dates)
                    residuals = 1000*(enu-enu_calc)
                    use.toText(dates,residuals,inputFolder+fileName,staName)
                   
                iterationNb += 1
            print('CHANGING LENGTH OF RAMP')

            
        jumpsFile = open(currentFolder+'/JPL_offsets/'+staName+'_jumps.dat','a') #APPENDS THE NEW DATE INTO THE FILE
        for jump in overallJumps:
            jumpsFile.write(dt.datetime.strftime(pred.fromday(jump)-twelveHours,'%Y-%m-%dT%H:%M:%S')+'\n') #add jump to file
        jumpsFile.close()
        
        slowSlipsFile = open(currentFolder+'/JPL_offsets/'+staName+'_slowslips.dat','a') #APPENDS THE NEW DATE INTO THE FILE
        for slowSlip in overallSlowSlips:
            slowSlipsFile.write(dt.datetime.strftime(pred.fromday(slowSlip[0])-twelveHours,'%Y-%m-%dT%H:%M:%S')+'\t'+dt.datetime.strftime(pred.fromday(slowSlip[1])-twelveHours,'%Y-%m-%dT%H:%M:%S')+'\n') #add jump to file
        slowSlipsFile.close()
        
        print('FINAL COUNT: '+str(len(overallJumps)+len(overallSlowSlips))+' EVENTS ('+str(len(overallJumps))+' JUMPS, '+str(len(overallSlowSlips))+' SLOW SLIPS) FOR '+staName)
        fileNb+=1

'''def FTest(valueT,value1,value2,order,degreesOfFreedom):
    return ((valueT-(value1+value2))*degreesOfFreedom)/(order*(value1+value2))

folder = "/home/cdrouadaine/timeseries/timeseries"
t = time.time()
JPL_STP1(folder)
print(time.time()-t)
'''
