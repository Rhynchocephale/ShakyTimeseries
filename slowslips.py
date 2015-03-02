import numpy as np
import os
import datetime as dt
import stn_pred_model as pred



def detectSlowslip(folder)

    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    inputFolder = folder+'/smallOffsetsRemoved/'
    for f in os.listdir(inputFolder):
        if os.path.isfile(os.path.join(inputFolder,f)) and f.endswith('.dat'):
            fileList.append(f.split('_')[0])
    fileList.sort
    
    subfolder ='/'
    
    for staName in fileList:
        jumpsFile = open(inputFolder+staName+'_jumps.dat','r')
        residualsFile = open(folder+'/residuals/'+staName+'_residuals.dat','r')
        rawFile = open(folder+subfolder+staName+'_igs08_xyz.dat','r')
        
        residualsDate = dt.date.min
        residualsFile.readline() #header
        rawFile.readline() #header
        theEndIsNear = False
        while not theEndIsNear:
            jumpLine = jumpsFile.readline()
            if not jumpLine:
                theEndIsNear = True #last iteration before end, no jump date to end the interval
            jumpDate = dt.datetime.strptime(jumpLine,'%Y-%m-%d').date()
            residualsPos = np.array([])
            residualsDateList = np.array([])
            
            while theEndIsNear or residualsDate < jumpDate:
                residualsLine = residualsFile.readline().split()
                if not residualsLine:
                    break                
                rawLine = rawFile.readline().split()
                
                residualsDate =  dt.datetime.strptime(residualsLine[1].split('T')[0],'%Y-%m-%d').date()
                residualsDateList.hstack(pred.asday(residualsDate))
                
                residualsPos.hstack(map(float,residualsLine[2:]))
                rawPos.hstack(map(float,rawLine[2:]))
            
            avgSquaredResiduals = np.mean(residualsPos**2,axis=1)
            
            if avgSquaredResiduals > whatever: #if too different than model
                p, V = np.polyfit(residualsDateList,rawPos,1)
                velocity = p[1]
        
        xyz = use.getxyz(currentFolder+subfolder+staName+'_igs08_xyz.dat')
        m = pred.model(station=staName, xyz=xyz)
        m.loadTimeSeries(filename=currentFolder+subfolder+staName+'_igs08_xyz.dat')
        
        
        
        residualsFile.readline()
