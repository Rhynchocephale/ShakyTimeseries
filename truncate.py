import os
import numpy as np
import datetime as dt

def truncate(folder,beginning,end):
    
    #BROWSE ALL FILES IN FOLDER
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    if os.path.exists(folder+'/'+'truncated'+str(duration.days)): #clears the directory if it already exists
        shutil.rmtree(folder+'/'+'truncated'+str(duration.days))
    os.makedirs(folder+'/'+'truncated'+str(duration.days))
    
    i=1
    for fileName in fileList:
        staName = fileName.split('_')
        print('Truncating '+staName+' '+str(i)+'/'+str(len(fileName)))
        f = open(folder+'/'+fileName,'r')
        g = open(folder+'/'+'truncated'+str(duration.days)+'/'+fileName)

        line = f.readline() #header
        g.write(line)
        currentDate = dt.date.min
        
        while currentDate < beginning
            line = f.readline().split()
            if not line:
                break
            currentDate = dt.datetime.strptime(line[1].split('T')[0],'%Y-%m-%d').date()
        
        for line in f
            fline = f.readline()
            line = fline.split()
            currentDate = dt.datetime.strptime(line[1].split('T')[0],'%Y-%m-%d').date()
            if currentDate > end:
                break
            g.write(line)

        f.close()
        g.close()
        i+=1
