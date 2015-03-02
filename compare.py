import numpy as np
import os

#compares the differences between original & smoothed
def compareFiles(folder1,folder2):
    
    currentFolder = os.path.dirname(os.path.abspath(__file__))

    if not os.path.exists(currentFolder+'/compare/'+folder1.split('/')[-2]+'&'+folder2.split('/')[-2]+'/'): #creates the /transverseMercator directory if it does not already exists
        os.makedirs(currentFolder+'/compare/'+folder1.split('/')[-2]+'&'+folder2.split('/')[-2]+'/')

    fileList = []
    for f in os.listdir(folder1): #get list of files in current folder
        if os.path.isfile(os.path.join(folder1,f)) and f.endswith('.dat'):
            fileList.append(f.split('_')[0])
    
    mm = 10 #if correction (in mm) is bigger than this, notes it
    staCount = 1
    lenFileList = len(fileList)
    for sta in fileList:
        print(sta+' '+'{:.0f}/{:.0f}'.format(staCount,lenFileList))
        
        f = open(folder1 + sta + '_igs08_xyz.dat','r')
        g = open(folder2 + sta + '_igs08_xyz_smoothed.dat','r')
        h = open(currentFolder+'/compare/'+folder1.split('/')[-2]+'&'+folder2.split('/')[-2]+'/' + sta + '_compare.dat','w')
        
        f.readline()
        g.readline()
        h.write('header\n')
        u=0
        v=0
        w=0
        while True:
            fline=f.readline() 
            gline=g.readline()
            if not fline: break
            
            fline = fline.split()
            gline = gline.split()
            fxyz = map(float,fline[2:5])
            gxyz = map(float,gline[2:5])
            
            hxyz = (np.array(fxyz)-np.array(gxyz))*1000
            #print("Huge numbers here! %f, %f, %f") %(hxyz[0],hxyz[1],hxyz[2])
            u+=any(hxyz > mm)
            v+=1
            w+=(0,1) [np.count_nonzero(hxyz > mm)>1]
            
            h.write("{:.1f} \t{:.1f} \t{:.1f}\n".format(hxyz[0],hxyz[1],hxyz[2]))

        print "Corr > %imm: %i/%i days: %f%% of total, %i double: %f%% of days, %f%% of total\n" %(mm,u,v,(float(u)/v)*100,w,(float(w)/u)*100,(float(w)/v)*100)

        f.close()
        g.close()
        h.close()

        staCount+=1
    
compareFiles('/home/cdrouadaine/timeseries/timeseries/','/home/cdrouadaine/timeseries/timeseries/smoothed/')

