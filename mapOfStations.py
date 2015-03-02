import numpy as np
import matplotlib.pyplot as plt
import os

def mapStations(folder):
    #BROWSE ALL FILES IN FOLDER AND SUBSTRACT STUFF
    #only takes the files, not the folders, and checks that the file is a .dat
    fileList = []
    for f in os.listdir(folder):
        if os.path.isfile(os.path.join(folder,f)) and f.endswith('.dat'):
            fileList.append(f)
    fileList.sort()
    
    listOfPos = np.empty((1,3))
    listOfNames=[]
    for fileName in fileList:
        f = open(folder+'/'+fileName,'r')
        f.readline() #header
        xyz = np.array(map(float,f.readline().split('\t')[2:5]))
        staName = np.array([fileName.split('_')[0]])
        listOfPos = np.vstack((listOfPos,xyz))
        listOfNames.append(staName)
        f.close()
    
    listOfPos = np.delete(listOfPos, (0), axis=0)
    #print(listOfPos)
    #plt.scatter(listOfPos[:,1],data[:,2],marker='o')
    plt.plot(listOfPos[:,0],listOfPos[:,1],'+b')
    for label, x, y in zip(listOfNames, listOfPos[:, 0], listOfPos[:, 1]):
        plt.annotate(
        label, 
        xy = (x, y), xytext = (-20, 20),
        textcoords = 'offset points', ha = 'right', va = 'bottom',
        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
    plt.savefig('/home/cdrouadaine/timeseries/timeseries/'+'gotcha.png',bbox_inches='tight')
    plt.show()
    
mapStations("/home/cdrouadaine/timeseries/timeseries")
