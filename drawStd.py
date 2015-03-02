from remAvgRes import *
import numpy as np
import os
import os.path
import pylab as pl
import time
import coordConv


def drawStd():
    currentFolder = os.path.dirname(os.path.abspath(__file__))
    makeAllTrusted(currentFolder + '/trustEveryone.txt')
    enu = remAvgRes(currentFolder + '/trustEveryone.txt')
    enu = np.std(enu,0)
    #y=range(0,enu.shape[0])
     
    '''
    pl.title('Standard deviation')
        
    pl.subplot(3,1, 1)  
    pl.plot(enu[:,0],',')
    pl.ylabel('e')
    
    pl.subplot(3, 1, 2)
    pl.plot(enu[:,1],',')
    pl.ylabel('n')
    
    pl.subplot(3, 1, 3)
    pl.plot(enu[:,2],',')
    pl.ylabel('u')
    
    pl.subplot(4, 1, 4)
    pl.plot(enu_corr[:,3],'r,')
    pl.ylabel('Sta nb')
    
    pl.show()
    '''
    
t = time.time()
#coordConv.all_xyz2trans('/home/cdrouadaine/timeseries/timeseries/')
drawStd()
#coordConv.all_trans2xyz('/home/cdrouadaine/timeseries/timeseries/smoothed/')
elapsed = time.time() - t
print(elapsed)

