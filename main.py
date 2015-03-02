import remAvgRes_withoutSplines
from usePredModel import *
from coordConv import *
#import rpy2.robjects as robjects
import os
import shutil
import time
import JPL_STP1 as JPL
from drawStuff import *

def main(folder):
    t = time.time()
    all_xyz2enu(folder)

    #------take enu, give enu. All of them--------  
    remAvgRes_withoutSplines.remAvgRes(folder)
    all_enu2xyz(folder)
    
    getResiduals(folder)
    return
    JPL.JPL_STP1(folder,hole=True)
    #-------take xyz, give xyz------------
    putOffsets(folder,False)
    #getSlowSlips(folder)
    #all_trans2xyz(folder)
    drawStuff()
    print(time.time()-t)
    
main("/home/cdrouadaine/timeseries/timeseries")
'''
folder = '/home/cdrouadaine/timeseries/timeseries-complete'
makeAllTrusted(folder,folder+'/stationsToCompute.txt')
'''
