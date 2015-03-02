import os
import os.path
import numpy as np
import numpy.ma as ma
import pylab as pl
import scipy
import scipy.stats

def lineCount(filePath):
    count = 0
    thefile = open(filePath, 'rb')   
    while True:
        buffer = thefile.read(8192*1024)
        if not buffer: break
        count += buffer.count('\n')
    thefile.close()
    return count


def corrEverything():
    currentFolder = os.path.dirname(os.path.abspath(__file__))
    
    #list of files in folder
    fileList = []
    for f in os.listdir(currentFolder + '/timeseries/enu'):
        if os.path.isfile(os.path.join(currentFolder + '/timeseries/enu',f)) and f.endswith('.enu'):
            fileList.append(f)
            
    deltaT = lineCount(currentFolder + '/timeseries/enu/'+fileList[0])

    g = open('corr_res.corr','w')
    g.write('sta1 sta2 corrE      corrN      corrU      days distance_in_km\n') #header
    k=1
    length = len(fileList)
    total = scipy.misc.comb(length,2)
    for i in range(0,length):
        for j in range(i+1,length):
            print "%i/%i" %(k,total+1)
            corr = correlate(currentFolder+'/timeseries/enu/'+fileList[i],currentFolder+'/timeseries/enu/'+fileList[j],deltaT)
            g.write(fileList[i].split('.')[0]+' '+fileList[j].split('.')[0]+" {:.8f} {:.8f} {:.8f} {:.0f} {:.2f}\n".format(corr[0],corr[1],corr[2],corr[3],corr[4]))
            k+=1
    g.close()

    printCorr('corr_res.corr')

def correlate(staFile1, staFile2,deltaT):
        
    f1 = open(staFile1,'r')
    f2 = open(staFile2,'r')
    M = np.zeros(shape=(deltaT,3,2))
    Mmask = np.zeros(shape=(deltaT,3,2))
    
    line1 = np.array(map(float,f1.readline().split()[1:4]))
    line2 = np.array(map(float,f2.readline().split()[1:4]))
    dist = np.linalg.norm(line1-line2) / 1000
    
    lineNb = 0
    while True:
        line1 = f1.readline().split()
        line2 = f2.readline().split()
        if not line1: break

        M[lineNb,:,0] = map(float,line1[1:4])
        M[lineNb,:,1] = map(float,line2[1:4])
            
        lineNb += 1
    
    #matrix of distances between the corrections of the two stations to correlate
    M2 = np.zeros(deltaT)
    lineNb = 0
    '''CAREFUL, MULTIDIMENSIONAL THINGIE'''
    for row in M:
        if not any(M[lineNb,:,0]) or not any(M[lineNb,:,1]): #if either one full of zeros
            Mmask[lineNb,:,:].fill(1)
        else:
            M2[lineNb] = [abs(np.linalg.norm(M[lineNb,:,0])-np.linalg.norm(M[lineNb,:,1]))][0]
        lineNb += 1
    
    M2 = np.argsort(M2)[::-1][0:(deltaT-np.count_nonzero(Mmask[:,0,0]))*5/100]  #takes the bigger 5% of values
    #removes corresponding rows from M
    for row in M2:
        Mmask[row,:,:] = 1
    
    M = ma.array(M,mask=Mmask)
    
    f1.close()
    f2.close()

    corr = [0,0,0,Mmask.shape[0]-np.count_nonzero(Mmask[:,0,0]),dist] #corr in e, corr in n, corr in u, number of correlated days, distance between stations
    for i in range(0,3):
        #pl.plot(M[:,i,0],M[:,i,1],'+')
        #pl.show()
        corr[i] = scipy.stats.pearsonr(M[:,i,0],M[:,i,1])[0]
    
    return corr
    
def printCorr( filePath ):
    
    fileLen = lineCount(filePath)-1
    f = open(filePath,'r')
    f.readline() #skips header
    M = np.zeros(shape=(fileLen,5))
    lineNb = 0
    for line in f:
        line = line.split()
        corr = map(float,line[2:5])
        dist = float(line[6])
        days = float(line[5])
        
        M[lineNb,:] = [(corr[0]+corr[1])/2,(corr[0]+corr[1]+corr[2])/3,corr[2],dist,days]
    
        lineNb+=1
    
    z = M[:,4]

    pl.title('Evolution of correlation with distance')
    '''    
    pl.subplot(3,1,1)
    pl.scatter(M[:,3],M[:,0],c=z,s=1,marker=',',linewidth=0)
    pl.xlabel('dist (km)')
    pl.ylabel('en')
    pl.xlim([min(M[:,3]),max(M[:,3])+100])
    pl.ylim([-1,1])
    '''
    #pl.subplot(3,1,2)
    pl.scatter(M[:,3],M[:,1],c=z,s=2,marker=',',linewidth=0)
    pl.xlabel('dist (km)')
    pl.ylabel('enu')
    pl.xlim([min(M[:,3]),max(M[:,3])+100])
    pl.ylim([0,1])
    ax = pl.gca()
    ax.patch.set_facecolor('0.5')
    '''
    pl.subplot(3,1,3)
    pl.scatter(M[:,3],M[:,2],c=z,s=1,marker=',',linewidth=0)
    pl.xlabel('dist (km)')
    pl.ylabel('u')
    pl.xlim([min(M[:,3]),max(M[:,3])])
    pl.ylim([-1,1])
    '''
    pl.show() 

#corrEverything()
printCorr('/home/cdrouadaine/timeseries/corr_res.corr')
