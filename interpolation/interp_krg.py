import numpy as np
from variogramme import *

## Cette fonction permet de calculer les valeurs (ZO) aux points 
## d'interet par une interpolation par krigeage.



def interp_krg(X,Y,Zinterp,Zalti):

    [C0,a0,vario,gam0]=variogramme(X,Y,Zalti,Zinterp)              ## appel de la fonction variogramme

    tail=len(X)
###   gam0=C0*(1-exp(-(h)/a0));

    for j in range(0,tail):
        B[:,j]=C0*(1-np.exp(-(np.sqrt(np.power(X-X[j],2)+np.power(Y-Y[j],2)+np.power(Zalti-Zalti[j],2)/a0)))) ## remplissage de l'estimateur 

    C=np.ones((1,tail))                                ## initialisation vecteur ligne
    A=np.bmat([B,C.transpose()],[C,0])               ## concatenation des matrices B C', C et de l'element nul
    #O=AP
    invers=np.linalg.pinv(A.transpose().dot(A)).dot(A.transpose())            
    Z0=np.zeros((tail,));

    for i in range(0,tail):
        O=C0*(1-np.exp(-np.sqrt(np.power(X[i]-X,2)+ np.power(Y[i]-Y,2) +np.power(Zalti[i]-Z,2))/a0))     
        O=O.append(1)
        P=invers*O                                              ## obtention du vecteur des parametres (lambda)
        P=P[:len(P)-1]
        z=P.dot(Zinterp)
        Z0[i]=sum(z)

    return [Z0,vario,gam0]
