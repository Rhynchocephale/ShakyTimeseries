import numpy as np

##distance : Calcul de distances euclidiennes
##---ENTREES---
##xa : vecteur des coordonnees de X du premier jeu de points
##xb : vecteur des coordonnees de X du second jeu de points
##ya : vecteur des coordonnees de Y du premier jeu de points
##yb : vecteur des coordonnees de Y du second jeu de points
##---SORTIES---
##dist : matrice des distances euclidiennes entre les points des vecteurs
def distance(xa,xb,ya,yb):
    if isinstance(xb, np.float64):
        dist = [np.sqrt(np.power((xb - xa),2) + np.power((yb - ya),2))]
    else:
        n=len(xa)
        m=len(xb)
        dist=np.zeros(shape=(n,m))
        for i in range(0,m):
            dist[:,i] = np.sqrt(np.power((xb[i]- xa),2) + np.power((yb[i] - ya),2))
    
    return dist
